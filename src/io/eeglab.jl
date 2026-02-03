"""
EEGLAB .set file import functionality.

Provides functions to read EEGLAB .set files (MATLAB format) and convert them
to EegFun datatypes.
"""


"""
    load_eeglab(filepath::String)  → EpochData

Load EEGLAB .set file and convert to EegFun EpochData.

Currently supports epoched data only. Returns an `EpochData` object containing
all epochs from the .set file.

# Arguments
- `filepath::String`: Path to .set file

# Returns
- `EpochData`: Converted epoched data

# Examples
```julia
using EegFun

# Load epochs from EEGLAB
epochs = load_eeglab("participant1.set")

# Use with EegFun functions
erp = average_epochs(epochs)
plot_erp(erp)
```

# Notes
- Requires MAT.jl for reading MATLAB files
- Channel locations are extracted if available
- Currently does not import: ICA data, event structures, rejection info
"""
function load_eeglab(filepath::String)

    @info "Loading EEGLAB .set file: $filepath"

    # Read .set file using MAT.jl
    eeg = matread(filepath)

    # Detect data type
    n_trials = Int(eeg["trials"])

    if n_trials == 1
        error("Continuous data not yet supported. Only epoched data can be imported.")
    end

    # Convert to EpochData
    return _eeglab_to_epochdata(eeg, filepath)
end


"""
    _eeglab_to_epochdata(eeg::Dict, filepath::String) → EpochData

Convert EEGLAB structure to EegFun EpochData.

Internal function called by `load_eeglab`.
"""
function _eeglab_to_epochdata(eeg::Dict, filepath::String)

    # Extract basic dimensions
    n_channels = Int(eeg["nbchan"])
    n_timepoints = Int(eeg["pnts"])
    n_trials = Int(eeg["trials"])
    sample_rate = Int(round(eeg["srate"]))

    @info "Data: $n_channels channels × $n_timepoints timepoints × $n_trials trials"
    @info "Sample rate: $sample_rate Hz"

    # Extract data array [channels × timepoints × trials]
    data_array = eeg["data"]

    # Extract channel names
    ch_names = _extract_channel_names(eeg["chanlocs"])

    # Create time vector
    if haskey(eeg, "times") && !isempty(eeg["times"])
        times = vec(eeg["times"]) ./ 1000.0  # Convert ms to seconds
    else
        # Calculate from xmin and sample rate
        xmin = Float64(eeg["xmin"])
        dt = 1.0 / sample_rate
        times = xmin .+ (0:(n_timepoints-1)) .* dt
    end

    # Extract triggers from epochs (returns both integer codes and string info)
    trigger_codes, trigger_info = _extract_triggers(eeg["epoch"], n_trials)

    # Convert 3D array to Vector{DataFrame}
    epoch_dataframes = Vector{DataFrame}(undef, n_trials)

    for trial_idx = 1:n_trials
        # Extract trial data [channels × timepoints]
        trial_data = data_array[:, :, trial_idx]'  # Transpose to [timepoints × channels]

        # Create DataFrame with channel data
        df = DataFrame(trial_data, ch_names)

        # Insert time column at the beginning
        insertcols!(df, 1, :time => times)

        # Add trigger columns (mark trigger only at time=0)
        # Find index closest to time=0
        zero_idx = argmin(abs.(times))

        # Create trigger vectors (0 everywhere except at t=0)
        trigger_vec = zeros(Int, n_timepoints)
        trigger_vec[zero_idx] = trigger_codes[trial_idx]

        trigger_info_vec = fill("", n_timepoints)
        trigger_info_vec[zero_idx] = trigger_info[trial_idx]

        insertcols!(df, 2, :triggers => trigger_vec)
        insertcols!(df, 3, :trigger_info => trigger_info_vec)

        epoch_dataframes[trial_idx] = df
    end

    # Parse layout (channel locations)
    layout = _parse_channel_locations(eeg["chanlocs"], ch_names)

    # Extract condition name if available
    condition_name = get(eeg, "setname", basename(filepath))

    # Create EpochData
    epoch_data = EpochData(
        filepath,          # file
        1,                 # condition
        condition_name,    # condition_name
        epoch_dataframes,  # data
        layout,            # layout
        sample_rate,       # sample_rate
        AnalysisInfo(),    # I guess the user would want to add their own analysis_info
    )

    @info "Successfully converted to EpochData with $n_trials epochs"

    return epoch_data
end


"""
    _extract_channel_names(chanlocs) → Vector{Symbol}

Extract channel names from EEGLAB chanlocs structure.
"""
function _extract_channel_names(chanlocs)
    if isempty(chanlocs)
        error("No channel location information found in .set file")
    end

    # Extract labels
    labels = chanlocs["labels"]

    # Convert to Symbol vector
    if labels isa AbstractVector
        return Symbol.(labels)
    elseif labels isa AbstractMatrix
        # Sometimes it's stored as a matrix
        return Symbol.(vec(labels))
    else
        error("Unexpected chanlocs.labels format")
    end
end


"""
    _extract_triggers(epochs, n_trials::Int) → Tuple{Vector{Int}, Vector{String}}

Extract trigger codes from EEGLAB epoch structure.

Maps unique trigger strings to sequential integers (1, 2, 3, ...) and returns both
the integer codes and original strings for storage in :trigger and :trigger_info columns.
"""
function _extract_triggers(epochs, n_trials::Int)
    # First pass: collect all unique trigger strings
    trigger_strings = Vector{String}(undef, n_trials)
    unique_values = Set{String}()

    for i = 1:n_trials
        epoch = epochs[i]

        # Try to extract event type
        if haskey(epoch, "eventtype") && !isempty(epoch["eventtype"])
            event_types = epoch["eventtype"]

            # Take first event type if multiple, convert to String
            if event_types isa AbstractArray
                trigger_str = string(first(event_types))
            else
                trigger_str = string(event_types)
            end

            trigger_strings[i] = trigger_str
            push!(unique_values, trigger_str)
        else
            # No event type - use empty string
            trigger_strings[i] = ""
            push!(unique_values, "")
        end
    end

    # Create mapping from strings to sequential integers
    string_to_int = Dict{String,Int}()
    for (idx, value) in enumerate(sort(collect(unique_values)))
        string_to_int[value] = idx
    end

    # Second pass: create integer trigger codes
    trigger_codes = [string_to_int[s] for s in trigger_strings]

    n_unique = length(unique_values)
    @info "Extracted $n_unique unique trigger types: $(sort(collect(unique_values)))"

    return trigger_codes, trigger_strings
end


"""
    _parse_channel_locations(chanlocs, ch_names::Vector{Symbol}) → Layout

Parse EEGLAB channel locations into EegFun Layout.

Converts 3D Cartesian coordinates (X, Y, Z) to spherical (inc, azi) if available,
otherwise creates a default layout.
"""
function _parse_channel_locations(chanlocs, ch_names::Vector{Symbol})
    # Try to extract 3D coordinates
    has_coords = haskey(chanlocs, "X") && haskey(chanlocs, "Y") && haskey(chanlocs, "Z")

    if has_coords && !isempty(chanlocs["X"])
        X = vec(chanlocs["X"])
        Y = vec(chanlocs["Y"])
        Z = vec(chanlocs["Z"])

        # Check for any missing coordinates (NaN or empty)
        valid_coords = .!isnan.(X) .& .!isnan.(Y) .& .!isnan.(Z)

        if all(valid_coords)
            # All coordinates valid - convert Cartesian to spherical
            inc, azi = _cartesian_to_spherical(X, Y, Z)
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing)
        else
            # Some missing - create default layout
            @warn "Some channels missing 3D coordinates, creating default layout"
            return Layout(_create_default_layout(ch_names), nothing, nothing)
        end
    else
        # No coordinates - create default layout
        @warn "No 3D coordinates found in .set file, creating default layout"
        return Layout(_create_default_layout(ch_names), nothing, nothing)
    end
end


"""
    _cartesian_to_spherical(X, Y, Z) → (inc, azi)

Convert 3D Cartesian coordinates to spherical coordinates (inclination, azimuth).

EEGLAB uses: X (front-back), Y (left-right), Z (inferior-superior)
Converts to: inc (0-90°), azi (0-360°)
"""
function _cartesian_to_spherical(X::Vector, Y::Vector, Z::Vector)
    n = length(X)
    inc = zeros(Float64, n)
    azi = zeros(Float64, n)

    for i = 1:n
        x, y, z = X[i], Y[i], Z[i]

        # Calculate radius
        r = sqrt(x^2 + y^2 + z^2)

        if r == 0
            # Handle origin case
            inc[i] = 0.0
            azi[i] = 0.0
        else
            # Inclination (angle from vertical/zenith, in degrees)
            # 0° = top (Cz), 90° = horizontal
            inc[i] = acosd(z / r)

            # Azimuth (angle in horizontal plane, in degrees)
            # 0° = front (Fpz), 90° = right, 180° = back, 270° = left
            # atan2(y, x) gives angle in radians, convert to degrees
            azi_rad = atan(y, x)
            azi[i] = rad2deg(azi_rad)

            # Rotate by 90 degrees to align EEGLAB with EegFun coordinate system
            azi[i] += 90.0

            # Ensure azimuth is in [0, 360)
            if azi[i] < 0
                azi[i] += 360.0
            elseif azi[i] >= 360.0
                azi[i] -= 360.0
            end
        end
    end

    return inc, azi
end


"""
    _create_default_layout(ch_names::Vector{Symbol}) → DataFrame

Create a default layout when no electrode coordinates are available.
Places all electrodes at the center.
"""
function _create_default_layout(ch_names::Vector{Symbol})
    n = length(ch_names)
    return DataFrame(label = ch_names, inc = zeros(Float64, n), azi = zeros(Float64, n))
end
