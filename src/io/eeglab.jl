"""
EEGLAB .set file import functionality.

Provides functions to read EEGLAB .set files (MATLAB format) and convert them
to EegFun datatypes.

TODO: This can be considered work in progress. I am not familiar with eeglab *.set/*.fdt files, so 
the code here is a bit of a guesswork. Please report any issues or missing features.
But it does seem to work with the two example datasets I found in eeglab/sample_data.
"""


"""
    read_eeglab(filepath::String; preserve_radial_distance::Bool = true)

Load EEGLAB .set file and return EEG data. If ICA decomposition is present,
also returns ICA data as a tuple (eeg_data, ica_data).

# Arguments
- `filepath::String`: Path to .set file
- `preserve_radial_distance::Bool`: If true, preserves anatomical distances for electrodes 
  with incidence > 90° (default: true)

# Returns
- If ICA present: `(eeg_data, ica_data::InfoIca)`
- If no ICA: `eeg_data` (EpochData, ContinuousData, or ErpData)

# Examples
```julia
# Without ICA
eeg = read_eeglab("data.set")

# With ICA (automatically detected)
eeg, ica = read_eeglab("data_with_ica.set")
```
"""
function read_eeglab(filepath::String; preserve_radial_distance::Bool = true)

    @info "Loading EEGLAB .set file: $filepath"

    mat_data = matread(filepath)

    # Handle nested EEGLAB format (data inside ["EEG"]) or direct format
    eeg = haskey(mat_data, "EEG") ? mat_data["EEG"] : mat_data

    # Detect data type and load EEG data
    n_trials = Int(eeg["trials"])
    xmin = get(eeg, "xmin", 0.0)

    eeg_data = if n_trials == 1
        # Single trial: could be continuous or averaged ERP
        # Use xmin heuristic: xmin < 0 suggests ERP (has pre-stimulus baseline)
        if xmin < 0
            @info "Detected averaged ERP data (xmin=$xmin < 0)"
            _eeglab_to_erpdata(eeg, filepath, preserve_radial_distance)
        else
            @info "Detected continuous data (xmin=$xmin >= 0)"
            _eeglab_to_continuousdata(eeg, filepath, preserve_radial_distance)
        end
    else
        # Multiple trials: epoched data
        _eeglab_to_epochdata(eeg, filepath, preserve_radial_distance)
    end

    # Check for ICA decomposition - keys always exist but matrices may be empty (0,0)
    has_ica = haskey(eeg, "icaweights") && !isempty(eeg["icaweights"])

    if has_ica
        @info "ICA decomposition detected, loading ICA data"
        ch_names = _extract_channel_names(eeg["chanlocs"])
        layout = _parse_channel_locations(eeg["chanlocs"], ch_names)
        polar_to_cartesian_xy!(layout, preserve_radial_distance = preserve_radial_distance)
        ica_data = _extract_ica_info(eeg, filepath, ch_names, layout)
        return eeg_data, ica_data
    else
        return eeg_data
    end
end


"""
Helper to extract common components used by all EEGLAB converters.
Returns: (n_channels, n_timepoints, sample_rate, ch_names, data, times, layout)
"""
function _load_common_eeglab_components(eeg::Dict, filepath::String, preserve_radial_distance::Bool)
    n_channels = Int(eeg["nbchan"])
    n_timepoints = Int(eeg["pnts"])
    sample_rate = Int(eeg["srate"])

    @info "Data: $n_channels channels × $n_timepoints timepoints" *
          (haskey(eeg, "trials") && eeg["trials"] > 1 ? " × $(Int(eeg["trials"])) trials" : "")
    @info "Sample rate: $sample_rate Hz"

    ch_names = _extract_channel_names(eeg["chanlocs"])
    data = _read_eeglab_data(eeg, filepath)

    xmin = get(eeg, "xmin", 0.0)
    times = if haskey(eeg, "times")
        vec(eeg["times"]) ./ 1000.0
    else
        dt = 1.0 / sample_rate
        xmin .+ (0:(n_timepoints-1)) .* dt
    end

    layout = _parse_channel_locations(eeg["chanlocs"], ch_names)
    polar_to_cartesian_xy!(layout, preserve_radial_distance = preserve_radial_distance)

    return n_timepoints, sample_rate, ch_names, data, times, layout
end


"""
    _eeglab_to_epochdata(eeg::Dict, filepath::String, preserve_radial_distance::Bool) → EpochData

Convert EEGLAB structure to EegFun EpochData.

Internal function called by `load_eeglab`.
"""
function _eeglab_to_epochdata(eeg::Dict, filepath::String, preserve_radial_distance::Bool)
    # Load common components
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    n_trials = Int(eeg["trials"])

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

    # Extract condition name if available
    condition_name = get(eeg, "setname", basename(filepath))

    # Create AnalysisInfo (EEGLAB doesn't have reliable filter/reference metadata)
    analysis_info = AnalysisInfo()

    # Create EpochData
    epoch_data = EpochData(
        filepath,          # file
        1,                 # condition
        condition_name,    # condition_name
        epoch_dataframes,  # data
        layout,            # layout
        sample_rate,       # sample_rate
        analysis_info,     # analysis_info
    )

    @info "Successfully converted to EpochData with $n_trials epochs"

    return epoch_data
end


function _eeglab_to_continuousdata(eeg::Dict, filepath::String, preserve_radial_distance::Bool)
    # Load common components
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    # Transpose to [timepoints × channels]
    data_matrix = data_array'

    # Extract triggers/events if present
    trigger_vec = zeros(Int, n_timepoints)
    trigger_info_vec = fill("", n_timepoints)

    if haskey(eeg, "event") && !isempty(eeg["event"])
        events = eeg["event"]

        # EEGLAB stores events as struct-of-arrays: .names and .values
        type_idx = findfirst(==("type"), events.names)
        lat_idx = findfirst(==("latency"), events.names)
        pos_idx = findfirst(==("position"), events.names)

        if lat_idx !== nothing && type_idx !== nothing
            latencies = events.values[lat_idx]
            types = events.values[type_idx]
            positions = pos_idx !== nothing ? events.values[pos_idx] : nothing

            @info "Extracted $(length(latencies)) events"

            for i in eachindex(latencies)
                latency = Int(round(latencies[i]))
                (latency < 1 || latency > n_timepoints) && continue

                # Extract trigger code from position field (use 0 for events without position)
                trigger_code = if positions !== nothing && i <= length(positions)
                    val = positions[i]
                    # Position is either a number or an empty array
                    (val isa Number) ? Int(val) : 0
                else
                    0
                end

                trigger_vec[latency] = trigger_code
                trigger_info_vec[latency] = string(types[i])
            end
        end
    end

    # Create DataFrame
    df = DataFrame(data_matrix, ch_names)
    insertcols!(df, 1, :time => times)
    insertcols!(df, 2, :triggers => trigger_vec)
    insertcols!(df, 3, :trigger_info => trigger_info_vec)

    # Parse layout
    layout = _parse_channel_locations(eeg["chanlocs"], ch_names)
    polar_to_cartesian_xy!(layout, preserve_radial_distance = preserve_radial_distance)

    # Create AnalysisInfo
    analysis_info = AnalysisInfo()

    # Create ContinuousData
    continuous_data = ContinuousData(filepath, df, layout, sample_rate, analysis_info)

    @info "Successfully converted to ContinuousData"
    return continuous_data
end


function _eeglab_to_erpdata(eeg::Dict, filepath::String, preserve_radial_distance::Bool)

    _, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    # Transpose to [timepoints × channels]
    data_matrix = data_array'

    # Create DataFrame
    df = DataFrame(data_matrix, ch_names)
    insertcols!(df, 1, :time => times)

    # Create AnalysisInfo (we don't really know much about the data)
    analysis_info = AnalysisInfo()

    # Get condition name (is this correct?)
    condition_name = get(eeg, "setname", basename(filepath))

    # Create ErpData
    # TODO: We should probably think about the condition number
    erp_data = ErpData(filepath, 1, condition_name, df, layout, sample_rate, analysis_info)

    @info "Successfully converted to ErpData"
    return erp_data
end


"""
    _read_eeglab_data(eeg_dict::Dict, filepath::String) → Array

Load EEG data from EEGLAB structure. Handles both embedded data arrays
and external .fdt (floating-point data) files.

# Arguments
- `eeg_dict::Dict`: EEGLAB data dictionary
- `filepath::String`: Path to the .set file (for resolving .fdt path)

# Returns
- Data array (channels × timepoints) or (channels × timepoints × trials)
"""
function _read_eeglab_data(eeg_dict::Dict, filepath::String)
    data_field = eeg_dict["data"]

    # Check if data is embedded (array) or external (filename string)
    if data_field isa AbstractArray
        # Data is embedded in .set file
        return data_field
    elseif data_field isa AbstractString && (endswith(lowercase(data_field), ".fdt") || endswith(lowercase(data_field), ".dat"))

        # Data is in external binary file (.fdt or .dat)
        fdt_filename = data_field

        # Construct full path to .fdt file (same directory as .set file)
        set_dir = dirname(filepath)
        fdt_path = joinpath(set_dir, fdt_filename)

        if !isfile(fdt_path)
            error("External data file not found: $fdt_path")
        end

        @info "Loading external data from: $fdt_filename"

        # Read binary .fdt file
        n_channels = Int(eeg_dict["nbchan"])
        n_timepoints = Int(eeg_dict["pnts"])
        n_trials = Int(eeg_dict["trials"])

        # Total expected elements
        n_elements = n_channels * n_timepoints * n_trials

        # Read binary file as Float32 array
        data_vec = open(fdt_path, "r") do io
            read!(io, Vector{Float32}(undef, n_elements))
        end

        # Reshape to [channels × timepoints × trials] or [channels × timepoints]
        if n_trials > 1
            return reshape(data_vec, n_channels, n_timepoints, n_trials)
        else
            return reshape(data_vec, n_channels, n_timepoints)
        end
    else
        error("Unexpected data field type: $(typeof(data_field))")
    end
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
the integer codes and original strings for storage in :triggers and :trigger_info columns.
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

Uses EEGLAB's native polar coordinates (theta, radius) which are what topoplot uses internally.
Falls back to 3D spherical or Cartesian conversion if polar coordinates are missing.
"""
function _parse_channel_locations(chanlocs, ch_names::Vector{Symbol})

    # EEGLAB set files have values for various coordinate systems
    # Let's try theta/radius first
    has_polar = haskey(chanlocs, "theta") && haskey(chanlocs, "radius")
    if has_polar && !isempty(chanlocs["theta"])

        theta = vec(chanlocs["theta"])
        radius = vec(chanlocs["radius"])

        # Check for valid polar coordinates
        valid_polar = .!isnan.(theta) .& .!isnan.(radius)

        if all(valid_polar) # convert to EegFun spherical coordinates
            azi = 90.0 .- theta
            inc = radius .* 180.0
            azi = mod.(azi, 360.0)
            @info "Loaded $(length(ch_names)) channels using native polar coordinates (theta, radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing)
        end

    end

    # Do we also need this if no theta/radius?
    has_spherical = haskey(chanlocs, "sph_theta") && haskey(chanlocs, "sph_phi") && haskey(chanlocs, "sph_radius")
    if has_spherical && !isempty(chanlocs["sph_theta"])
        sph_theta = vec(chanlocs["sph_theta"])
        sph_phi = vec(chanlocs["sph_phi"])
        sph_radius = vec(chanlocs["sph_radius"])

        # Check for valid spherical coordinates (all three components)
        valid_spherical = .!isnan.(sph_theta) .& .!isnan.(sph_phi) .& .!isnan.(sph_radius) .& (sph_radius .> 0)

        if all(valid_spherical) # Convert to EegFun spherical coordinates
            inc = 90.0 .- sph_phi
            azi = 90.0 .+ sph_theta
            azi = mod.(azi, 360.0)
            @info "Loaded $(length(ch_names)) channels using 3D spherical coordinates (sph_theta, sph_phi, sph_radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing)
        end

    end

    # TODO?: plotting eeglab x,y,z does not look correct for my test file. 
    # Actually, it also looks incorrect in eeglab itself.

    # TODO?: would a Coordinate conversion package e.g. CoordionateTransforations.jl be useful here?

end

"""
    _extract_ica_info(eeg::Dict, filepath::String, ch_names::Vector{Symbol}, layout::Layout) → InfoIca

Extract ICA information from EEGLAB structure.

Converts EEGLAB ICA matrices to EegFun InfoIca format.
"""
function _extract_ica_info(eeg::Dict, filepath::String, ch_names::Vector{Symbol}, layout::Layout)
    # Extract ICA matrices
    unmixing = Matrix{Float64}(eeg["icaweights"])  # Unmixing matrix
    sphere = Matrix{Float64}(eeg["icasphere"])     # Sphering matrix
    mixing = Matrix{Float64}(eeg["icawinv"])       # Mixing matrix (inverse)

    # Get number of components
    n_components = size(unmixing, 1)

    @info "Extracting $n_components ICA components"

    # Generate component labels
    ica_label = [Symbol("IC$i") for i = 1:n_components]

    # Create variance vector (placeholder - would need activations to compute properly)
    variance = ones(Float64, n_components)

    # TODO: can we actually get all the same info from *.set files here?

    # Create InfoIca object
    return InfoIca(
        filepath,                           # filename
        unmixing,                           # unmixing
        mixing,                             # mixing
        sphere,                             # sphere
        variance,                           # variance (placeholder)
        1.0,                                # scale
        zeros(Float64, n_components),       # mean
        ica_label,                          # ica_label
        OrderedDict{Int,Matrix{Float64}}(), # removed_activations (empty)
        layout,                             # layout
        fill(false, n_components),          # is_sub_gaussian (all super-Gaussian by default)
    )
end
