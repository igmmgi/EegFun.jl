"""
EEGLAB .set file import functionality.

Provides functions to read EEGLAB .set files (MATLAB format) and convert them
to EegFun datatypes.

Note: This EEGLAB import functionality is considered experimental. It has been
validated against standard EEGLAB sample datasets, but some edge-case structures
may not be fully supported. Please report any issues or missing features.
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
        if isnothing(layout)
            @minimal_warning "No valid coordinate system found for ICA. Layout will contain only channel names."
            layout =
                Layout(DataFrame(label = ch_names, inc = zeros(length(ch_names)), azi = zeros(length(ch_names))), nothing, nothing, nothing)
        else
            polar_to_cartesian_xy!(layout, preserve_radial_distance = preserve_radial_distance)
        end
        ica_data = _extract_ica_info(eeg, filepath, layout)
        return eeg_data, ica_data
    else
        return eeg_data
    end
end


"""
Helper to extract common components used by all EEGLAB converters.
Returns: (n_channels, n_timepoints, sample_rate, ch_names, data, times, layout)
"""
function _load_common_eeglab_components(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    n_channels = Int(eeg["nbchan"])
    n_timepoints = Int(eeg["pnts"])
    sample_rate = Int(round(eeg["srate"]))

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
    if isnothing(layout)
        @minimal_warning "No complete coordinate system found. Layout will contain only channel names."
        layout =
            Layout(DataFrame(label = ch_names, inc = zeros(length(ch_names)), azi = zeros(length(ch_names))), nothing, nothing, nothing)
    else
        polar_to_cartesian_xy!(layout, preserve_radial_distance = preserve_radial_distance)
    end

    return n_timepoints, sample_rate, ch_names, data, times, layout
end


"""
    _eeglab_to_epochdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool) → EpochData

Convert EEGLAB structure to EegFun EpochData.

Internal function called by `load_eeglab`.
"""
function _eeglab_to_epochdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    # Load common components
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    n_trials = Int(eeg["trials"])

    # Extract triggers from epochs (returns both integer codes and string info)
    trigger_codes, trigger_info =
        haskey(eeg, "epoch") && !isempty(eeg["epoch"]) ? _extract_triggers(eeg["epoch"], n_trials) :
        (zeros(Int, n_trials), fill("", n_trials))

    # Convert 3D array to Vector{DataFrame}
    epoch_dataframes = Vector{DataFrame}(undef, n_trials)
    zero_idx = find_closest_time_index(times, 0.0)
    sample_col = collect(1:n_timepoints)

    for trial_idx = 1:n_trials
        # Extract trial data [channels × timepoints]
        trial_data = data_array[:, :, trial_idx]'  # Transpose to [timepoints × channels]

        # Create DataFrame with channel data
        df = DataFrame(trial_data, ch_names)

        # Create trigger vectors (0 everywhere except at t=0)
        trigger_vec = zeros(Int, n_timepoints)
        trigger_vec[zero_idx] = trigger_codes[trial_idx]

        trigger_info_vec = fill("", n_timepoints)
        trigger_info_vec[zero_idx] = trigger_info[trial_idx]

        insertcols!(
            df,
            1,
            :time => times,
            :sample => copy(sample_col),
            :epoch => fill(trial_idx, n_timepoints),
            :trigger => trigger_vec,
            :trigger_info => trigger_info_vec,
        )

        epoch_dataframes[trial_idx] = df
    end

    # Extract condition name if available
    condition_name = get(eeg, "setname", basename(filepath))

    # Create AnalysisInfo (EEGLAB doesn't have reliable filter/reference metadata)
    analysis_info = AnalysisInfo(sample_rate = Float64(sample_rate))

    # Create EpochData
    epoch_data = EpochData(
        filepath,          # file
        1,                 # condition
        condition_name,    # condition_name
        epoch_dataframes,  # data
        copy(layout),      # layout
        sample_rate,       # sample_rate
        analysis_info,     # analysis_info
    )

    @info "Successfully converted to EpochData with $n_trials epochs"

    return epoch_data
end


function _eeglab_event_code(typ, pos)::Int
    if typ isa Integer && !(typ isa Bool)
        return Int(typ)
    elseif typ isa Real && !(typ isa Bool) && !isnan(typ)
        return round(Int, typ)
    elseif !(typ isa AbstractArray && isempty(typ))
        parsed = _parse_marker_code(string(typ))
        parsed != 0 && return parsed
    end
    if pos isa Integer && !(pos isa Bool)
        return Int(pos)
    elseif pos isa Real && !(pos isa Bool) && !isnan(pos)
        return round(Int, pos)
    end
    return 0
end

function _eeglab_extract_continuous_events!(trigger_vec::Vector{Int}, trigger_info_vec::Vector{String}, events, n_timepoints::Int)
    if events isa AbstractArray
        # Array-of-dicts representation
        for ev in events
            ev isa AbstractDict || continue
            haskey(ev, "latency") || continue
            lat_raw = ev["latency"]
            lat_val = lat_raw isa AbstractArray ? (isempty(lat_raw) ? nothing : first(lat_raw)) : lat_raw
            (lat_val isa Number && !isnan(lat_val)) || continue
            latency = Int(round(lat_val))
            (1 <= latency <= n_timepoints) || continue

            typ = get(ev, "type", "")
            typ_val = typ isa AbstractArray ? (isempty(typ) ? "" : first(typ)) : typ
            pos = get(ev, "position", nothing)
            pos_val = pos isa AbstractArray ? (isempty(pos) ? nothing : first(pos)) : pos

            code = _eeglab_event_code(typ_val, pos_val)
            if code != 0 || trigger_vec[latency] == 0
                trigger_vec[latency] = code
            end
            info_str = string(typ_val)
            if code != 0 || isempty(trigger_info_vec[latency])
                trigger_info_vec[latency] = info_str
            end
        end
    else
        # Struct-of-arrays representation (MAT.MatlabStructArray or Dict of arrays)
        has_lat = haskey(events, "latency")
        has_type = haskey(events, "type")
        if has_lat && (has_type || haskey(events, "position"))
            latencies = vec(events["latency"])
            types = has_type ? vec(events["type"]) : nothing
            positions = haskey(events, "position") ? vec(events["position"]) : nothing

            @info "Extracted $(length(latencies)) events"

            for i in eachindex(latencies)
                lat_val = latencies[i]
                (lat_val isa Number && !isnan(lat_val)) || continue
                latency = Int(round(lat_val))
                (1 <= latency <= n_timepoints) || continue

                typ_val = (!isnothing(types) && i <= length(types)) ? types[i] : ""
                pos_val = (!isnothing(positions) && i <= length(positions)) ? positions[i] : nothing

                code = _eeglab_event_code(typ_val, pos_val)
                if code != 0 || trigger_vec[latency] == 0
                    trigger_vec[latency] = code
                end
                info_str = typ_val isa AbstractArray && isempty(typ_val) ? "" : string(typ_val)
                if code != 0 || isempty(trigger_info_vec[latency])
                    trigger_info_vec[latency] = info_str
                end
            end
        end
    end
end

function _eeglab_to_continuousdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    # Load common components
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    # Transpose to [timepoints × channels]
    data_matrix = data_array'

    # Extract triggers/events if present
    trigger_vec = zeros(Int, n_timepoints)
    trigger_info_vec = fill("", n_timepoints)

    if haskey(eeg, "event") && !isempty(eeg["event"])
        _eeglab_extract_continuous_events!(trigger_vec, trigger_info_vec, eeg["event"], n_timepoints)
    end

    # Create DataFrame
    df = DataFrame(data_matrix, ch_names)
    insertcols!(df, 1, :time => times, :sample => collect(1:n_timepoints), :trigger => trigger_vec, :trigger_info => trigger_info_vec)

    # Create AnalysisInfo
    analysis_info = AnalysisInfo(sample_rate = Float64(sample_rate))

    # Create ContinuousData
    continuous_data = ContinuousData(filepath, df, copy(layout), sample_rate, analysis_info)

    @info "Successfully converted to ContinuousData"
    return continuous_data
end


function _eeglab_to_erpdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)

    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    # Transpose to [timepoints × channels]
    data_matrix = data_array'

    # Create DataFrame
    df = DataFrame(data_matrix, ch_names)
    insertcols!(df, 1, :time => times, :sample => collect(1:n_timepoints))

    # Create AnalysisInfo (we don't really know much about the data)
    analysis_info = AnalysisInfo(sample_rate = Float64(sample_rate))

    # Get condition name
    condition_name = get(eeg, "setname", basename(filepath))

    # Create ErpData
    # Note: EEGLAB averaged ERPs don't natively store epoch counts in a standard way,
    # so we default to n_epochs=1 and condition=1.
    erp_data = ErpData(filepath, 1, condition_name, df, copy(layout), sample_rate, analysis_info, 1)

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
            @minimal_error("External data file not found: $fdt_path")
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
        @minimal_error("Unexpected data field type: $(typeof(data_field))")
    end
end


"""
    _extract_channel_names(chanlocs) → Vector{Symbol}


Extract channel names from EEGLAB chanlocs structure.
"""
function _extract_channel_names(chanlocs)
    if isempty(chanlocs)
        @minimal_error("No channel location information found in .set file")
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
        @minimal_error("Unexpected chanlocs.labels format")
    end
end


function _select_epoch_event_type(event_types, event_latencies)::String
    (event_types === nothing || (event_types isa AbstractArray && isempty(event_types))) && return ""
    if event_types isa AbstractArray
        types_vec = vec(event_types)
        if length(types_vec) > 1 && event_latencies isa AbstractArray && length(event_latencies) == length(types_vec)
            lats_vec = vec(event_latencies)
            best_idx = 1
            best_abs = Inf
            for (k, lat) in enumerate(lats_vec)
                lat_num = lat isa Number ? Float64(lat) : tryparse(Float64, string(lat))
                if !isnothing(lat_num) && !isnan(lat_num) && abs(lat_num) < best_abs
                    best_abs = abs(lat_num)
                    best_idx = k
                end
            end
            return string(types_vec[best_idx])
        end
        return string(first(types_vec))
    else
        return string(event_types)
    end
end

"""
    _extract_triggers(epochs, n_trials::Int) → Tuple{Vector{Int}, Vector{String}}

Extract trigger codes from EEGLAB epoch structure, selecting the event closest to `eventlatency == 0`
when an epoch contains multiple events. Handles both struct-of-arrays and array-of-dicts representations.
"""
function _extract_triggers(epochs, n_trials::Int)
    trigger_strings = Vector{String}(undef, n_trials)

    if epochs isa AbstractDict
        # Dict struct-of-arrays representation
        has_et = haskey(epochs, "eventtype")
        has_lat = haskey(epochs, "eventlatency")
        ets = has_et ? vec(epochs["eventtype"]) : nothing
        lats = has_lat ? vec(epochs["eventlatency"]) : nothing
        for i = 1:n_trials
            et_i = (!isnothing(ets) && i <= length(ets)) ? ets[i] : nothing
            lat_i = (!isnothing(lats) && i <= length(lats)) ? lats[i] : nothing
            trigger_strings[i] = _select_epoch_event_type(et_i, lat_i)
        end
    else
        # MAT.MatlabStructArray or Vector{Dict} representation
        for i = 1:n_trials
            epoch = epochs[i]
            et_i = haskey(epoch, "eventtype") ? epoch["eventtype"] : nothing
            lat_i = haskey(epoch, "eventlatency") ? epoch["eventlatency"] : nothing
            trigger_strings[i] = _select_epoch_event_type(et_i, lat_i)
        end
    end

    # Parse numeric codes directly when present; assign deterministic sequential codes to purely non-numeric labels
    parsed_codes = [_parse_marker_code(s) for s in trigger_strings]
    max_num = isempty(parsed_codes) ? 0 : max(0, maximum(parsed_codes))
    non_numeric = sort(unique([s for (s, c) in zip(trigger_strings, parsed_codes) if !isempty(strip(s)) && c == 0]))
    label_to_int = Dict{String,Int}(lbl => (max_num + idx) for (idx, lbl) in enumerate(non_numeric))

    trigger_codes = Vector{Int}(undef, n_trials)
    for i = 1:n_trials
        s = trigger_strings[i]
        if isempty(strip(s))
            trigger_codes[i] = 0
        elseif parsed_codes[i] != 0
            trigger_codes[i] = parsed_codes[i]
        else
            trigger_codes[i] = label_to_int[s]
        end
    end

    unique_vals = sort(unique(trigger_strings))
    @info "Extracted $(length(unique_vals)) unique trigger types: $(unique_vals)"

    return trigger_codes, trigger_strings
end


"""
    _parse_channel_locations(chanlocs, ch_names::Vector{Symbol}) → Layout

Parse EEGLAB channel locations into EegFun Layout.

Uses EEGLAB's native polar coordinates (theta, radius) which are what topoplot uses internally.
Falls back to 3D spherical or Cartesian conversion if polar coordinates are missing.
"""
function _parse_channel_locations(chanlocs, ch_names::AbstractVector{Symbol})

    # EEGLAB set files have values for various coordinate systems
    # Let's try theta/radius first
    has_polar = haskey(chanlocs, "theta") && haskey(chanlocs, "radius")
    if has_polar && !isempty(chanlocs["theta"])

        # Vector extraction mapping empty matrices to NaN
        theta = map(x -> x isa Number ? Float64(x) : NaN, vec(chanlocs["theta"]))
        radius = map(x -> x isa Number ? Float64(x) : NaN, vec(chanlocs["radius"]))

        # Check for valid polar coordinates
        valid_polar = .!isnan.(theta) .& .!isnan.(radius)

        if any(valid_polar) # convert to EegFun spherical coordinates
            azi = 90.0 .- theta
            inc = radius .* 180.0
            azi = mod.(azi, 360.0)
            @info "Loaded $(length(ch_names)) channels using native polar coordinates (theta, radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing, nothing)
        end

    end

    # Do we also need this if no theta/radius?
    has_spherical = haskey(chanlocs, "sph_theta") && haskey(chanlocs, "sph_phi") && haskey(chanlocs, "sph_radius")
    if has_spherical && !isempty(chanlocs["sph_theta"])
        # Vector extraction mapping empty matrices to NaN
        sph_theta = map(x -> x isa Number ? Float64(x) : NaN, vec(chanlocs["sph_theta"]))
        sph_phi = map(x -> x isa Number ? Float64(x) : NaN, vec(chanlocs["sph_phi"]))
        sph_radius = map(x -> x isa Number ? Float64(x) : NaN, vec(chanlocs["sph_radius"]))

        # Check for valid spherical coordinates (all three components)
        valid_spherical = .!isnan.(sph_theta) .& .!isnan.(sph_phi) .& .!isnan.(sph_radius) .& (sph_radius .> 0)

        if any(valid_spherical) # Convert to EegFun spherical coordinates
            inc = 90.0 .- sph_phi
            azi = 90.0 .+ sph_theta
            azi = mod.(azi, 360.0)
            @info "Loaded $(length(ch_names)) channels using 3D spherical coordinates (sph_theta, sph_phi, sph_radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing, nothing)
        end

    end

    # Note: Plotting EEGLAB Cartesian x,y,z coordinates directly often appears distorted
    # due to different coordinate system conventions between EEGLAB and EegFun.
    # We rely on polar/spherical coordinates above which are more robustly standardized.

end

"""
    _extract_ica_info(eeg::Dict, filepath::String, layout::Layout) → InfoIca

Extract ICA information from EEGLAB structure.

Converts EEGLAB ICA matrices to EegFun InfoIca format.
"""
function _extract_ica_info(eeg::AbstractDict, filepath::String, layout::Layout)
    # Extract ICA matrices: in EEGLAB, full unmixing = icaweights * icasphere
    W = Matrix{Float64}(eeg["icaweights"])
    has_sphere = haskey(eeg, "icasphere") && !isempty(eeg["icasphere"])
    sphere = has_sphere ? Matrix{Float64}(eeg["icasphere"]) : [i == j ? 1.0 : 0.0 for i = 1:size(W, 2), j = 1:size(W, 2)]
    unmixing = (has_sphere && size(W, 2) == size(sphere, 1)) ? (W * sphere) : W
    mixing =
        (haskey(eeg, "icawinv") && !isempty(eeg["icawinv"])) ? Matrix{Float64}(eeg["icawinv"]) :
        Matrix{Float64}(LinearAlgebra.pinv(unmixing))

    # Get number of components and ICA channels
    n_components = size(unmixing, 1)
    n_ica_channels = size(unmixing, 2)

    @info "Extracting $n_components ICA components across $n_ica_channels channels"

    # Subset copied layout to icachansind if present
    ica_layout = copy(layout)
    if haskey(eeg, "icachansind") && !isempty(eeg["icachansind"])
        raw_inds = vec(eeg["icachansind"])
        if all(x -> x isa Number && !isnan(x), raw_inds)
            ch_inds = Int.(round.(raw_inds))
            all_lbls = ica_layout.data.label
            if all(idx -> 1 <= idx <= length(all_lbls), ch_inds) && length(ch_inds) <= length(all_lbls)
                sel_labels = all_lbls[ch_inds]
                if sel_labels != all_lbls
                    ica_layout = subset_layout(ica_layout; channel_selection = channels(sel_labels))
                end
            end
        end
    end

    # Generate component labels
    ica_label = [Symbol("IC$i") for i = 1:n_components]

    # Create variance vector (placeholder - would need activations to compute properly)
    variance = ones(Float64, n_components)

    # Create InfoIca object
    return InfoIca(
        filepath,                           # filename
        unmixing,                           # unmixing (W * sphere)
        mixing,                             # mixing
        sphere,                             # sphere
        variance,                           # variance (placeholder)
        1.0,                                # scale
        zeros(Float64, n_ica_channels),     # mean (length = size(unmixing, 2))
        ica_label,                          # ica_label
        Dict{Int,Matrix{Float64}}(),        # removed_activations (empty)
        ica_layout,                         # layout (copied & subsetted to icachansind)
        fill(false, n_components),          # is_sub_gaussian (all super-Gaussian by default)
    )
end
