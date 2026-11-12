# EEGLAB .set file import functionality (experimental).

function _extract_channel_cols!(cols::Vector{AbstractVector}, offset::Int, mat::AbstractMatrix{T}, n_ch::Int, n_pts::Int) where {T}
    ch_vecs = Vector{Vector{T}}(undef, n_ch)
    @inbounds for c = 1:n_ch
        v = Vector{T}(undef, n_pts)
        ch_vecs[c] = v
        cols[offset+c] = v
    end
    tb = 64
    @inbounds for t0 = 1:tb:n_pts
        t1 = min(t0 + tb - 1, n_pts)
        for c = 1:n_ch
            dst = ch_vecs[c]
            @simd for t = t0:t1
                dst[t] = mat[c, t]
            end
        end
    end
    return cols
end

"""
    read_eeglab(filepath::String; preserve_radial_distance::Bool = true)

Load an EEGLAB `.set` file and return `ContinuousData`, `EpochData`, or `ErpData`
(or `(eeg_data, ica_data::InfoIca)` when ICA decomposition weights are present).

# Arguments
- `filepath::String`: Path to `.set` file
- `preserve_radial_distance::Bool`: Preserve radial distances for electrodes with `inc > 90°` (default: `true`)
"""
function read_eeglab(filepath::String; preserve_radial_distance::Bool = true)
    @info "Loading EEGLAB .set file: $filepath"
    mat_data = matread(filepath)
    eeg = haskey(mat_data, "EEG") ? mat_data["EEG"] : mat_data

    n_trials = Int(eeg["trials"])
    xmin = get(eeg, "xmin", 0.0)

    eeg_data = if n_trials == 1
        if xmin < 0
            @info "Detected averaged ERP data (xmin=$xmin < 0)"
            _eeglab_to_erpdata(eeg, filepath, preserve_radial_distance)
        else
            @info "Detected continuous data (xmin=$xmin >= 0)"
            _eeglab_to_continuousdata(eeg, filepath, preserve_radial_distance)
        end
    else
        _eeglab_to_epochdata(eeg, filepath, preserve_radial_distance)
    end

    if haskey(eeg, "icaweights") && !isempty(eeg["icaweights"])
        @info "ICA decomposition detected, loading ICA data"
        return eeg_data, _extract_ica_info(eeg, filepath, eeg_data.layout)
    end
    return eeg_data
end

function _load_common_eeglab_components(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    n_channels = Int(eeg["nbchan"])
    n_timepoints = Int(eeg["pnts"])
    sample_rate = round(Int, eeg["srate"])

    @info "Data: $n_channels channels × $n_timepoints timepoints" *
          (haskey(eeg, "trials") && eeg["trials"] > 1 ? " × $(Int(eeg["trials"])) trials" : "")
    @info "Sample rate: $sample_rate Hz"

    ch_names = _extract_channel_names(eeg["chanlocs"])
    data = _read_eeglab_data(eeg, filepath)

    xmin = Float64(get(eeg, "xmin", 0.0))
    times = if haskey(eeg, "times")
        [Float64(x) / 1000.0 for x in vec(eeg["times"])]
    else
        dt = 1.0 / sample_rate
        [xmin + (i - 1) * dt for i = 1:n_timepoints]
    end

    layout = _parse_channel_locations(eeg["chanlocs"], ch_names)
    if isnothing(layout)
        @minimal_warning "No complete coordinate system found. Layout will contain only channel names."
        layout =
            Layout(DataFrame(label = ch_names, inc = zeros(length(ch_names)), azi = zeros(length(ch_names))), nothing, nothing, nothing)
    else
        polar_to_cartesian_xy!(layout; preserve_radial_distance = preserve_radial_distance)
    end

    return n_timepoints, sample_rate, ch_names, data, times, layout
end

function _build_eeglab_epochs(
    data_array::AbstractArray{T,3},
    ch_names::Vector{Symbol},
    times::Vector{Float64},
    trigger_codes::Vector{Int},
    trigger_info::Vector{String},
    zero_idx::Int,
) where {T}
    n_ch, n_pts, n_trials = size(data_array)
    col_names = vcat([:time, :sample, :epoch, :trigger, :trigger_info], ch_names)
    sample_template = collect(1:n_pts)
    epoch_dfs = Vector{DataFrame}(undef, n_trials)

    @inbounds for trial_idx = 1:n_trials
        cols = Vector{AbstractVector}(undef, 5 + n_ch)
        cols[1] = copy(times)
        cols[2] = copy(sample_template)
        cols[3] = fill(trial_idx, n_pts)
        trig_vec = zeros(Int, n_pts)
        trig_vec[zero_idx] = trigger_codes[trial_idx]
        cols[4] = trig_vec
        info_vec = fill("", n_pts)
        info_vec[zero_idx] = trigger_info[trial_idx]
        cols[5] = info_vec
        _extract_channel_cols!(cols, 5, view(data_array, :, :, trial_idx), n_ch, n_pts)
        epoch_dfs[trial_idx] = DataFrame(cols, copy(col_names); copycols = false)
    end
    return epoch_dfs
end

function _eeglab_to_epochdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    _, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)
    n_trials = Int(eeg["trials"])

    trigger_codes, trigger_info =
        (haskey(eeg, "epoch") && !isempty(eeg["epoch"])) ? _extract_triggers(eeg["epoch"], n_trials) :
        (zeros(Int, n_trials), fill("", n_trials))

    zero_idx = find_closest_time_index(times, 0.0)
    epoch_dataframes = _build_eeglab_epochs(data_array, ch_names, times, trigger_codes, trigger_info, zero_idx)

    condition_name = get(eeg, "setname", basename(filepath))
    analysis_info = AnalysisInfo(sample_rate = Float64(sample_rate))
    epoch_data = EpochData(filepath, 1, condition_name, epoch_dataframes, layout, sample_rate, analysis_info)
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

@inline _unwrap_mat_scalar(v, default) = v isa AbstractArray ? (isempty(v) ? default : first(v)) : v

function _record_eeglab_event!(trigger_vec::Vector{Int}, trigger_info_vec::Vector{String}, lat_val, typ_val, pos_val, n_timepoints::Int)
    (lat_val isa Number && !isnan(lat_val)) || return
    latency = round(Int, lat_val)
    (1 <= latency <= n_timepoints) || return

    code = _eeglab_event_code(typ_val, pos_val)
    if code != 0 || trigger_vec[latency] == 0
        trigger_vec[latency] = code
    end
    info_str = (typ_val isa AbstractArray && isempty(typ_val)) ? "" : string(typ_val)
    if code != 0 || isempty(trigger_info_vec[latency])
        trigger_info_vec[latency] = info_str
    end
end

function _record_eeglab_struct_events!(
    trigger_vec::Vector{Int},
    trigger_info_vec::Vector{String},
    latencies::AbstractVector,
    types,
    positions,
    n_timepoints::Int,
)
    for i in eachindex(latencies)
        typ_val = (!isnothing(types) && i <= length(types)) ? types[i] : ""
        pos_val = (!isnothing(positions) && i <= length(positions)) ? positions[i] : nothing
        _record_eeglab_event!(trigger_vec, trigger_info_vec, latencies[i], typ_val, pos_val, n_timepoints)
    end
end

function _eeglab_extract_continuous_events!(trigger_vec::Vector{Int}, trigger_info_vec::Vector{String}, events, n_timepoints::Int)
    if events isa AbstractArray
        for ev in events
            (ev isa AbstractDict && haskey(ev, "latency")) || continue
            lat_val = _unwrap_mat_scalar(ev["latency"], nothing)
            typ_val = _unwrap_mat_scalar(get(ev, "type", ""), "")
            pos_val = _unwrap_mat_scalar(get(ev, "position", nothing), nothing)
            _record_eeglab_event!(trigger_vec, trigger_info_vec, lat_val, typ_val, pos_val, n_timepoints)
        end
    elseif haskey(events, "latency") && (haskey(events, "type") || haskey(events, "position"))
        latencies = _as_eeglab_vec(events["latency"])
        types = haskey(events, "type") ? _as_eeglab_vec(events["type"]) : nothing
        positions = haskey(events, "position") ? _as_eeglab_vec(events["position"]) : nothing
        @info "Extracted $(length(latencies)) events"
        _record_eeglab_struct_events!(trigger_vec, trigger_info_vec, latencies, types, positions, n_timepoints)
    end
end

_as_eeglab_vec(x) = x isa AbstractArray ? vec(x) : [x]

function _eeglab_to_continuousdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    trigger_vec = zeros(Int, n_timepoints)
    trigger_info_vec = fill("", n_timepoints)
    if haskey(eeg, "event") && !isempty(eeg["event"])
        _eeglab_extract_continuous_events!(trigger_vec, trigger_info_vec, eeg["event"], n_timepoints)
    end

    n_ch = length(ch_names)
    cols = Vector{AbstractVector}(undef, 4 + n_ch)
    cols[1] = times
    cols[2] = collect(1:n_timepoints)
    cols[3] = trigger_vec
    cols[4] = trigger_info_vec
    _extract_channel_cols!(cols, 4, data_array, n_ch, n_timepoints)
    df = DataFrame(cols, vcat([:time, :sample, :trigger, :trigger_info], ch_names); copycols = false)

    continuous_data = ContinuousData(filepath, df, layout, sample_rate, AnalysisInfo(sample_rate = Float64(sample_rate)))
    @info "Successfully converted to ContinuousData"
    return continuous_data
end

function _eeglab_to_erpdata(eeg::AbstractDict, filepath::String, preserve_radial_distance::Bool)
    n_timepoints, sample_rate, ch_names, data_array, times, layout = _load_common_eeglab_components(eeg, filepath, preserve_radial_distance)

    n_ch = length(ch_names)
    cols = Vector{AbstractVector}(undef, 2 + n_ch)
    cols[1] = times
    cols[2] = collect(1:n_timepoints)
    _extract_channel_cols!(cols, 2, data_array, n_ch, n_timepoints)
    df = DataFrame(cols, vcat([:time, :sample], ch_names); copycols = false)

    condition_name = get(eeg, "setname", basename(filepath))
    erp_data = ErpData(filepath, 1, condition_name, df, layout, sample_rate, AnalysisInfo(sample_rate = Float64(sample_rate)), 1)
    @info "Successfully converted to ErpData"
    return erp_data
end

function _read_eeglab_data(eeg_dict::AbstractDict, filepath::String)
    data_field = eeg_dict["data"]
    data_field isa AbstractArray && return data_field

    if data_field isa AbstractString && (endswith(lowercase(data_field), ".fdt") || endswith(lowercase(data_field), ".dat"))
        fdt_path = joinpath(dirname(filepath), data_field)
        !isfile(fdt_path) && @minimal_error("External data file not found: $fdt_path")
        @info "Loading external data from: $data_field"

        n_channels = Int(eeg_dict["nbchan"])
        n_timepoints = Int(eeg_dict["pnts"])
        n_trials = Int(eeg_dict["trials"])
        data_vec = open(fdt_path, "r") do io
            read!(io, Vector{Float32}(undef, n_channels * n_timepoints * n_trials))
        end
        return n_trials > 1 ? reshape(data_vec, n_channels, n_timepoints, n_trials) : reshape(data_vec, n_channels, n_timepoints)
    end
    @minimal_error("Unexpected data field type: $(typeof(data_field))")
end

function _extract_channel_names(chanlocs)
    isempty(chanlocs) && @minimal_error("No channel location information found in .set file")
    labels = chanlocs["labels"]
    labels isa Union{AbstractString,Symbol} && return [Symbol(labels)]
    labels isa AbstractVector && return Symbol.(labels)
    labels isa AbstractMatrix && return Symbol.(vec(labels))
    @minimal_error("Unexpected chanlocs.labels format")
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
    end
    return string(event_types)
end

function _extract_triggers(epochs, n_trials::Int)
    trigger_strings = Vector{String}(undef, n_trials)

    if epochs isa AbstractDict
        ets = haskey(epochs, "eventtype") ? vec(epochs["eventtype"]) : nothing
        lats = haskey(epochs, "eventlatency") ? vec(epochs["eventlatency"]) : nothing
        for i = 1:n_trials
            et_i = (!isnothing(ets) && i <= length(ets)) ? ets[i] : nothing
            lat_i = (!isnothing(lats) && i <= length(lats)) ? lats[i] : nothing
            trigger_strings[i] = _select_epoch_event_type(et_i, lat_i)
        end
    else
        for i = 1:n_trials
            epoch = epochs[i]
            et_i = haskey(epoch, "eventtype") ? epoch["eventtype"] : nothing
            lat_i = haskey(epoch, "eventlatency") ? epoch["eventlatency"] : nothing
            trigger_strings[i] = _select_epoch_event_type(et_i, lat_i)
        end
    end

    parsed_codes = [_parse_marker_code(s) for s in trigger_strings]
    max_num = isempty(parsed_codes) ? 0 : max(0, maximum(parsed_codes))
    non_numeric = sort!(unique!([s for (s, c) in zip(trigger_strings, parsed_codes) if !isempty(strip(s)) && c == 0]))
    label_to_int = Dict{String,Int}(lbl => (max_num + idx) for (idx, lbl) in enumerate(non_numeric))

    trigger_codes = Vector{Int}(undef, n_trials)
    for i = 1:n_trials
        s = trigger_strings[i]
        trigger_codes[i] = isempty(strip(s)) ? 0 : (parsed_codes[i] != 0 ? parsed_codes[i] : label_to_int[s])
    end

    unique_vals = sort!(unique(trigger_strings))
    @info "Extracted $(length(unique_vals)) unique trigger types: $(unique_vals)"
    return trigger_codes, trigger_strings
end

function _parse_channel_locations(chanlocs, ch_names::AbstractVector{Symbol})
    if haskey(chanlocs, "theta") && haskey(chanlocs, "radius") && !isempty(chanlocs["theta"])
        theta = [x isa Number ? Float64(x) : NaN for x in vec(chanlocs["theta"])]
        radius = [x isa Number ? Float64(x) : NaN for x in vec(chanlocs["radius"])]
        if any(i -> !isnan(theta[i]) && !isnan(radius[i]), eachindex(theta, radius))
            inc = radius .* 180.0
            azi = mod.(90.0 .- theta, 360.0)
            @info "Loaded $(length(ch_names)) channels using native polar coordinates (theta, radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing, nothing)
        end
    end

    if haskey(chanlocs, "sph_theta") && haskey(chanlocs, "sph_phi") && haskey(chanlocs, "sph_radius") && !isempty(chanlocs["sph_theta"])
        sph_theta = [x isa Number ? Float64(x) : NaN for x in vec(chanlocs["sph_theta"])]
        sph_phi = [x isa Number ? Float64(x) : NaN for x in vec(chanlocs["sph_phi"])]
        sph_radius = [x isa Number ? Float64(x) : NaN for x in vec(chanlocs["sph_radius"])]
        if any(i -> !isnan(sph_theta[i]) && !isnan(sph_phi[i]) && !isnan(sph_radius[i]) && sph_radius[i] > 0, eachindex(sph_theta))
            inc = 90.0 .- sph_phi
            azi = mod.(90.0 .+ sph_theta, 360.0)
            @info "Loaded $(length(ch_names)) channels using 3D spherical coordinates (sph_theta, sph_phi, sph_radius)"
            return Layout(DataFrame(label = ch_names, inc = inc, azi = azi), nothing, nothing, nothing)
        end
    end
    return nothing
end

function _extract_ica_info(eeg::AbstractDict, filepath::String, layout::Layout)
    W = Matrix{Float64}(eeg["icaweights"])
    has_sphere = haskey(eeg, "icasphere") && !isempty(eeg["icasphere"])
    sphere = has_sphere ? Matrix{Float64}(eeg["icasphere"]) : Matrix{Float64}(LinearAlgebra.I, size(W, 2), size(W, 2))
    unmixing = (has_sphere && size(W, 2) == size(sphere, 1)) ? (W * sphere) : W
    mixing =
        (haskey(eeg, "icawinv") && !isempty(eeg["icawinv"])) ? Matrix{Float64}(eeg["icawinv"]) :
        Matrix{Float64}(LinearAlgebra.pinv(unmixing))

    n_components, n_ica_channels = size(unmixing)
    @info "Extracting $n_components ICA components across $n_ica_channels channels"

    ica_layout = copy(layout)
    if haskey(eeg, "icachansind") && !isempty(eeg["icachansind"])
        raw_inds = vec(eeg["icachansind"])
        if all(x -> x isa Number && !isnan(x), raw_inds)
            ch_inds = [round(Int, x) for x in raw_inds]
            all_lbls = ica_layout.data.label
            if all(idx -> 1 <= idx <= length(all_lbls), ch_inds) && length(ch_inds) <= length(all_lbls)
                sel_labels = all_lbls[ch_inds]
                if sel_labels != all_lbls
                    ica_layout = subset_layout(ica_layout; channel_selection = channels(sel_labels))
                end
            end
        end
    end

    return InfoIca(
        filepath,
        unmixing,
        mixing,
        sphere,
        ones(Float64, n_components),
        1.0,
        zeros(Float64, n_ica_channels),
        [Symbol("IC$i") for i = 1:n_components],
        Dict{Int,Matrix{Float64}}(),
        ica_layout,
        fill(false, n_components),
    )
end
