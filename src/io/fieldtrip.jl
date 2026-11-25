# FieldTrip `.mat` data structure import functionality.

######################################
########### User Functions ###########
######################################

"""
    read_fieldtrip(filepath::AbstractString, layout::Layout)

Load a FieldTrip `.mat` file and convert it to `ContinuousData`, `EpochData`, or `ErpData`.

# Arguments
- `filepath::AbstractString`: Path to FieldTrip `.mat` file
- `layout::Layout`: Channel layout
"""
function read_fieldtrip(filepath::AbstractString, layout::Layout)
    filepath_str = String(filepath)
    @info "Loading FieldTrip .mat file: $filepath_str"
    mat_data = matread(filepath_str)

    ft = nothing
    is_erp = false
    for (key, value) in mat_data
        value isa AbstractDict || continue
        if haskey(value, "avg") && haskey(value, "time")
            ft = value
            is_erp = true
            @info "Found FieldTrip ERP data in variable: $key"
            break
        elseif haskey(value, "trial") && haskey(value, "fsample")
            ft = value
            @info "Found FieldTrip data in variable: $key"
            break
        end
    end

    isnothing(ft) &&
        @minimal_error("Could not find FieldTrip data structure in file. Expected a dict with 'trial'+'fsample' or 'avg'+'time' fields.",)

    if is_erp
        time_vec = vec(ft["time"])
        length(time_vec) > 1 ||
            @minimal_error("Invalid ERP data: time vector has only $(length(time_vec)) point(s). Cannot infer sample rate.")
        sample_rate = round(Int, 1.0 / (time_vec[2] - time_vec[1]))
        @info "Detected averaged ERP data, inferred sample rate: $sample_rate Hz"
        return _fieldtrip_to_erpdata(ft, filepath_str, layout, sample_rate)
    end

    n_trials = length(ft["trial"])
    sample_rate = round(Int, ft["fsample"])
    @info "Found $n_trials trial(s), sample rate: $sample_rate Hz"
    if n_trials == 1
        @info "Detected continuous data"
        return _fieldtrip_to_continuousdata(ft, filepath_str, layout, sample_rate)
    end
    @info "Detected epoched data with $n_trials trials"
    return _fieldtrip_to_epochdata(ft, filepath_str, layout, sample_rate)
end

######################################
######### Internal Functions #########
######################################

function _extract_ft_channel_names(ch_labels::AbstractVector)::Vector{Symbol}
    if !isempty(ch_labels) && lowercase(string(ch_labels[end])) == "status"
        @info "Detected BDF Status channel, removing it"
        return [Symbol(ch_labels[i]) for i = 1:(length(ch_labels)-1)]
    end
    return [Symbol(c) for c in ch_labels]
end

function _extract_time_vector(time_data)::Vector{Float64}
    raw = (time_data isa AbstractMatrix && eltype(time_data) <: Number && size(time_data, 1) == 1) ? vec(time_data) : vec(time_data[1])
    return raw isa Vector{Float64} ? raw : [Float64(x) for x in raw]
end

_as_ft_vec(x) = x isa AbstractArray ? vec(x) : Any[x]

function _set_ft_event!(trigger_vec::Vector{Int}, trigger_info_vec::Vector{String}, n_timepoints::Int, raw_sample, raw_val, raw_typ)
    s_raw = raw_sample isa AbstractArray ? (isempty(raw_sample) ? nothing : first(raw_sample)) : raw_sample
    (s_raw isa Number && isfinite(s_raw)) || return nothing
    s = round(Int, s_raw)
    (1 <= s <= n_timepoints) || return nothing

    if !isnothing(raw_val)
        v = raw_val isa AbstractArray ? (isempty(raw_val) ? nothing : first(raw_val)) : raw_val
        if v isa Number && isfinite(v)
            trigger_vec[s] = round(Int, v)
        elseif v isa AbstractString
            code = _parse_marker_code(v)
            code != 0 && (trigger_vec[s] = code)
        end
    end

    if !isnothing(raw_typ)
        trigger_info_vec[s] = raw_typ isa AbstractArray ? (isempty(raw_typ) ? "" : string(first(raw_typ))) : string(raw_typ)
    end
    return nothing
end

function _extract_ft_events!(trigger_vec::Vector{Int}, trigger_info_vec::Vector{String}, events, n_timepoints::Int)
    if events isa AbstractDict
        haskey(events, "sample") || return nothing
        samples = _as_ft_vec(events["sample"])
        @info "Extracting $(length(samples)) events from cfg.event"
        vals = haskey(events, "value") ? _as_ft_vec(events["value"]) : nothing
        types = haskey(events, "type") ? _as_ft_vec(events["type"]) : nothing
        for i in eachindex(samples)
            val = (!isnothing(vals) && i <= length(vals)) ? vals[i] : nothing
            typ = (!isnothing(types) && i <= length(types)) ? types[i] : nothing
            _set_ft_event!(trigger_vec, trigger_info_vec, n_timepoints, samples[i], val, typ)
        end
    elseif events isa AbstractArray
        @info "Extracting $(length(events)) events from cfg.event"
        for ev in events
            ev isa AbstractDict || continue
            _set_ft_event!(
                trigger_vec,
                trigger_info_vec,
                n_timepoints,
                get(ev, "sample", nothing),
                get(ev, "value", nothing),
                get(ev, "type", nothing),
            )
        end
    elseif hasproperty(events, :names) && hasproperty(events, :values)
        sample_idx = findfirst(==("sample"), events.names)
        isnothing(sample_idx) && return nothing
        samples = _as_ft_vec(events.values[sample_idx])
        value_idx = findfirst(==("value"), events.names)
        type_idx = findfirst(==("type"), events.names)
        vals = !isnothing(value_idx) ? _as_ft_vec(events.values[value_idx]) : nothing
        types = !isnothing(type_idx) ? _as_ft_vec(events.values[type_idx]) : nothing
        for i in eachindex(samples)
            val = (!isnothing(vals) && i <= length(vals)) ? vals[i] : nothing
            typ = (!isnothing(types) && i <= length(types)) ? types[i] : nothing
            _set_ft_event!(trigger_vec, trigger_info_vec, n_timepoints, samples[i], val, typ)
        end
    end
    return nothing
end

function _fieldtrip_to_continuousdata(ft::AbstractDict, filepath::String, layout::Layout, sample_rate::Int)
    ch_names = _extract_ft_channel_names(vec(ft["label"]))
    data_matrix = ft["trial"][1]
    time_raw = vec(ft["time"][1])
    time_vec = time_raw isa Vector{Float64} ? copy(time_raw) : [Float64(x) for x in time_raw]

    n_channels = length(ch_names)
    n_timepoints = size(data_matrix, 2)
    @info "Data: $n_channels channels × $n_timepoints timepoints"

    trigger_vec = zeros(Int, n_timepoints)
    trigger_info_vec = fill("", n_timepoints)

    if haskey(ft, "cfg") && ft["cfg"] isa AbstractDict && haskey(ft["cfg"], "event") && !isempty(ft["cfg"]["event"])
        _extract_ft_events!(trigger_vec, trigger_info_vec, ft["cfg"]["event"], n_timepoints)
    end

    cols = Vector{AbstractVector}(undef, 4 + n_channels)
    cols[1] = time_vec
    cols[2] = collect(1:n_timepoints)
    cols[3] = trigger_vec
    cols[4] = trigger_info_vec
    _extract_channel_cols!(cols, 4, data_matrix, n_channels, n_timepoints)
    df = DataFrame(cols, vcat([:time, :sample, :trigger, :trigger_info], ch_names); copycols = false)

    continuous_data = ContinuousData(filepath, df, copy(layout), sample_rate, AnalysisInfo(sample_rate = Float64(sample_rate)))
    @info "Successfully converted to ContinuousData"
    return continuous_data
end

function _fieldtrip_to_epochdata(ft::AbstractDict, filepath::String, layout::Layout, sample_rate::Int)
    ch_names = _extract_ft_channel_names(vec(ft["label"]))
    trials = ft["trial"]
    n_trials = length(trials)
    n_channels = length(ch_names)
    n_timepoints = size(trials[1], 2)
    @info "Data: $n_channels channels × $n_timepoints timepoints × $n_trials trials"

    times = _extract_time_vector(ft["time"])
    has_ti = haskey(ft, "trialinfo") && !isempty(ft["trialinfo"])
    trial_triggers = has_ti ? vec(ft["trialinfo"]) : nothing

    epoch_dataframes = Vector{DataFrame}(undef, n_trials)
    time_zero_idx = _find_closest_time_index(times, 0.0)
    col_names = vcat([:time, :sample, :epoch, :trigger, :trigger_info], ch_names)
    sample_template = collect(1:n_timepoints)

    @inbounds for trial_idx = 1:n_trials
        cols = Vector{AbstractVector}(undef, 5 + n_channels)
        cols[1] = copy(times)
        cols[2] = copy(sample_template)
        cols[3] = fill(trial_idx, n_timepoints)
        trig_vec = zeros(Int, n_timepoints)
        info_vec = fill("", n_timepoints)
        if !isnothing(trial_triggers) && trial_idx <= length(trial_triggers)
            code = round(Int, trial_triggers[trial_idx])
            if code != 0
                trig_vec[time_zero_idx] = code
                info_vec[time_zero_idx] = string(code)
            end
        end
        cols[4] = trig_vec
        cols[5] = info_vec
        _extract_channel_cols!(cols, 5, trials[trial_idx], n_channels, n_timepoints)
        epoch_dataframes[trial_idx] = DataFrame(cols, copy(col_names); copycols = false)
    end

    epoch_data = EpochData(
        filepath,
        1,
        basename(filepath),
        epoch_dataframes,
        copy(layout),
        sample_rate,
        AnalysisInfo(sample_rate = Float64(sample_rate)),
    )
    @info "Successfully converted to EpochData with $n_trials epochs"
    return epoch_data
end

function _fieldtrip_to_erpdata(ft::AbstractDict, filepath::String, layout::Layout, sample_rate::Int)
    ch_names = _extract_ft_channel_names(vec(ft["label"]))
    data_matrix = ft["avg"]
    time_raw = vec(ft["time"])
    time_vec = time_raw isa Vector{Float64} ? copy(time_raw) : [Float64(x) for x in time_raw]

    n_channels = length(ch_names)
    n_timepoints = size(data_matrix, 2)
    @info "Data: $n_channels channels × $n_timepoints timepoints"

    cols = Vector{AbstractVector}(undef, 4 + n_channels)
    cols[1] = time_vec
    cols[2] = collect(1:n_timepoints)
    cols[3] = zeros(Int, n_timepoints)
    cols[4] = fill("", n_timepoints)
    _extract_channel_cols!(cols, 4, data_matrix, n_channels, n_timepoints)
    df = DataFrame(cols, vcat([:time, :sample, :trigger, :trigger_info], ch_names); copycols = false)

    erp_data = ErpData(filepath, 1, basename(filepath), df, copy(layout), sample_rate, AnalysisInfo(sample_rate = Float64(sample_rate)), 1)
    @info "Successfully converted to ErpData"
    return erp_data
end
