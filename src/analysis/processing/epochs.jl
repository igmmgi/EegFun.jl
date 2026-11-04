"""
    condition_parse_epoch(config::Dict) -> Vector{EpochCondition}
    condition_parse_epoch(toml_file::String) -> Vector{EpochCondition}

Parse epoch condition configurations from a TOML file path or a pre-parsed dictionary.
"""
condition_parse_epoch(toml_file::String) = condition_parse_epoch(TOML.parsefile(toml_file))

_expand_arg(x::AbstractVector) = isempty(x) ? x : mapreduce(_expand_arg_element, vcat, x)
_expand_arg(x::Union{Number,Symbol}) = [x]
_expand_arg(x::UnitRange) = collect(x)
_expand_arg(x) = x

_expand_arg_element(v::UnitRange) = collect(v)
_expand_arg_element(v) = [v]

const _TriggerSeqElem = Union{Int,Symbol,UnitRange{Int},Vector{Int}}

function allcomb(args...)
    expanded_args = map(_expand_arg, args)
    return vec([collect(_TriggerSeqElem, c) for c in Iterators.product(expanded_args...)])
end

function _parse_allcomb_string(condition_name::String, expr::String)
    m = match(r"^allcomb\((.*)\)$"s, strip(expr))
    isnothing(m) &&
        @minimal_error("Invalid allcomb expression for condition '$condition_name': expected 'allcomb(arg1, arg2, ...)'. Got: '$expr'")

    args = String[]
    depth, buf = 0, IOBuffer()
    for c in strip(m.captures[1])
        if c in ('[', '(')
            depth += 1
            write(buf, c)
        elseif c in (']', ')')
            depth -= 1
            write(buf, c)
        elseif c == ',' && depth == 0
            push!(args, strip(String(take!(buf))))
        else
            write(buf, c)
        end
    end
    rem = strip(String(take!(buf)))
    !isempty(rem) && push!(args, rem)
    isempty(args) && @minimal_error("allcomb() called with no arguments for condition '$condition_name'")

    parse_int(s) =
        (n = tryparse(Int, strip(s)); isnothing(n) && @minimal_error("Cannot parse '$s' as integer in allcomb for '$condition_name'"); n)
    parsed = map(args) do s
        s = strip(s)
        (m2 = match(r"^:([A-Za-z_]\w*)$", s); !isnothing(m2)) && return Symbol(m2.captures[1])
        (m2 = match(r"^(-?\d+):(-?\d+)$", s); !isnothing(m2)) && return parse(Int, m2.captures[1]):parse(Int, m2.captures[2])
        (m2 = match(r"^\[\[(.*)\]\]$", s); !isnothing(m2)) && return [parse_int.(split(m2.captures[1], ","))]
        (m2 = match(r"^\[(.*)\]$", s); !isnothing(m2)) && return parse_int.(split(m2.captures[1], ","))
        return parse_int(s)
    end
    return allcomb(parsed...)
end

function condition_parse_epoch(config::Dict)
    epochs_section = get(config, "epochs", Dict())
    conditions = EpochCondition[]

    for condition_config in get(epochs_section, "conditions", [])
        name = condition_config["name"]
        trigger_sequences_raw = get(condition_config, "trigger_sequences", nothing)
        isnothing(trigger_sequences_raw) && @minimal_error("trigger_sequences must be specified for condition '$name'")

        if trigger_sequences_raw isa String
            if startswith(strip(trigger_sequences_raw), "allcomb")
                trigger_sequences_raw = _parse_allcomb_string(name, trigger_sequences_raw)
            else
                @minimal_error(
                    "trigger_sequences for condition '$name' is a plain string, which is not supported. " *
                    "Use a TOML array-of-arrays (e.g. trigger_sequences = [[31, 51, 101], [39, 51, 101]]), " *
                    "or an allcomb() expression (e.g. trigger_sequences = \"allcomb([31, 39], 51:54, 101)\").",
                )
            end
        end

        trigger_sequences = Vector{Vector{_TriggerSeqElem}}()
        for sequence in trigger_sequences_raw
            seq_vec = sequence isa Vector ? sequence : [sequence]
            converted = _TriggerSeqElem[x isa String ? Symbol(x) : x isa Vector ? Vector{Int}(x) : x for x in seq_vec]
            push!(trigger_sequences, converted)
        end

        reference_index = get(condition_config, "reference_index", 1)
        timing_pairs_raw = get(condition_config, "timing_pairs", nothing)
        if isnothing(timing_pairs_raw)
            timing_pairs, min_interval, max_interval = nothing, nothing, nothing
        else
            timing_pairs = [(pair[1], pair[2]) for pair in timing_pairs_raw]
            min_interval = get(condition_config, "min_interval", nothing)
            max_interval = get(condition_config, "max_interval", nothing)
            (isnothing(min_interval) || isnothing(max_interval)) &&
                @minimal_error("Both min_interval and max_interval must be specified when timing_pairs is specified for condition '$name'")
        end

        mask_after_trigger = get(condition_config, "mask_after_trigger", nothing)
        mask_before_trigger = get(condition_config, "mask_before_trigger", nothing)
        mask_triggers = get(condition_config, "mask_triggers", nothing)
        mask_between_raw = get(condition_config, "mask_between_triggers", nothing)
        mask_between_triggers = isnothing(mask_between_raw) ? nothing : [(m[1], m[2]) for m in mask_between_raw]

        seq_length = length(trigger_sequences[1])
        all(seq -> length(seq) == seq_length, trigger_sequences) ||
            @minimal_error("All trigger sequences must have the same length for condition '$name'")
        (1 <= reference_index <= seq_length) || @minimal_error("reference_index must be between 1 and $seq_length for condition '$name'")

        if !isnothing(timing_pairs)
            min_interval < max_interval || @minimal_error("min_interval must be < max_interval for condition '$name'")
            for (start_idx, end_idx) in timing_pairs
                (1 <= start_idx <= seq_length && 1 <= end_idx <= seq_length) ||
                    @minimal_error("timing_pairs contains invalid indices for sequence of length $seq_length in condition '$name'")
                start_idx < end_idx || @minimal_error("timing_pairs must have start_idx < end_idx for condition '$name'")
            end
        end

        (!isnothing(mask_after_trigger) && !isnothing(mask_before_trigger)) &&
            @minimal_error("Cannot specify both 'mask_after_trigger' and 'mask_before_trigger' constraints for condition '$name'")

        push!(
            conditions,
            EpochCondition(
                name,
                trigger_sequences,
                reference_index,
                timing_pairs,
                min_interval,
                max_interval,
                mask_before_trigger,
                mask_after_trigger,
                mask_triggers,
                mask_between_triggers,
            ),
        )
    end

    return conditions
end

_interval_endpoints(interval::Tuple{Real,Real}) = (Float64(interval[1]), Float64(interval[2]))
_interval_endpoints(interval::AbstractVector{<:Real}) =
    length(interval) == 2 ? (Float64(interval[1]), Float64(interval[2])) :
    @minimal_error("Epoch interval must have exactly 2 elements, got $(length(interval))")
_interval_endpoints(interval::TimeSelection) = (interval.start, interval.stop)
_interval_endpoints(interval) = @minimal_error("Epoch interval must be a Tuple, Vector, or TimeSelection, got $(typeof(interval))")

function _validate_epoch_interval_params(dat::ContinuousData, epoch_interval)
    t_start, t_end = _interval_endpoints(epoch_interval)
    t_start <= t_end || @minimal_error("Epoch interval start ($t_start) must be less than or equal to end ($t_end)")
    hasproperty(dat.data, :trigger) || @minimal_error("Data must have a trigger column")
    hasproperty(dat.data, :time) || @minimal_error("Data must have a time column")
    !isempty(dat.data.time) || @minimal_error("Time column cannot be empty")
    !isempty(dat.data.trigger) || @minimal_error("Triggers column cannot be empty")
    issorted(dat.data.time) || @minimal_error("Time column must be sorted in ascending order")
    return nothing
end

function _mark_windows_at_indices!(
    dat::ContinuousData,
    reference_indices::Vector{Int},
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection},
    channel_out::Symbol,
)::Int
    times = dat.data.time
    col_out = dat.data[!, channel_out]
    n_samp = length(times)
    w_start_offset, w_end_offset = _interval_endpoints(time_window)

    n_marked = 0
    for idx in reference_indices
        if idx < 1 || idx > n_samp
            @minimal_warning "Reference index $idx is out of bounds, skipping"
            continue
        end
        ref_t = times[idx]
        s_idx = searchsortedfirst(times, ref_t + w_start_offset)
        e_idx = searchsortedlast(times, ref_t + w_end_offset)
        if s_idx <= e_idx
            @inbounds fill!(@view(col_out[s_idx:e_idx]), true)
            n_marked += (e_idx - s_idx + 1)
        end
    end
    return n_marked
end

function _init_bool_column!(df::DataFrame, col::Symbol)
    if hasproperty(df, col) && df[!, col] isa Vector{Bool}
        fill!(df[!, col], false)
    else
        df[!, col] = fill(false, nrow(df))
    end
    return nothing
end

"""
    mark_epoch_intervals!(dat::ContinuousData, triggers_of_interest::AbstractVector{<:Integer}, time_window; channel_out = :epoch_interval)
    mark_epoch_intervals!(dat::ContinuousData, channel_in::Symbol, time_window; channel_out = Symbol(string(channel_in) * "_window"))
    mark_epoch_intervals!(dat::ContinuousData, channels_in::AbstractVector{Symbol}, time_window; channel_out = :event_window)
    mark_epoch_intervals!(dat::ContinuousData, epoch_conditions, time_window; channel_out = :epoch_interval)

Mark samples within `time_window` around trigger events, boolean flag columns, or `EpochCondition` matches in-place.
"""
function mark_epoch_intervals!(
    dat::ContinuousData,
    triggers_of_interest::AbstractVector{<:Integer},
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection};
    channel_out::Symbol = :epoch_interval,
)
    _validate_epoch_interval_params(dat, time_window)
    _init_bool_column!(dat.data, channel_out)

    trig_col = dat.data.trigger
    all_reference_indices = Int[]
    for trigger in triggers_of_interest
        found = false
        @inbounds for i in eachindex(trig_col)
            if trig_col[i] == trigger
                push!(all_reference_indices, i)
                found = true
            end
        end
        found || @minimal_warning("Trigger $trigger not found in data")
    end

    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(all_reference_indices)) trigger locations"
    return dat
end

function mark_epoch_intervals!(
    dat::ContinuousData,
    channel_in::Symbol,
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection};
    channel_out::Symbol = Symbol(string(channel_in) * "_window"),
)
    _validate_epoch_interval_params(dat, time_window)
    hasproperty(dat.data, channel_in) || @minimal_error("Column $(channel_in) not found in data")
    eltype(dat.data[!, channel_in]) <: Bool || @minimal_error("Column $(channel_in) must be a Bool column")

    _init_bool_column!(dat.data, channel_out)
    reference_indices = findall(dat.data[!, channel_in])
    if isempty(reference_indices)
        @minimal_warning "No true values found in column $(channel_in)"
        return dat
    end

    n_marked = _mark_windows_at_indices!(dat, reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(reference_indices)) $(channel_in) events"
    return dat
end

function mark_epoch_intervals!(
    dat::ContinuousData,
    channels_in::AbstractVector{Symbol},
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection};
    channel_out::Symbol = :event_window,
)
    _validate_epoch_interval_params(dat, time_window)
    for channel_in in channels_in
        hasproperty(dat.data, channel_in) || @minimal_error("Column $(channel_in) not found in data")
        eltype(dat.data[!, channel_in]) <: Bool || @minimal_error("Column $(channel_in) must be a Bool column")
    end

    _init_bool_column!(dat.data, channel_out)
    reference_indices = Int[]
    for channel_in in channels_in
        append!(reference_indices, findall(dat.data[!, channel_in]))
    end
    sort!(unique!(reference_indices))

    if isempty(reference_indices)
        @minimal_warning "No true values found in any of the columns $(channels_in)"
        return dat
    end

    n_marked = _mark_windows_at_indices!(dat, reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(reference_indices)) events from $(length(channels_in)) columns"
    return dat
end

function _apply_trigger_mask(triggers::Vector, condition::EpochCondition)
    if isnothing(condition.mask_triggers) &&
       isnothing(condition.mask_between_triggers) &&
       isnothing(condition.mask_before_trigger) &&
       isnothing(condition.mask_after_trigger)
        return triggers, 1:length(triggers)
    end

    hidden = falses(length(triggers))
    if !isnothing(condition.mask_triggers)
        for t in condition.mask_triggers
            @inbounds for i in eachindex(triggers)
                triggers[i] == t && (hidden[i] = true)
            end
        end
    end

    if !isnothing(condition.mask_between_triggers)
        for (start_t, end_t) in condition.mask_between_triggers
            in_mask = false
            @inbounds for i in eachindex(triggers)
                if !in_mask
                    if triggers[i] == start_t
                        in_mask = true
                        hidden[i] = true
                    end
                elseif triggers[i] == end_t
                    in_mask = false
                else
                    hidden[i] = true
                end
            end
        end
    end

    if !isnothing(condition.mask_before_trigger)
        idx = findfirst(==(condition.mask_before_trigger), triggers)
        !isnothing(idx) && (hidden[1:idx] .= true)
    end

    if !isnothing(condition.mask_after_trigger)
        idx = findlast(==(condition.mask_after_trigger), triggers)
        !isnothing(idx) && (hidden[idx:end] .= true)
    end

    valid_indices = findall(.!hidden)
    return triggers[valid_indices], valid_indices
end

function _map_valid_indices(raw_matches::Vector{Vector{Int}}, valid_indices)
    (valid_indices isa UnitRange && first(valid_indices) == 1) && return raw_matches
    return [Int[valid_indices[i] for i in m] for m in raw_matches]
end

function _filter_matches(matches::Vector{Vector{Int}}, condition::EpochCondition, triggers, time_data)
    if !isnothing(condition.timing_pairs) && !isnothing(condition.min_interval) && !isnothing(condition.max_interval)
        min_int, max_int = condition.min_interval, condition.max_interval
        pairs = condition.timing_pairs
        matches = filter(matches) do m
            for (start_idx, end_idx) in pairs
                interval = time_data[m[end_idx]] - time_data[m[start_idx]]
                (interval < min_int || interval > max_int) && return false
            end
            return true
        end
    end
    return matches
end

function mark_epoch_intervals!(
    dat::ContinuousData,
    epoch_conditions::AbstractVector{<:EpochCondition},
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection};
    channel_out::Symbol = :epoch_interval,
)
    _validate_epoch_interval_params(dat, time_window)
    _init_bool_column!(dat.data, channel_out)

    all_reference_indices = Int[]
    n_trig = length(dat.data.trigger)
    for condition in epoch_conditions
        filtered_triggers, valid_indices = _apply_trigger_mask(dat.data.trigger, condition)
        matches = _map_valid_indices(search_sequence(filtered_triggers, condition.trigger_sequences), valid_indices)
        if isempty(matches)
            @minimal_warning "No triggers found for condition '$(condition.name)'"
            continue
        end

        matches = _filter_matches(matches, condition, dat.data.trigger, dat.data.time)
        if isempty(matches)
            @minimal_warning "No trigger sequences found that meet timing constraints for condition '$(condition.name)'"
            continue
        end

        for m in matches
            if condition.reference_index > length(m)
                @minimal_warning "Reference index $(condition.reference_index) exceeds sequence length $(length(m)) for condition '$(condition.name)'"
                continue
            end
            ref_idx = m[condition.reference_index]
            if ref_idx <= n_trig
                push!(all_reference_indices, ref_idx)
            else
                @minimal_warning "Reference index $(condition.reference_index) for condition '$(condition.name)' is out of bounds"
            end
        end
    end

    n_marked = _mark_windows_at_indices!(dat, all_reference_indices, time_window, channel_out)
    @info "Marked $n_marked samples across $(length(all_reference_indices)) reference points"
    return dat
end

mark_epoch_intervals!(
    dat::ContinuousData,
    epoch_condition::EpochCondition,
    time_window::Union{Tuple{Real,Real},AbstractVector{<:Real},TimeSelection};
    channel_out::Symbol = :epoch_interval,
) = mark_epoch_intervals!(dat, [epoch_condition], time_window; channel_out = channel_out)

function _shifted_time_slice(time_col::Vector{Float64}, pre::Int, win_len::Int, zero_idx::Int)::Vector{Float64}
    t0 = @inbounds time_col[zero_idx]
    t_vec = Vector{Float64}(undef, win_len)
    @inbounds @simd for i = 1:win_len
        t_vec[i] = time_col[pre+i-1] - t0
    end
    return t_vec
end

function _shifted_time_slice(time_col::AbstractVector, pre::Int, win_len::Int, zero_idx::Int)::Vector{Float64}
    t0 = Float64(time_col[zero_idx])
    t_vec = Vector{Float64}(undef, win_len)
    @inbounds for i = 1:win_len
        t_vec[i] = Float64(time_col[pre+i-1]) - t0
    end
    return t_vec
end

@inline function _slice_col(src, rng::UnitRange{Int})::AbstractVector
    src isa Vector{Float64} && return (src::Vector{Float64})[rng]
    src isa Vector{Int} && return (src::Vector{Int})[rng]
    src isa Vector{Bool} && return (src::Vector{Bool})[rng]
    return src[rng]
end

"""
    extract_epochs(dat::ContinuousData, condition::Integer, epoch_condition::EpochCondition, epoch_interval) -> EpochData
    extract_epochs(dat::ContinuousData, epoch_conditions::AbstractVector{<:EpochCondition}, epoch_interval) -> Vector{EpochData}
    extract_epochs(dat::ContinuousData, epoch_condition::EpochCondition, epoch_interval) -> Vector{EpochData}

Extract time-locked epochs from `dat` matching `epoch_condition` over `epoch_interval` `(start_time, end_time)` in seconds.
"""
function extract_epochs(dat::ContinuousData, condition::Integer, epoch_condition::EpochCondition, epoch_interval)
    start_time, end_time = _interval_endpoints(epoch_interval)
    start_time <= end_time || @minimal_error("Epoch interval start ($start_time) must be less than or equal to end ($end_time)")

    filtered_triggers, valid_indices = _apply_trigger_mask(dat.data.trigger, epoch_condition)
    matches = _map_valid_indices(search_sequence(filtered_triggers, epoch_condition.trigger_sequences), valid_indices)
    isempty(matches) && @minimal_error("None of the trigger sequences $(epoch_condition.trigger_sequences) found!")

    matches = _filter_matches(matches, epoch_condition, dat.data.trigger, dat.data.time)
    isempty(matches) && @minimal_error("No trigger sequences found that meet timing constraints for condition '$(epoch_condition.name)'")

    ref_pos = epoch_condition.reference_index
    matches = filter(m -> ref_pos <= length(m), matches)
    isempty(matches) && @minimal_error("Reference index $ref_pos exceeds sequence length for condition '$(epoch_condition.name)'")

    start_offset = round(Int, start_time * dat.sample_rate)
    end_offset = round(Int, end_time * dat.sample_rate)

    cols_to_keep = setdiff(propertynames(dat.data), [:condition, :condition_name, :file])
    meta_cols_kept = setdiff(meta_labels(dat), [:condition, :condition_name, :file])
    epoch_col_pos = clamp(length(meta_cols_kept) + 1, 1, length(cols_to_keep) + 1)

    n_kept = length(cols_to_keep)
    n_out_cols = n_kept + 1
    out_names = Vector{Symbol}(undef, n_out_cols)
    src_cols = Vector{AbstractVector}(undef, n_out_cols)
    for k = 1:(epoch_col_pos-1)
        out_names[k] = cols_to_keep[k]
        src_cols[k] = dat.data[!, columnindex(dat.data, cols_to_keep[k])]
    end
    out_names[epoch_col_pos] = :epoch
    src_cols[epoch_col_pos] = Int[]
    for k = epoch_col_pos:n_kept
        out_names[k+1] = cols_to_keep[k]
        src_cols[k+1] = dat.data[!, columnindex(dat.data, cols_to_keep[k])]
    end

    n_total = n_samples(dat)
    time_col = dat.data[!, :time]
    epochs = sizehint!(DataFrame[], length(matches))

    for (epoch, m) in enumerate(matches)
        zero_idx = m[ref_pos]
        pre = zero_idx + start_offset
        post = zero_idx + end_offset
        if pre < 1 || post > n_total
            @minimal_warning "Epoch $epoch extends beyond data bounds (pre=$pre, post=$post, data_length=$n_total) - skipping"
            continue
        end

        rng = pre:post
        win_len = post - pre + 1
        out_cols = Vector{AbstractVector}(undef, n_out_cols)
        for c = 1:n_out_cols
            col_sym = out_names[c]
            if c == epoch_col_pos
                out_cols[c] = fill(epoch, win_len)
            elseif col_sym == :time
                out_cols[c] = _shifted_time_slice(time_col, pre, win_len, zero_idx)
            else
                out_cols[c] = _slice_col(src_cols[c], rng)
            end
        end
        push!(epochs, DataFrame(out_cols, out_names; copycols = false))
    end

    return EpochData(dat.file, condition, epoch_condition.name, epochs, copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))
end

extract_epochs(dat::ContinuousData, epoch_conditions::AbstractVector{<:EpochCondition}, epoch_interval) =
    EpochData[extract_epochs(dat, idx, cond, epoch_interval) for (idx, cond) in enumerate(epoch_conditions)]

extract_epochs(dat::ContinuousData, epoch_condition::EpochCondition, epoch_interval) =
    extract_epochs(dat, [epoch_condition], epoch_interval)

function _accumulate_epoch_average!(
    avg_cols::Vector{Vector{Float64}},
    data::Vector{DataFrame},
    ch_col_indices::Vector{Int},
    n_timepoints::Int,
)
    n_ch = length(ch_col_indices)
    for epoch in data
        @inbounds for ch_idx = 1:n_ch
            dst = avg_cols[ch_idx]
            col = epoch[!, ch_col_indices[ch_idx]]
            if col isa Vector{Float64}
                col_f64 = col::Vector{Float64}
                @simd for i = 1:n_timepoints
                    dst[i] += col_f64[i]
                end
            else
                for i = 1:n_timepoints
                    dst[i] += col[i]
                end
            end
        end
    end
    return nothing
end

"""
    average_epochs(dat::EpochData) -> ErpData
    average_epochs(dat::AbstractVector{<:EpochData}) -> Vector{ErpData}
    average_epochs(file_pattern::String; input_dir = pwd(), participant_selection = participants(),
                   condition_selection = conditions(), output_dir = nothing)

Average epoched EEG data across trials to produce `ErpData` (in-memory or batch JLD2 processing).
"""
function average_epochs(dat::EpochData)
    isempty(dat.data) && @minimal_error("Cannot average empty EpochData")

    first_epoch = first(dat.data)
    numeric_columns = Symbol[col for col in propertynames(first_epoch) if (T = eltype(first_epoch[!, col]); T <: Number && T != Bool)]
    eeg_channels = intersect(numeric_columns, [channel_labels(dat); extra_labels(dat)])
    isempty(eeg_channels) && @minimal_error("No EEG channels found to average")

    n_timepoints = nrow(first_epoch)
    for (idx, epoch) in enumerate(dat.data)
        nrow(epoch) == n_timepoints ||
            @minimal_error("Epoch $idx has $(nrow(epoch)) timepoints, expected $n_timepoints. All epochs must have the same length.")
    end

    keep_time = :time in meta_labels(dat) && hasproperty(first_epoch, :time)
    n_eps = n_epochs(dat)
    n_chans = length(eeg_channels)
    avg_cols = Vector{Float64}[zeros(Float64, n_timepoints) for _ = 1:n_chans]
    ch_col_indices = Int[columnindex(first_epoch, ch) for ch in eeg_channels]

    _accumulate_epoch_average!(avg_cols, dat.data, ch_col_indices, n_timepoints)

    inv_n_eps = 1.0 / n_eps
    for avg_col in avg_cols
        @inbounds @simd for i = 1:n_timepoints
            avg_col[i] *= inv_n_eps
        end
    end

    if keep_time
        out_names = Vector{Symbol}(undef, n_chans + 1)
        out_cols = Vector{AbstractVector}(undef, n_chans + 1)
        out_names[1] = :time
        out_cols[1] = copy(first_epoch[!, :time])
        for i = 1:n_chans
            out_names[i+1] = eeg_channels[i]
            out_cols[i+1] = avg_cols[i]
        end
        erp = DataFrame(out_cols, out_names; copycols = false)
    else
        erp = DataFrame(AbstractVector[c for c in avg_cols], eeg_channels; copycols = false)
    end

    return ErpData(dat.file, dat.condition, dat.condition_name, erp, copy(dat.layout), dat.sample_rate, copy(dat.analysis_info), n_eps)
end

average_epochs(dat::AbstractVector{<:EpochData}) = average_epochs.(dat)

"""
    reject_epochs!(dat::EpochData, info::EpochRejectionInfo) -> EpochData
    reject_epochs!(dat::AbstractVector{<:EpochData}, info::AbstractVector{<:EpochRejectionInfo})
    reject_epochs!(dat::EpochData, bad_columns::Union{Symbol,AbstractVector{Symbol}}) -> EpochData
    reject_epochs!(dat::AbstractVector{<:EpochData}, bad_columns::Union{Symbol,AbstractVector{Symbol}})

Remove bad epochs from `dat` in-place using `EpochRejectionInfo` or one or more boolean flag columns.
"""
function reject_epochs!(dat::EpochData, info::EpochRejectionInfo)::EpochData
    if isempty(info.rejected)
        @info "No epochs to reject"
        return dat
    end

    n_eps = n_epochs(dat)
    rejected_indices = unique(Int[r.epoch for r in info.rejected])
    keep_mask = trues(n_eps)
    for ep in rejected_indices
        1 <= ep <= n_eps && (keep_mask[ep] = false)
    end
    dat.data = dat.data[keep_mask]
    @info "Condition $(dat.condition) ($(dat.condition_name)) - Rejected $(length(rejected_indices)) of $(info.info.n) epochs."
    return dat
end

function reject_epochs!(dat::AbstractVector{<:EpochData}, info::AbstractVector{<:EpochRejectionInfo})
    length(dat) == length(info) ||
        @minimal_error("Length of EpochData vector ($(length(dat))) must match EpochRejectionInfo vector ($(length(info)))")
    for (d, inf) in zip(dat, info)
        reject_epochs!(d, inf)
    end
    return dat
end

function reject_epochs!(dat::EpochData, bad_columns::AbstractVector{Symbol})::EpochData
    isempty(dat.data) && @minimal_error("Cannot remove bad epochs from empty EpochData")
    isempty(bad_columns) && return dat

    first_epoch = first(dat.data)
    col_indices = Vector{Int}(undef, length(bad_columns))
    for (i, col) in enumerate(bad_columns)
        hasproperty(first_epoch, col) || @minimal_error("Column '$col' not found in epoch data")
        col_indices[i] = columnindex(first_epoch, col)
        col_type = eltype(first_epoch[!, col_indices[i]])
        col_type != Bool && @minimal_warning("Column '$col' is not Boolean type ($(col_type)). Non-zero values will be treated as 'true'")
    end

    n_eps = n_epochs(dat)
    good_epochs = sizehint!(DataFrame[], n_eps)
    n_removed = 0
    for epoch_df in dat.data
        has_bad = false
        for col_idx in col_indices
            col_vals = epoch_df[!, col_idx]
            if eltype(col_vals) == Bool ? any(col_vals) : any(!iszero, col_vals)
                has_bad = true
                break
            end
        end
        has_bad ? (n_removed += 1) : push!(good_epochs, epoch_df)
    end

    if n_removed > 0
        @info "Condition $(dat.condition) ($(dat.condition_name)) removed $n_removed of $n_eps epochs ($(round(100*n_removed/n_eps, digits=1))%)"
    end
    dat.data = good_epochs
    return dat
end

reject_epochs!(dat::EpochData, bad_column::Symbol)::EpochData = reject_epochs!(dat, [bad_column])

function reject_epochs!(dat::AbstractVector{<:EpochData}, bad_column::Symbol)
    for d in dat
        reject_epochs!(d, bad_column)
    end
    return dat
end

function reject_epochs!(dat::AbstractVector{<:EpochData}, bad_columns::AbstractVector{Symbol})
    for d in dat
        reject_epochs!(d, bad_columns)
    end
    return dat
end

"""
    reject_epochs(dat::EpochData, info::EpochRejectionInfo) -> EpochData
    reject_epochs(dat::AbstractVector{<:EpochData}, info::AbstractVector{<:EpochRejectionInfo}) -> Vector{EpochData}
    reject_epochs(dat::EpochData, bad_columns::Union{Symbol,AbstractVector{Symbol}}) -> EpochData
    reject_epochs(dat::AbstractVector{<:EpochData}, bad_columns::Union{Symbol,AbstractVector{Symbol}}) -> Vector{EpochData}

Non-mutating version of [`reject_epochs!`](@ref).
"""
reject_epochs(dat::EpochData, info::EpochRejectionInfo)::EpochData = reject_epochs!(copy(dat), info)

function reject_epochs(dat::AbstractVector{<:EpochData}, info::AbstractVector{<:EpochRejectionInfo})
    length(dat) == length(info) ||
        @minimal_error("Length of EpochData vector ($(length(dat))) must match EpochRejectionInfo vector ($(length(info)))")
    return reject_epochs.(dat, info)
end

reject_epochs(dat::EpochData, bad_columns::AbstractVector{Symbol})::EpochData = reject_epochs!(copy(dat), bad_columns)
reject_epochs(dat::EpochData, bad_column::Symbol) = reject_epochs(dat, [bad_column])
reject_epochs(dat::AbstractVector{<:EpochData}, bad_column::Symbol) = [reject_epochs(d, bad_column) for d in dat]
reject_epochs(dat::AbstractVector{<:EpochData}, bad_columns::AbstractVector{Symbol}) = [reject_epochs(d, bad_columns) for d in dat]

"""
    epochs_count(dat::ContinuousData, toml_file::String; print_table::Bool = true) -> DataFrame
    epochs_count(dat::ContinuousData, conditions::Vector{EpochCondition}; print_table::Bool = true) -> DataFrame

Dry-run epoch counting for each condition without extracting waveform slices.
"""
epochs_count(dat::ContinuousData, toml_file::String; print_table::Bool = true) =
    epochs_count(dat, condition_parse_epoch(TOML.parsefile(toml_file)); print_table = print_table)

function epochs_count(dat::ContinuousData, conditions::Vector{EpochCondition}; print_table::Bool = true)
    n_cond = length(conditions)
    counts = Vector{Int}(undef, n_cond)
    names = Vector{String}(undef, n_cond)

    for (idx, condition) in enumerate(conditions)
        names[idx] = condition.name
        filtered_triggers, valid_indices = _apply_trigger_mask(dat.data.trigger, condition)
        matches = _map_valid_indices(search_sequence(filtered_triggers, condition.trigger_sequences), valid_indices)
        matches = _filter_matches(matches, condition, dat.data.trigger, dat.data.time)
        counts[idx] = count(m -> condition.reference_index <= length(m), matches)
    end

    df = DataFrame(condition = collect(1:n_cond), condition_name = names, n_epochs = counts)
    print_table && pretty_table(stdout, df; alignment = [:l, :l, :r], title = "Dry Run: Epoch Counts")
    return df
end

function _build_base_epochs_df(epochs::AbstractVector{<:EpochData})::DataFrame
    return DataFrame(
        file = [filename(epoch) for epoch in epochs],
        condition = [condition_number(epoch) for epoch in epochs],
        condition_name = [condition_name(epoch) for epoch in epochs],
    )
end

"""
    epochs_table(epochs::AbstractVector{<:EpochData}; print_table::Bool = true, kwargs...) -> DataFrame
    epochs_table(epochs_raw::AbstractVector{<:EpochData}, epochs_corrected::AbstractVector{<:EpochData};
                 print_table::Bool = true, kwargs...) -> DataFrame

Format and optionally print a summary (or before/after rejection comparison) table of epoch counts per condition.
"""
function epochs_table(epochs::AbstractVector{<:EpochData}; print_table::Bool = true, kwargs...)
    isempty(epochs) && @minimal_error("epochs vector cannot be empty")
    df = _build_base_epochs_df(epochs)
    df.n_epochs = [n_epochs(epoch) for epoch in epochs]
    print_table && pretty_table(stdout, df; alignment = [:l, :r, :l, :r], kwargs...)
    return df
end

function epochs_table(
    epochs_raw::AbstractVector{<:EpochData},
    epochs_corrected::AbstractVector{<:EpochData};
    print_table::Bool = true,
    kwargs...,
)
    length(epochs_raw) != length(epochs_corrected) && @minimal_error("epochs_raw and epochs_corrected must have same length")
    df = _build_base_epochs_df(epochs_raw)
    df.n_epochs_raw = [n_epochs(epoch) for epoch in epochs_raw]
    df.n_epochs_corrected = [n_epochs(epoch) for epoch in epochs_corrected]
    df.percentage = round.((df.n_epochs_corrected ./ df.n_epochs_raw) .* 100; digits = 1)
    print_table && pretty_table(stdout, df; alignment = [:l, :r, :l, :r, :r, :r], kwargs...)
    return df
end

"""
    log_epochs_table(epochs...; print_table::Bool = false, kwargs...) -> DataFrame

Log an [`epochs_table`](@ref) via `@info` and return the summary `DataFrame`.
"""
function log_epochs_table(epochs...; print_table::Bool = false, kwargs...)
    df = epochs_table(epochs...; print_table = print_table)
    table_output = sprint() do output_io
        io_context = IOContext(output_io, :displaysize => displaysize(stdout))
        pretty_table(io_context, df; alignment = [:l, :r, :l, :r, :r, :r], kwargs...)
    end
    @info "\n\n$table_output\n"
    return df
end

_validate_epochs_pattern(pattern::String) =
    contains(pattern, "epochs") ? nothing : "average_epochs only works with epoch data. File pattern must contain 'epochs', got: '$pattern'"

_default_average_output_dir(input_dir::String, pattern::String) = joinpath(input_dir, "averaged_$(clean_pattern(pattern))")

function _process_average_file(filepath::String, output_path::String, condition_selection::Function)
    filename = basename(filepath)
    epochs_data = read_data(filepath)
    isnothing(epochs_data) && return BatchResult(false, filename, "No data variables found")
    epochs_data isa AbstractVector{<:EpochData} || return BatchResult(false, filename, "Invalid data type: expected Vector{EpochData}")

    epochs_data = _condition_select(epochs_data, condition_selection)
    erps_data = average_epochs.(epochs_data)
    jldsave(output_path; data = erps_data, format_version = EEGFUN_FORMAT_VERSION)
    return BatchResult(true, filename, "Averaged $(length(erps_data)) condition(s)")
end

function average_epochs(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection = participants(),
    condition_selection = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    _with_batch_logging("epoch averaging", "average_epochs.log", input_dir, output_dir) do ctx
        @log_call "average_epochs"
        error_msg = _validate_epochs_pattern(file_pattern)
        !isnothing(error_msg) && @minimal_error(error_msg)
        output_dir = _prepare_output_dir!(ctx, _default_average_output_dir(input_dir, file_pattern))
        process_fn = (input_path, output_path) -> _process_average_file(input_path, output_path, condition_selection)
        batch_process(
            process_fn,
            file_pattern,
            input_dir,
            output_dir,
            participant_selection,
            "Averaging";
            filename_modifier = fn -> replace(fn, "epochs" => "erps"),
        )
    end
end
