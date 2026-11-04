function _clean_triggers(trigger_data::Vector{<:Integer})::Vector{<:Integer}
    n = length(trigger_data)
    cleaned = zeros(eltype(trigger_data), n)
    n == 0 && return cleaned

    @inbounds begin
        trigger_data[1] > 0 && (cleaned[1] = trigger_data[1])
        for i = 2:n
            v = trigger_data[i]
            (v > 0 && v != trigger_data[i-1]) && (cleaned[i] = v)
        end
    end
    return cleaned
end

"""
    trigger_count(dat::ContinuousData) -> TriggerInfo
    trigger_count(df::DataFrame) -> TriggerInfo
    trigger_count(dat::BiosemiDataFormat.BiosemiData) -> TriggerInfo
    trigger_count(dat::EuropeanDataFormat.EdfData) -> TriggerInfo
    trigger_count(dat::BrainVisionDataFormat.BrainVisionData) -> TriggerInfo
    trigger_count(dat::ExtensibleDataFormat.XdfData) -> TriggerInfo

Count occurrences of each non-zero trigger value in continuous or raw EEG data,
returning a formatted [`TriggerInfo`](@ref) summary table.
"""
function trigger_count(df::DataFrame)::TriggerInfo
    hasproperty(df, :trigger) || @minimal_error("DataFrame must have a trigger column")
    t_info = hasproperty(df, :trigger_info) ? df.trigger_info : nothing
    return _trigger_count_impl([df.trigger], ["count"]; trigger_info = t_info)
end
trigger_count(dat::ContinuousData)::TriggerInfo = trigger_count(dat.data)

function _trigger_count_impl(
    trigger_datasets::Vector{<:Vector{<:Integer}},
    column_names::Vector{String};
    trigger_info::Union{Nothing,Vector{String}} = nothing,
)
    n_ds = length(trigger_datasets)
    counts_per_ds = [Dict{Int,Int}() for _ = 1:n_ds]

    for (ds_idx, dataset) in enumerate(trigger_datasets)
        counts = counts_per_ds[ds_idx]
        @inbounds for val in dataset
            if val != 0
                iv = Int(val)
                counts[iv] = get(counts, iv, 0) + 1
            end
        end
    end

    non_zero_triggers = Int[]
    if n_ds == 1
        append!(non_zero_triggers, keys(counts_per_ds[1]))
    elseif n_ds > 1
        seen = Set{Int}()
        for counts in counts_per_ds, k in keys(counts)
            k in seen || (push!(seen, k); push!(non_zero_triggers, k))
        end
    end
    sort!(non_zero_triggers)

    has_info = !isnothing(trigger_info)
    n_out_cols = 1 + (has_info ? 1 : 0) + length(column_names)
    out_names = Vector{Symbol}(undef, n_out_cols)
    out_cols = Vector{AbstractVector}(undef, n_out_cols)
    out_names[1] = :trigger
    out_cols[1] = non_zero_triggers

    col_offset = 1
    if has_info
        col_offset += 1
        out_names[col_offset] = :trigger_info
        if isempty(non_zero_triggers)
            out_cols[col_offset] = String[]
        else
            trigger_info_map = Dict{Int,String}()
            ds1 = trigger_datasets[1]
            tinfo = trigger_info::Vector{String}
            @inbounds for i in eachindex(ds1, tinfo)
                t = Int(ds1[i])
                if t != 0 && !haskey(trigger_info_map, t)
                    trigger_info_map[t] = tinfo[i]
                end
            end
            out_cols[col_offset] = String[get(trigger_info_map, t, "") for t in non_zero_triggers]
        end
    end

    for (i, cname) in enumerate(column_names)
        out_names[col_offset+i] = Symbol(cname)
        out_cols[col_offset+i] = isempty(non_zero_triggers) ? Int[] : Int[get(counts_per_ds[i], t, 0) for t in non_zero_triggers]
    end

    return TriggerInfo(DataFrame(out_cols, out_names; copycols = false))
end

function Base.show(io::IO, info::TriggerInfo)
    if isempty(info.data)
        println(io, "No non-zero triggers found in the data.")
        return
    end
    col_names = names(info.data)
    title = ("raw_count" in col_names && "cleaned_count" in col_names) ? "Trigger Count Summary (Raw vs Cleaned)" : "Trigger Count Summary"
    n_cols = length(col_names)
    alignment = hasproperty(info.data, :trigger_info) ? [:r, :l, fill(:r, n_cols - 2)...] : fill(:r, n_cols)
    pretty_table(io, info.data; title = title, alignment = alignment)
    println(io)
end

function trigger_count(dat::BiosemiDataFormat.BiosemiData)::TriggerInfo
    cleaned_triggers = _clean_triggers(dat.triggers.raw)
    return _trigger_count_impl([dat.triggers.raw, cleaned_triggers], ["raw_count", "cleaned_count"])
end

function trigger_count(dat::EuropeanDataFormat.EdfData)::TriggerInfo
    n_samp = size(dat.data, 1)
    sr = dat.header.sample_rate[1]
    trigger = zeros(Int, n_samp)
    t_info = fill("", n_samp)

    if !isnothing(dat.triggers) && !isempty(dat.triggers.onset)
        unique_annotations = sort(unique(dat.triggers.annotation))
        value_to_trigger = Dict(ann => i for (i, ann) in enumerate(unique_annotations))
        for i in eachindex(dat.triggers.onset)
            sample_idx = round(Int, dat.triggers.onset[i] * sr) + 1
            if 1 <= sample_idx <= n_samp
                ann = dat.triggers.annotation[i]
                trigger[sample_idx] = value_to_trigger[ann]
                t_info[sample_idx] = ann
            end
        end
    end
    return _trigger_count_impl([trigger], ["count"]; trigger_info = t_info)
end

function trigger_count(dat::BrainVisionDataFormat.BrainVisionData)::TriggerInfo
    (isnothing(dat.markers) || isempty(dat.markers)) && return _trigger_count_impl([Int[]], ["count"])
    trigger, t_info = _extract_triggers_from_markers(dat.markers, size(dat.data, 1))
    return _trigger_count_impl([trigger], ["count"]; trigger_info = t_info)
end

function trigger_count(dat::ExtensibleDataFormat.XdfData)::TriggerInfo
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    isempty(eeg_streams) && return _trigger_count_impl([Int[]], ["count"])
    eeg_stream = eeg_streams[1]
    n_samp = size(eeg_stream.time_series, 1)
    trigger = zeros(Int, n_samp)
    t_info = fill("", n_samp)

    marker_streams = [s for s in values(dat.streams) if s.header.type == "Markers"]
    if !isempty(marker_streams)
        marker_stream = marker_streams[1]
        for (i, t) in enumerate(marker_stream.timestamps)
            idx = find_closest_time_index(eeg_stream.timestamps, t)
            marker_val = marker_stream.time_series[i, 1]
            marker_str = string(marker_val)
            int_val = tryparse(Int, marker_str)
            trigger[idx] = !isnothing(int_val) ? int_val : 1
            t_info[idx] = marker_str
        end
    end
    return _trigger_count_impl([trigger], ["count"]; trigger_info = t_info)
end

"""
    search_sequence(array, sequences::Vector{<:Vector}; ignore_values = [0], sort_indices = true) -> Vector{Vector{Int}}
    search_sequence(array, sequence::Vector; ignore_values = [0], sort_indices = true) -> Vector{Vector{Int}}
    search_sequence(array, value::Integer; ignore_values = [0]) -> Vector{Int}
    search_sequence(array, range::UnitRange; ignore_values = [0], sort_indices = true) -> Vector{Int}
    search_sequence(array, ranges::Vector{<:UnitRange}; ignore_values = [0], sort_indices = true) -> Vector{Int}

Search `array` for trigger onsets, ranges, or multi-element trigger sequences (supporting `:any` wildcards,
`UnitRange`s, and set-match `Vector{Int}` elements).
"""
function search_sequence(array, sequences::Vector{<:Vector}; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)
    (isempty(array) || isempty(sequences)) && return Vector{Int}[]
    length(sequences) == 1 &&
        return _search_single_sequence(array, sequences[1]; ignore_values = ignore_values, sort_indices = sort_indices)

    all_matches = Vector{Int}[]
    for sequence in sequences
        append!(all_matches, _search_single_sequence(array, sequence; ignore_values = ignore_values, sort_indices = false))
    end
    unique!(m -> m[1], all_matches)
    sort_indices && sort!(all_matches, by = first)
    return all_matches
end

function _search_single_sequence(array, sequence::Vector; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)
    (isempty(array) || isempty(sequence)) && return Vector{Int}[]

    if all(x -> x isa UnitRange, sequence)
        indices = search_sequence(array, UnitRange{Int}[x for x in sequence]; ignore_values = ignore_values, sort_indices = sort_indices)
        return Vector{Int}[[idx] for idx in indices]
    end

    seq_len = length(sequence)
    if seq_len == 1
        indices = _search_single_trigger(array, sequence[1], ignore_values)
        return Vector{Int}[[idx] for idx in indices]
    end

    idx_start_positions = _search_single_trigger(array, sequence[1], ignore_values)
    matches = Vector{Int}[]
    max_idx = length(array) - seq_len + 1
    for idx in idx_start_positions
        idx > max_idx && continue
        positions = _collect_sequence_positions(array, sequence, idx, ignore_values)
        !isnothing(positions) && push!(matches, positions)
    end
    sort_indices && !issorted(matches, by = first) && sort!(matches, by = first)
    return matches
end

search_sequence(array, sequence::Vector; ignore_values::Vector{Int} = [0], sort_indices::Bool = true) =
    _search_single_sequence(array, sequence; ignore_values = ignore_values, sort_indices = sort_indices)

function search_sequence(array, value::Integer; ignore_values::Vector{Int} = [0])
    indices = Int[]
    (isempty(array) || value in ignore_values) && return indices
    @inbounds begin
        array[1] == value && push!(indices, 1)
        for i = 2:length(array)
            (array[i] == value && array[i-1] != value) && push!(indices, i)
        end
    end
    return indices
end

search_sequence(array, range::UnitRange; ignore_values::Vector{Int} = [0], sort_indices::Bool = true) =
    search_sequence(array, [range]; ignore_values = ignore_values, sort_indices = sort_indices)

function search_sequence(array, ranges::Vector{<:UnitRange}; ignore_values::Vector{Int} = [0], sort_indices::Bool = true)
    isempty(ranges) && return Int[]
    all_indices = Int[]
    for value in union(ranges...)
        append!(all_indices, search_sequence(array, value; ignore_values = ignore_values))
    end
    sort_indices && sort!(all_indices)
    return all_indices
end

_search_single_trigger(array, trigger::Integer, ignore_values::Vector{Int} = [0]) =
    search_sequence(array, trigger; ignore_values = ignore_values)
_search_single_trigger(array, trigger::UnitRange, ignore_values::Vector{Int} = [0]) = search_sequence(array, [trigger])
_search_single_trigger(array, trigger::Symbol, ignore_values::Vector{Int} = [0]) = @minimal_error("Single wildcard sequences not supported")
function _search_single_trigger(array, trigger::Vector{<:Integer}, ignore_values::Vector{Int} = [0])
    isempty(trigger) && return Int[]
    out = Int[]
    for v in trigger
        append!(out, search_sequence(array, Int(v); ignore_values = ignore_values))
    end
    sort!(out)
    return out
end

@inline _is_ignored(v, ignore_values::Vector{Int}, ignore_zero_only::Bool) = ignore_zero_only ? (v == 0) : (v in ignore_values)

function _collect_sequence_positions(array, sequence, start_idx::Int, ignore_values::Vector{Int})
    n = length(array)
    seq_len = length(sequence)
    ignore_zero_only = length(ignore_values) == 1 && @inbounds(ignore_values[1]) == 0

    # Check step 2 before allocating output vector so failed candidate starts allocate 0 bytes
    current_idx = start_idx + 1
    @inbounds while current_idx <= n && _is_ignored(array[current_idx], ignore_values, ignore_zero_only)
        current_idx += 1
    end
    current_idx > n && return nothing
    _matches_expected(@inbounds(array[current_idx]), @inbounds(sequence[2])) || return nothing

    positions = Vector{Int}(undef, seq_len)
    @inbounds positions[1] = start_idx
    @inbounds positions[2] = current_idx

    @inbounds for step = 3:seq_len
        current_idx += 1
        while current_idx <= n && _is_ignored(array[current_idx], ignore_values, ignore_zero_only)
            current_idx += 1
        end
        current_idx > n && return nothing
        _matches_expected(array[current_idx], sequence[step]) || return nothing
        positions[step] = current_idx
    end
    return positions
end

_matches_expected(actual::Integer, expected::Integer) = actual == expected
_matches_expected(actual::Real, expected::Integer) = actual == expected
_matches_expected(actual::Integer, expected::Real) = actual == expected
_matches_expected(actual::Real, expected::Real) = actual == expected
_matches_expected(actual::Real, expected::Symbol) = expected == :any
_matches_expected(actual::Real, expected::UnitRange) = actual in expected
_matches_expected(actual::Real, expected::Vector{<:Integer}) = actual in expected
_matches_expected(actual, expected) = @minimal_error("Unsupported sequence type: $expected")

const _TriggerSeqEntry = NamedTuple{(:sequence, :t0),Tuple{Vector{Int},Int}}

"""
    trigger_info(dat::EpochData) -> Vector{NamedTuple{(:sequence, :t0),Tuple{Vector{Int},Int}}}
    trigger_info(dat_vec::AbstractVector{<:EpochData}) -> Dict{String,Vector{...}}
    trigger_info(file_pattern::String; input_dir = pwd(), participant_selection = participants())

Return the unique non-zero trigger sequences and `t == 0` trigger values observed per condition in epoched data.
"""
function trigger_info(dat::EpochData)
    unique_seqs = _TriggerSeqEntry[]
    for epoch in dat
        if hasproperty(epoch, :trigger) && hasproperty(epoch, :time)
            zero_idx = find_closest_time_index(epoch.time, 0.0)
            t0_trigger = epoch.trigger[zero_idx]
            seq = Int[x for x in epoch.trigger if x != 0]
            nt = (sequence = seq, t0 = t0_trigger)
            nt in unique_seqs || push!(unique_seqs, nt)
        end
    end
    return unique_seqs
end

function trigger_info(dat_vec::AbstractVector{<:EpochData})
    info = Dict{String,Vector{_TriggerSeqEntry}}()
    for dat in dat_vec
        info[dat.condition_name] = trigger_info(dat)
    end
    return info
end

function trigger_info(file_pattern::String; input_dir::String = pwd(), participant_selection::Function = participants())
    if endswith(file_pattern, ".jld2")
        direct_path = isfile(file_pattern) ? file_pattern : joinpath(input_dir, file_pattern)
        if isfile(direct_path)
            dat = read_data(direct_path)
            return isnothing(dat) ? nothing : trigger_info(dat isa Vector ? dat : [dat])
        end
    end

    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    if isempty(files)
        @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
        return nothing
    end

    @info "Scanning $(length(files)) files for unique trigger sequences..."
    merged_info = Dict{String,Vector{_TriggerSeqEntry}}()
    for file in files
        dat = read_data(joinpath(input_dir, file))
        isnothing(dat) && continue
        info = trigger_info(dat isa EpochData ? [dat] : dat)
        for (cond, seqs) in info
            if !haskey(merged_info, cond)
                merged_info[cond] = seqs
            else
                existing = merged_info[cond]
                for seq in seqs
                    seq in existing || push!(existing, seq)
                end
            end
        end
    end
    return merged_info
end
