# Channel, epoch, condition, participant, and time-interval selection predicates and `subset`.

######################################
######### Types & Constants ##########
######################################

struct _AllSamples <: Function end

struct _ChannelNamesSelector <: Function
    channel_names::Vector{Symbol}
end

struct _ChannelNumbersSelector <: Function
    channel_numbers::Vector{Int}
end

struct _SampleIndicesSelector <: Function
    indices::Vector{Int}
end

struct _EpochNumbersSelector <: Function
    epoch_numbers::Vector{Int}
end

######################################
########### User Functions ###########
######################################

"""
    subset(dat::EegData; channel_selection = channels(), sample_selection = samples(), interval_selection = times(), epoch_selection = epochs(), include_extra = false)
    subset(datasets::AbstractVector{<:EegData}; condition_selection = conditions(), kwargs...)

Return a subsetted copy of `dat` (or vector of condition datasets) filtered by channel, sample/time interval, epoch, and/or condition predicates.
"""
function subset(
    dat::T;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::T where {T<:SingleDataFrameEeg}
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, combined_sel, include_extra)
    dat_subset = _subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat, dat_subset, layout_subset)
end

function subset(
    dat::TimeFreqData;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::TimeFreqData
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, combined_sel, include_extra)
    power_subset = _subset_dataframe(dat.data_power, selected_channels, selected_samples)
    phase_subset = _subset_dataframe(dat.data_phase, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

function subset(
    dat::T;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::T where {T<:MultiDataFrameEeg}
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    dat_subset = _subset_dataframes(dat.data, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, dat_subset, layout_subset)
end

function subset(
    dat::TimeFreqEpochData;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::TimeFreqEpochData
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    power_subset = _subset_dataframes(dat.data_power, selected_epochs, selected_channels, selected_samples)
    phase_subset = _subset_dataframes(dat.data_phase, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

function subset(
    datasets::AbstractVector{T};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Vector{T} where {T<:Union{ErpData,TimeFreqData,SpectrumData}}
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    return T[
        subset(
            datasets[i];
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            include_extra = include_extra,
        ) for i in selected_conditions
    ]
end

function subset(
    datasets::AbstractVector{T};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::Vector{T} where {T<:Union{EpochData,TimeFreqEpochData}}
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    return T[
        subset(
            datasets[i];
            channel_selection = channel_selection,
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            epoch_selection = epoch_selection,
            include_extra = include_extra,
        ) for i in selected_conditions
    ]
end

"""
    channels()
    channels(names::Symbol...)
    channels(names::AbstractVector{Symbol})
    channels(indices)
    channels(predicate::Function)

Create a channel selection predicate for `channel_selection` arguments.
"""
channels() = x -> fill(true, length(x))

channels(channel_names::AbstractVector{Symbol}) = _ChannelNamesSelector(collect(Symbol, channel_names))

channels(channel_name::Symbol) = channels([channel_name])

channels(channel_names::Symbol...) = channels(collect(Symbol, channel_names))

channels(channel_names::AbstractString...) = channels(Symbol.(collect(channel_names)))

channels(channel_number::Integer) = channels([Int(channel_number)])

channels(channel_numbers::AbstractVector{<:Integer}) = _ChannelNumbersSelector(collect(Int, channel_numbers))

channels(channel_numbers::AbstractRange{<:Integer}) = _ChannelNumbersSelector(collect(Int, channel_numbers))

channels(channel_ranges::AbstractVector{<:AbstractRange{<:Integer}}) = _ChannelNumbersSelector(collect(Int, union(channel_ranges...)))

channels(rx::Regex) = x -> Bool[occursin(rx, string(ch)) for ch in x]

channels(predicate::Function) = predicate

channels(channel_names::AbstractVector{<:AbstractString}) = channels(Symbol.(channel_names))

channels(channel_name::AbstractString) = channels([Symbol(channel_name)])

"""
    channels_not(args...)

Create a channel selection predicate excluding the specified channels (inverse of `channels`).
"""
channels_not(args...) = _negate_selector(channels(args...))

"""
    components()
    components(indices...)
    components(predicate::Function)

Create an ICA component selection predicate.
"""
components() = x -> fill(true, length(x))

components(component_numbers::AbstractVector{<:Integer}) =
    let nums = collect(Int, component_numbers)
        x -> Bool[i in nums for i = 1:length(x)]
    end

components(component_numbers::AbstractRange{<:Integer}) = components(collect(Int, component_numbers))

components(component_number::Integer) = components([Int(component_number)])

components(component_ranges::AbstractVector{<:AbstractRange{<:Integer}}) = components(collect(Int, union(component_ranges...)))

components(component_numbers::Integer...) = components(collect(Int, component_numbers))

components(predicate::Function) = predicate

"""
    components_not(args...)

Create an ICA component selection predicate excluding the specified components (inverse of `components`).
"""
components_not(args...) = _negate_selector(components(args...))

"""
    samples()
    samples(column::Symbol)
    samples(predicate::Function)

Create a row/sample selection predicate for `sample_selection` arguments.
"""
samples() = _AllSamples()

samples(column::Symbol) = x -> x[!, column]

samples(predicate::Function) = predicate

"""
    samples_or(columns::Vector{Symbol})

Select rows where any boolean column in `columns` is `true`.
"""
samples_or(columns::Vector{Symbol}) = x -> _combine_sample_cols(x, columns, true, false)

"""
    samples_and(columns::Vector{Symbol})

Select rows where all boolean columns in `columns` are `true`.
"""
samples_and(columns::Vector{Symbol}) = x -> _combine_sample_cols(x, columns, false, false)

"""
    samples_not(column::Symbol)

Select rows where boolean `column` is `false`.
"""
samples_not(column::Symbol) = x -> .!(x[!, column])

"""
    samples_or_not(columns::Vector{Symbol})

Select rows where none of the boolean `columns` are `true`.
"""
samples_or_not(columns::Vector{Symbol}) = x -> _combine_sample_cols(x, columns, true, true)

"""
    samples_and_not(columns::Vector{Symbol})

Select rows where not all boolean `columns` are `true`.
"""
samples_and_not(columns::Vector{Symbol}) = x -> _combine_sample_cols(x, columns, false, true)

"""
    times()
    times(t::Real)
    times(start::Real, stop::Real)
    times(interval::Tuple{Real,Real})

Construct a time-interval specification (`AllSelection` or `TimeSelection`) for `interval_selection`.
"""
times() = AllSelection()

times(t::Real) = TimeSelection(Float64(t), Float64(t))

function times(start::Real, stop::Real)
    s, e = Float64(start), Float64(stop)
    s <= e || @minimal_error("Time interval start ($s) must be <= stop ($e)")
    return TimeSelection(s, e)
end

times(interval::Tuple{Real,Real}) = times(interval[1], interval[2])

"""
    epochs()
    epochs(indices...)
    epochs(predicate::Function)

Create an epoch selection predicate for `epoch_selection` arguments.
"""
epochs() = x -> fill(true, length(x))

epochs(epoch_numbers::AbstractVector{<:Integer}) = _EpochNumbersSelector(collect(Int, epoch_numbers))

epochs(epoch_numbers::AbstractRange{<:Integer}) = epochs(collect(Int, epoch_numbers))

epochs(epoch_number::Integer) = epochs([Int(epoch_number)])

epochs(epoch_ranges::AbstractVector{<:AbstractRange{<:Integer}}) = epochs(collect(Int, union(epoch_ranges...)))

epochs(epoch_numbers::Integer...) = epochs(collect(Int, epoch_numbers))

epochs(predicate::Function) = x -> map(predicate, x)

"""
    epochs_not(args...)

Create an epoch selection predicate excluding the specified epochs (inverse of `epochs`).
"""
epochs_not(args...) = _negate_selector(epochs(args...))

"""
    participants()
    participants(ids...)
    participants(predicate::Function)

Create a participant ID selection predicate for batch operations.
"""
participants() = x -> fill(true, length(x))

participants(participant_ids::AbstractVector{<:Integer}) =
    let ids = collect(Int, participant_ids)
        x -> Bool[id in ids for id in x]
    end

participants(participant_ids::AbstractRange{<:Integer}) = participants(collect(Int, participant_ids))

participants(participant_id::Integer) = participants([Int(participant_id)])

participants(participant_ranges::AbstractVector{<:AbstractRange{<:Integer}}) = participants(collect(Int, union(participant_ranges...)))

participants(participant_ids::Integer...) = participants(collect(Int, participant_ids))

participants(predicate::Function) = predicate

"""
    participants_not(args...)

Create a participant ID selection predicate excluding the specified IDs (inverse of `participants`).
"""
participants_not(args...) = _negate_selector(participants(args...))

"""
    conditions()
    conditions(indices_or_names...)
    conditions(predicate::Function)

Create a condition selection predicate by condition number, condition name, or custom function.
"""
conditions() = x -> fill(true, length(x))

conditions(condition_indices::AbstractVector{<:Integer}) =
    let idxs = collect(Int, condition_indices)
        x -> _match_condition_indices(x, idxs)
    end

conditions(condition_indices::AbstractRange{<:Integer}) = conditions(collect(Int, condition_indices))

conditions(condition_index::Integer) = conditions([Int(condition_index)])

conditions(condition_ranges::AbstractVector{<:AbstractRange{<:Integer}}) = conditions(collect(Int, union(condition_ranges...)))

conditions(condition_indices::Integer...) = conditions(collect(Int, condition_indices))

conditions(cond_names::AbstractString...) = conditions(collect(String, cond_names))

conditions(rx::Regex) = x -> Bool[occursin(rx, condition_name(dat)) for dat in x]

conditions(predicate::Function) = predicate

conditions(condition_names::AbstractVector{<:AbstractString}) =
    let nms = collect(String, condition_names)
        x -> _match_condition_names(x, nms)
    end

conditions(cond_name::AbstractString) = conditions([String(cond_name)])

"""
    conditions_not(args...)

Create a condition selection predicate excluding the specified conditions (inverse of `conditions`).
"""
conditions_not(args...) = _negate_selector(conditions(args...))

"""
    get_selected_channels(dat, channel_selection; include_meta = true, include_extra = true) -> Vector{Symbol}

Resolve `channel_selection` on `dat` into ordered column symbols, optionally prepending metadata columns and including extra columns.
"""
get_selected_channels(dat, ::Colon; include_meta::Bool = true, include_extra::Bool = true) =
    get_selected_channels(dat, channels(); include_meta = include_meta, include_extra = include_extra)

get_selected_channels(dat, channel_selection; include_meta::Bool = true, include_extra::Bool = true) =
    get_selected_channels(dat, _to_channel_func(channel_selection); include_meta = include_meta, include_extra = include_extra)

function get_selected_channels(dat, channel_selection::Function; include_meta::Bool = true, include_extra::Bool = true)
    channel_selection === Colon() && (channel_selection = channels())
    metadata_cols = include_meta ? meta_labels(dat) : Symbol[]
    selectable_cols = include_extra ? vcat(channel_labels(dat), extra_labels(dat)) : channel_labels(dat)

    selected = selectable_cols[channel_selection(selectable_cols)]

    if channel_selection isa _ChannelNamesSelector
        ordered = _handle_channel_names_order(channel_selection.channel_names, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    elseif channel_selection isa _ChannelNumbersSelector
        ordered = _handle_channel_numbers_order(channel_selection.channel_numbers, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    end

    return isempty(metadata_cols) ? selected : vcat(metadata_cols, selected)
end

"""
    get_selected_components(ica_result::InfoIca, component_selection) -> Vector{Int}

Resolve `component_selection` into 1-based ICA component indices.
"""
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    component_selection === Colon() && return collect(1:length(ica_result.ica_label))
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end

get_selected_components(ica_result::InfoIca, ::Union{Nothing,Colon}) = collect(1:length(ica_result.ica_label))

"""
    get_selected_samples(dat, sample_selection) -> Vector{Int}

Resolve `sample_selection` into 1-based row indices.
"""
get_selected_samples(dat::SingleDataFrameEeg, ::Union{Nothing,AllSelection,Colon,_AllSamples}) = collect(1:nrow(dat.data))

get_selected_samples(dat::SingleDataFrameEeg, s::_SampleIndicesSelector) = _resolve_sample_indices(nrow(dat.data), s.indices)

function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    sample_selection === Colon() && return collect(1:nrow(dat.data))
    res = sample_selection(dat.data)
    return res isa AbstractVector{Bool} ? findall(res) : res
end

get_selected_samples(dat::MultiDataFrameEeg, ::Union{Nothing,AllSelection,Colon,_AllSamples}) =
    isempty(dat.data) ? Int[] : collect(1:nrow(dat.data[1]))

get_selected_samples(dat::MultiDataFrameEeg, s::_SampleIndicesSelector) =
    isempty(dat.data) ? Int[] : _resolve_sample_indices(nrow(dat.data[1]), s.indices)

function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    isempty(dat.data) && return Int[]
    sample_selection === Colon() && return collect(1:nrow(dat.data[1]))
    res = sample_selection(dat.data[1])
    return res isa AbstractVector{Bool} ? findall(res) : res
end

get_selected_samples(dat::DataFrame, ::Union{Nothing,AllSelection,Colon,_AllSamples}) = collect(1:nrow(dat))

get_selected_samples(dat::DataFrame, s::_SampleIndicesSelector) = _resolve_sample_indices(nrow(dat), s.indices)

function get_selected_samples(dat::DataFrame, sample_selection::Function)
    sample_selection === Colon() && return collect(1:nrow(dat))
    res = sample_selection(dat)
    return res isa AbstractVector{Bool} ? findall(res) : res
end

get_selected_samples(dat::Union{SingleDataFrameEeg,MultiDataFrameEeg,DataFrame}, sel::Union{TimeSelection,Tuple{Real,Real}}) =
    get_selected_samples(dat, _interval_to_samples(sel))

get_selected_samples(dat::Union{SingleDataFrameEeg,MultiDataFrameEeg,DataFrame}, sample_selection) =
    get_selected_samples(dat, _to_sample_func(sample_selection))

"""
    get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection) -> Vector{Int}

Resolve `epoch_selection` into 1-based epoch indices.
"""
get_selected_epochs(dat::MultiDataFrameEeg, ::Union{Nothing,Colon}) = collect(1:length(dat.data))

get_selected_epochs(dat::MultiDataFrameEeg, s::_EpochNumbersSelector) = _resolve_sample_indices(length(dat.data), s.epoch_numbers)

get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection) = get_selected_epochs(dat, _to_epoch_func(epoch_selection))

function get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection::Function)
    epoch_selection === Colon() && return collect(1:length(dat.data))
    res = epoch_selection(dat.data)
    return res isa AbstractVector{Bool} ? findall(res) : res
end

"""
    get_selected_conditions(datasets::AbstractVector{<:EegFunData}, condition_selection) -> Vector{Int}

Resolve `condition_selection` into 1-based condition dataset indices.
"""
get_selected_conditions(datasets::AbstractVector{<:EegFunData}, ::Union{Nothing,Colon}) = collect(1:length(datasets))

get_selected_conditions(datasets::AbstractVector{<:EegFunData}, condition_selection) =
    get_selected_conditions(datasets, _to_condition_func(condition_selection))

function get_selected_conditions(datasets::AbstractVector{<:EegFunData}, condition_selection::Function)
    condition_selection === Colon() && return collect(1:length(datasets))
    res = condition_selection(datasets)
    return res isa AbstractVector{Bool} ? findall(res) : res
end

"""
    epoch_to_continuous(dat::MultiDataFrameEeg, epoch_idx::Int) -> ContinuousData

Extract epoch `epoch_idx` from `dat` as an independent `ContinuousData` object.
"""
function epoch_to_continuous(dat::MultiDataFrameEeg, epoch_idx::Int)
    1 <= epoch_idx <= length(dat.data) || @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    return ContinuousData(dat.file, copy(dat.data[epoch_idx]; copycols = true), copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))
end

######################################
######### Internal Functions #########
######################################

_subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame =
    df[selected_samples, selected_channels]

_subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame} = DataFrame[_subset_dataframe(dataframes[ep], selected_channels, selected_samples) for ep in selected_epochs]

(::_AllSamples)(x) = fill(true, nrow(x))

(s::_ChannelNamesSelector)(x) = x .∈ Ref(s.channel_names)

(s::_ChannelNumbersSelector)(x) = Bool[i in s.channel_numbers for i = 1:length(x)]

function (s::_EpochNumbersSelector)(x)
    n = length(x)
    mask = fill(false, n)
    @inbounds for idx in s.epoch_numbers
        1 <= idx <= n && (mask[idx] = true)
    end
    return mask
end

function _resolve_sample_indices(n_rows::Int, indices::Vector{Int})::Vector{Int}
    @inbounds for idx in indices
        (1 <= idx <= n_rows) || return Int[i for i in indices if 1 <= i <= n_rows]
    end
    return copy(indices)
end

function (s::_SampleIndicesSelector)(x::DataFrame)
    n = nrow(x)
    mask = falses(n)
    @inbounds for idx in s.indices
        1 <= idx <= n && (mask[idx] = true)
    end
    return mask
end

function _subset_common(dat, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting channels and samples"
    meta_cols = meta_labels(dat)
    chosen = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    eeg_set = Set{Symbol}(channel_labels(dat))
    eeg_selected = Symbol[ch for ch in chosen if ch in eeg_set]
    extra_selected = include_extra ? Symbol[ch for ch in chosen if ch ∉ eeg_set] : Symbol[]
    selected_channels = vcat(meta_cols, eeg_selected, extra_selected)
    selected_samples = get_selected_samples(dat, sample_selection)
    layout_subset = subset_layout(dat.layout; channel_selection = _ChannelNamesSelector(eeg_selected))
    return selected_channels, selected_samples, layout_subset
end

function _subset_common(dat::MultiDataFrameEeg, epoch_selection, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting epochs, channels and samples"
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    return selected_epochs, selected_channels, selected_samples, layout_subset
end

function _time_interval_mask(t::AbstractVector{<:Real}, start_time::Float64, stop_time::Float64)
    n_rows = length(t)
    mask = BitVector(undef, n_rows)
    n_sel = 0
    first_t = 0.0
    last_t = 0.0
    @inbounds for i = 1:n_rows
        ti = Float64(t[i])
        hit = ti >= start_time && ti <= stop_time
        mask[i] = hit
        if hit
            n_sel == 0 && (first_t = ti)
            last_t = ti
            n_sel += 1
        end
    end

    if n_sel > 0
        @debug "times(): Requested $(start_time)s–$(stop_time)s, actual sample range $(round(first_t; digits = 4))s–$(round(last_t; digits = 4))s ($(n_sel) sample$(n_sel > 1 ? "s" : ""))"
    elseif start_time == stop_time && n_rows > 0
        best_idx = 1
        best_dist = abs(Float64(t[1]) - start_time)
        @inbounds for i = 2:n_rows
            d = abs(Float64(t[i]) - start_time)
            if d < best_dist
                best_dist = d
                best_idx = i
            end
        end
        snapped_t = t[best_idx]
        @debug "times(): Requested $(start_time)s, snapped to nearest sample at $(round(snapped_t; digits = 4))s"
        @inbounds for i = 1:n_rows
            mask[i] = (t[i] == snapped_t)
        end
    else
        @minimal_warning "times(): No samples found in range $(start_time)s–$(stop_time)s"
    end
    return mask
end

function _interval_to_samples(sel::TimeSelection)
    start_time, stop_time = sel.start, sel.stop
    start_time <= stop_time || @minimal_error("Time interval start ($start_time) must be <= stop ($stop_time)")
    return x -> _time_interval_mask(x[!, :time], start_time, stop_time)
end

_interval_to_samples(::AllSelection) = samples()

_interval_to_samples(::Nothing) = samples()

_interval_to_samples(t::Tuple{Real,Real}) = _interval_to_samples(times(t[1], t[2]))

function _interval_to_samples(v::AbstractVector{<:Real})
    length(v) == 2 || @minimal_error("Time interval vector must have exactly 2 elements [start, stop], got $(length(v))")
    return _interval_to_samples(times(v[1], v[2]))
end

_interval_to_samples(f::Function) = f

function _combine_interval_sample(interval::Interval, sample_selection)
    sample_func = _to_sample_func(sample_selection)
    (interval isa Nothing || interval isa AllSelection) && return sample_func
    interval_sel = _interval_to_samples(interval)
    sample_func isa _AllSamples && return interval_sel
    if sample_func isa _SampleIndicesSelector
        idxs = sample_func.indices
        return x -> begin
            mask = interval_sel(x)
            n = length(mask)
            Int[i for i in idxs if 1 <= i <= n && mask[i]]
        end
    end
    return x -> interval_sel(x) .& sample_func(x)
end

_create_subset(dat::ContinuousData, ds, ls) = ContinuousData(dat.file, ds, ls, dat.sample_rate, copy(dat.analysis_info))

_create_subset(dat::ErpData, ds, ls) =
    ErpData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info), dat.n_epochs)

_create_subset(dat::EpochData, ds, ls) =
    EpochData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info))

_create_subset(dat::SpectrumData, ds, ls) =
    SpectrumData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, dat.method, copy(dat.analysis_info))

_create_subset(dat::TimeFreqData, ds_power::DataFrame, ds_phase::DataFrame, ls) = TimeFreqData(
    dat.file,
    dat.condition,
    dat.condition_name,
    ds_power,
    ds_phase,
    ls,
    dat.sample_rate,
    dat.method,
    dat.baseline,
    copy(dat.analysis_info),
)

_create_subset(dat::TimeFreqEpochData, ds_power::Vector{DataFrame}, ds_phase::Vector{DataFrame}, ls) = TimeFreqEpochData(
    dat.file,
    dat.condition,
    dat.condition_name,
    ds_power,
    ds_phase,
    ls,
    dat.sample_rate,
    dat.method,
    dat.baseline,
    copy(dat.analysis_info),
)

_negate_selector(f::Function) = x -> begin
    res = f(x)
    if res isa AbstractVector{Bool}
        return .!res
    elseif res isa AbstractVector{<:Integer}
        mask = fill(true, length(x))
        for idx in res
            1 <= idx <= length(x) && (mask[idx] = false)
        end
        return mask
    else
        return .!(res)
    end
end

@inline function _copy_bools!(out::Vector{Bool}, col::AbstractVector)
    @inbounds for i in eachindex(out, col)
        out[i] = Bool(col[i])
    end
end

@inline function _or_bools!(out::Vector{Bool}, col::AbstractVector)
    @inbounds for i in eachindex(out, col)
        out[i] |= Bool(col[i])
    end
end

@inline function _and_bools!(out::Vector{Bool}, col::AbstractVector)
    @inbounds for i in eachindex(out, col)
        out[i] &= Bool(col[i])
    end
end

function _combine_sample_cols(df::DataFrame, columns::Vector{Symbol}, is_or::Bool, negate::Bool)
    n = nrow(df)
    isempty(columns) && return fill(negate || !is_or, n)
    out = Vector{Bool}(undef, n)
    _copy_bools!(out, df[!, columns[1]])
    if is_or
        for k = 2:length(columns)
            _or_bools!(out, df[!, columns[k]])
        end
    else
        for k = 2:length(columns)
            _and_bools!(out, df[!, columns[k]])
        end
    end
    if negate
        @inbounds for i = 1:n
            out[i] = !out[i]
        end
    end
    return out
end

function _match_condition_indices(x, condition_indices::AbstractVector{Int})::Vector{Int}
    isempty(x) && return Int[]
    if all(d -> hasfield(typeof(d), :condition) && d.condition isa Integer && d.condition > 0, x)
        cond_nums = Int[Int(d.condition) for d in x]
        req = unique(condition_indices)
        present_req = Int[c for c in req if c in cond_nums]
        isempty(present_req) && return Int[]
        first_pos = Int[findfirst(==(c), cond_nums)::Int for c in present_req]
        if issorted(first_pos)
            keep = Set{Int}(present_req)
            return Int[i for (i, cn) in enumerate(cond_nums) if cn in keep]
        end
        return Int[i for c in present_req for (i, cn) in enumerate(cond_nums) if cn == c]
    end
    return intersect(condition_indices, 1:length(x))
end

function _match_condition_names(x, condition_names::Vector{String})::Vector{Int}
    isempty(x) && return Int[]
    names_in_x = String[condition_name(dat) for dat in x]
    req = unique(condition_names)
    present_req = String[n for n in req if n in names_in_x]
    isempty(present_req) && return Int[]
    first_pos = Int[findfirst(==(n), names_in_x)::Int for n in present_req]
    if issorted(first_pos)
        keep = Set{String}(present_req)
        return Int[i for (i, nm) in enumerate(names_in_x) if nm in keep]
    end
    return Int[i for n in present_req for (i, nm) in enumerate(names_in_x) if nm == n]
end

_to_sample_func(::Nothing) = samples()

_to_sample_func(::Colon) = samples()

_to_sample_func(f::Function) = f === Colon() ? samples() : f

_to_sample_func(col::Symbol) = samples(col)

_to_sample_func(mask::AbstractVector{Bool}) =
    x ->
        length(mask) == nrow(x) ? collect(Bool, mask) :
        @minimal_error("Boolean sample mask length ($(length(mask))) must match number of samples ($(nrow(x)))")

_to_sample_func(idx::Integer) = _SampleIndicesSelector(Int[Int(idx)])

_to_sample_func(idx::AbstractVector{<:Integer}) = _SampleIndicesSelector(collect(Int, idx))

_to_sample_func(idx::AbstractRange{<:Integer}) = _SampleIndicesSelector(collect(Int, idx))

function _to_sample_func(idx)
    set = Set(idx)
    return x -> Bool[i in set for i = 1:nrow(x)]
end

_to_channel_func(::Nothing) = channels()

_to_channel_func(::Colon) = channels()

_to_channel_func(f::Function) = f === Colon() ? channels() : f

_to_channel_func(ch::Symbol) = channels(ch)

_to_channel_func(chans::AbstractVector{Symbol}) = channels(collect(chans))

_to_channel_func(idx::AbstractString) = channels(Symbol(idx))

_to_channel_func(idx::AbstractVector{<:AbstractString}) = channels(Symbol.(idx))

_to_channel_func(mask::AbstractVector{Bool}) = channels(findall(mask))

_to_channel_func(idx::Integer) = channels(Int(idx))

_to_channel_func(idx::AbstractVector{<:Integer}) = channels(collect(Int, idx))

_to_channel_func(idx::AbstractRange) = channels(collect(Int, idx))

_to_channel_func(rx::Regex) = channels(rx)

_to_channel_func(idx) = channels(idx)

_to_epoch_func(::Nothing) = epochs()

_to_epoch_func(::Colon) = epochs()

_to_epoch_func(f::Function) = f === Colon() ? epochs() : f

_to_epoch_func(mask::AbstractVector{Bool}) = epochs(findall(mask))

_to_epoch_func(idx::Integer) = epochs(Int(idx))

_to_epoch_func(idx::AbstractVector{<:Integer}) = epochs(collect(Int, idx))

_to_epoch_func(idx::AbstractRange) = epochs(collect(Int, idx))

_to_epoch_func(idx) = epochs(idx)

_to_condition_func(::Nothing) = conditions()

_to_condition_func(::Colon) = conditions()

_to_condition_func(f::Function) = f === Colon() ? conditions() : f

_to_condition_func(mask::AbstractVector{Bool}) =
    x ->
        length(mask) == length(x) ? collect(Bool, mask) :
        @minimal_error("Boolean condition mask length ($(length(mask))) must match number of datasets ($(length(x)))")

_to_condition_func(idx::Integer) = conditions(Int(idx))

_to_condition_func(idx::AbstractVector{<:Integer}) = conditions(collect(Int, idx))

_to_condition_func(idx::AbstractRange) = conditions(collect(Int, idx))

_to_condition_func(name::AbstractString) = conditions(String(name))

_to_condition_func(names::AbstractVector{<:AbstractString}) = conditions(collect(String, names))

_to_condition_func(rx::Regex) = conditions(rx)

_to_condition_func(idx) = conditions(idx)

_to_participant_func(::Nothing) = participants()

_to_participant_func(::Colon) = participants()

_to_participant_func(f::Function) = f === Colon() ? participants() : f

_to_participant_func(mask::AbstractVector{Bool}) =
    x ->
        length(mask) == length(x) ? collect(Bool, mask) :
        @minimal_error("Boolean participant mask length ($(length(mask))) must match number of files ($(length(x)))")

_to_participant_func(idx::Integer) = participants(Int(idx))

_to_participant_func(idx::AbstractVector{<:Integer}) = participants(collect(Int, idx))

_to_participant_func(idx::AbstractRange) = participants(collect(Int, idx))

_to_participant_func(idx) = participants(idx)

function _handle_channel_names_order(user_order::AbstractVector{Symbol}, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    ordered = sizehint!(Symbol[], length(user_order))
    all_selected = true
    for ch in user_order
        if ch ∉ selectable_cols
            @minimal_warning "Channel $(ch) not found in data!"
        elseif ch ∈ ordered
            @minimal_warning "Channel $(ch) already specified!"
        else
            push!(ordered, ch)
            all_selected &= (ch ∈ selected)
        end
    end
    return (!isempty(ordered) && all_selected) ? ordered : nothing
end

function _handle_channel_numbers_order(user_order_numbers, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    n_cols = length(selectable_cols)
    seen = sizehint!(Int[], length(user_order_numbers))
    ordered = sizehint!(Symbol[], length(user_order_numbers))
    all_selected = true
    for i in user_order_numbers
        if i < 1 || i > n_cols
            @minimal_warning "Channel index $(i) out of range (valid: 1:$n_cols)!"
        elseif i ∈ seen
            @minimal_warning "Channel index $(i) already specified!"
        else
            push!(seen, i)
            ch = selectable_cols[i]
            ch ∉ ordered && push!(ordered, ch)
            all_selected &= (ch ∈ selected)
        end
    end
    return (!isempty(ordered) && all_selected) ? ordered : nothing
end
