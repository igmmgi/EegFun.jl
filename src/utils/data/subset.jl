_subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame =
    df[selected_samples, selected_channels]

_subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame} = DataFrame[_subset_dataframe(dataframes[ep], selected_channels, selected_samples) for ep in selected_epochs]

struct _AllSamples <: Function end
(::_AllSamples)(x) = fill(true, nrow(x))

function _subset_common(dat, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting channels and samples"
    selected_channels = get_selected_channels(dat, channel_selection; include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    layout_subset = subset_layout(dat.layout; channel_selection = channel_selection)
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
        @debug "times(): Requested $(start_time)s, snapped to nearest sample at $(round(t[best_idx]; digits = 4))s"
        mask[best_idx] = true
    else
        @minimal_warning "times(): No samples found in range $(start_time)s–$(stop_time)s"
    end
    return mask
end

function _interval_to_samples(sel::TimeSelection)
    start_time, stop_time = sel.start, sel.stop
    return x -> _time_interval_mask(x[!, :time], start_time, stop_time)
end

_interval_to_samples(::AllSelection) = samples()
_interval_to_samples(::Nothing) = samples()
_interval_to_samples(t::Tuple{Real,Real}) = _interval_to_samples(TimeSelection(t[1], t[2]))
function _interval_to_samples(v::AbstractVector{<:Real})
    length(v) == 2 || @minimal_error("Time interval vector must have exactly 2 elements [start, stop], got $(length(v))")
    return _interval_to_samples(TimeSelection(v[1], v[2]))
end
_interval_to_samples(f::Function) = f

function _combine_interval_sample(interval::Interval, sample_selection)
    sample_func = _to_sample_func(sample_selection)
    (interval isa Nothing || interval isa AllSelection) && return sample_func
    interval_sel = _interval_to_samples(interval)
    sample_func isa _AllSamples && return interval_sel
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

# === SELECTION PREDICATE GENERATORS ===

"""
    channels()
    channels(names::Symbol...)
    channels(names::AbstractVector{Symbol})
    channels(indices)
    channels(predicate::Function)

Create a channel selection predicate for `channel_selection` arguments.
"""
channels() = x -> fill(true, length(x))
channels(channel_names::AbstractVector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_names::AbstractVector{<:AbstractString}) = channels(Symbol.(channel_names))
channels(channel_name::Symbol) = channels([channel_name])
channels(channel_name::AbstractString) = channels([Symbol(channel_name)])
channels(channel_names::Symbol...) = channels(collect(channel_names))
channels(channel_names::AbstractString...) = channels(Symbol.(channel_names)...)
channels(channel_number::Int) = channels([channel_number])
channels(channel_numbers::AbstractVector{Int}) = x -> Bool[i in channel_numbers for i = 1:length(x)]
channels(channel_numbers::UnitRange{Int}) = channels(collect(Int, channel_numbers))
channels(channel_ranges::Vector{UnitRange{Int}}) = channels(collect(Int, union(channel_ranges...)))
channels(predicate::Function) = predicate

# Note: keep captured variable named `f` (inspected by grand_average condition checking).
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
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> Bool[i in component_numbers for i = 1:length(x)]
components(component_number::Int) = components([component_number])
components(component_ranges::Vector{UnitRange{Int}}) = components(union(component_ranges...))
components(component_numbers::Int...) = components(collect(component_numbers))
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
times(start::Real, stop::Real) = TimeSelection(Float64(start), Float64(stop))
times(interval::Tuple{Real,Real}) = TimeSelection(Float64(interval[1]), Float64(interval[2]))

"""
    epochs()
    epochs(indices...)
    epochs(predicate::Function)

Create an epoch selection predicate for `epoch_selection` arguments.
"""
epochs() = x -> fill(true, length(x))
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> Bool[i in epoch_numbers for i = 1:length(x)]
epochs(epoch_number::Int) = epochs([epoch_number])
epochs(epoch_ranges::Vector{UnitRange{Int}}) = epochs(union(epoch_ranges...))
epochs(epoch_numbers::Int...) = epochs(collect(epoch_numbers))
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
participants(participant_ids::AbstractVector{Int}) = x -> Bool[id in participant_ids for id in x]
participants(participant_ids::UnitRange{Int}) = participants(collect(Int, participant_ids))
participants(participant_id::Int) = participants([participant_id])
participants(participant_ranges::Vector{UnitRange{Int}}) = participants(collect(Int, union(participant_ranges...)))
participants(participant_ids::Int...) = participants(collect(Int, participant_ids))
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

Create a condition selection predicate by index, name, or custom function.
"""
conditions() = x -> fill(true, length(x))
# Note: keep captured variable named `condition_indices` (inspected by grand_average).
conditions(condition_indices::AbstractVector{Int}) = x -> intersect(condition_indices, 1:length(x))
conditions(condition_indices::UnitRange{Int}) = conditions(collect(Int, condition_indices))
conditions(condition_index::Int) = conditions([condition_index])
conditions(condition_ranges::Vector{UnitRange{Int}}) = conditions(collect(Int, union(condition_ranges...)))
conditions(condition_indices::Int...) = conditions(collect(Int, condition_indices))
conditions(condition_names::Vector{String}) =
    x -> Int[idx for n in condition_names for idx in (findfirst(dat -> condition_name(dat) == n, x),) if !isnothing(idx)]
conditions(cond_name::String) = conditions([cond_name])
conditions(cond_names::String...) = conditions(collect(cond_names))
conditions(predicate::Function) = predicate

"""
    conditions_not(args...)

Create a condition selection predicate excluding the specified conditions (inverse of `conditions`).
"""
conditions_not(args...) = _negate_selector(conditions(args...))

# === SELECTOR CONVERSION & APPLICATION ===

_to_sample_func(::Nothing) = samples()
_to_sample_func(::Colon) = samples()
_to_sample_func(f::Function) = f === Colon() ? samples() : f
_to_sample_func(col::Symbol) = samples(col)
_to_sample_func(mask::AbstractVector{Bool}) =
    x ->
        length(mask) == nrow(x) ? collect(Bool, mask) :
        @minimal_error("Boolean sample mask length ($(length(mask))) must match number of samples ($(nrow(x)))")
_to_sample_func(idx::Integer) = x -> Bool[i == idx for i = 1:nrow(x)]
_to_sample_func(idx) = x -> Bool[i in idx for i = 1:nrow(x)]

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
_to_condition_func(mask::AbstractVector{Bool}) = conditions(findall(mask))
_to_condition_func(idx::Integer) = conditions(Int(idx))
_to_condition_func(idx::AbstractVector{<:Integer}) = conditions(collect(Int, idx))
_to_condition_func(idx::AbstractRange) = conditions(collect(Int, idx))
_to_condition_func(name::AbstractString) = conditions(String(name))
_to_condition_func(names::AbstractVector{<:AbstractString}) = conditions(collect(String, names))
_to_condition_func(idx) = conditions(idx)

_to_participant_func(::Nothing) = participants()
_to_participant_func(::Colon) = participants()
_to_participant_func(f::Function) = f === Colon() ? participants() : f
_to_participant_func(mask::AbstractVector{Bool}) = x -> collect(Bool, mask)
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

    selection_type = typeof(channel_selection)
    if hasfield(selection_type, :channel_names)
        ordered = _handle_channel_names_order(getfield(channel_selection, :channel_names), selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    elseif hasfield(selection_type, :channel_numbers)
        ordered = _handle_channel_numbers_order(getfield(channel_selection, :channel_numbers), selectable_cols, selected)
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
function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    sample_selection === Colon() && return collect(1:nrow(dat.data))
    return findall(sample_selection(dat.data))
end

get_selected_samples(dat::MultiDataFrameEeg, ::Union{Nothing,AllSelection,Colon,_AllSamples}) =
    isempty(dat.data) ? Int[] : collect(1:nrow(dat.data[1]))
function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    isempty(dat.data) && return Int[]
    sample_selection === Colon() && return collect(1:nrow(dat.data[1]))
    return findall(sample_selection(dat.data[1]))
end

get_selected_samples(dat::DataFrame, ::Union{Nothing,AllSelection,Colon,_AllSamples}) = collect(1:nrow(dat))
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    sample_selection === Colon() && return collect(1:nrow(dat))
    return findall(sample_selection(dat))
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
