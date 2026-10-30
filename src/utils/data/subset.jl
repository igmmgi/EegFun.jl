# === DATAFRAME SUBSETTING UTILITIES ===

"""
    _subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> DataFrame

Create a subset of a DataFrame by selecting specific channels and samples.

# Arguments
- `df::DataFrame`: The DataFrame to subset
- `selected_channels::Vector{Symbol}`: Column names to include
- `selected_samples::Vector{Int}`: Row indices to include

# Returns
- `DataFrame`: Subset DataFrame with selected channels and samples
"""
function _subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame
    return df[selected_samples, selected_channels]
end

"""
    _subset_dataframes(dataframes::Vector{DataFrame}, selected_epochs::Vector{Int}, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> Vector{DataFrame}

Create subsets of multiple DataFrames by selecting specific epochs, channels, and samples.

# Arguments
- `dataframes::Vector{DataFrame}`: Vector of DataFrames to subset
- `selected_epochs::Vector{Int}`: Indices of epochs (DataFrames) to include
- `selected_channels::Vector{Symbol}`: Column names to include in each DataFrame
- `selected_samples::Vector{Int}`: Row indices to include in each DataFrame

# Returns
- `Vector{DataFrame}`: Vector of subset DataFrames
"""
function _subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame}
    return [_subset_dataframe(dataframes[ep], selected_channels, selected_samples) for ep in selected_epochs]
end



# === INTERNAL SUBSET HELPERS ===

"""
    _subset_common(dat, channel_selection, sample_selection, include_extra)

Internal helper for common subsetting operations on EEG data.

# Arguments
- `dat`: EEG data object
- `channel_selection::Function`: Channel selection predicate
- `sample_selection::Function`: Sample selection predicate
- `include_extra::Bool`: Whether to include extra columns

# Returns
- `Tuple`: (selected_channels, selected_samples, layout_subset)
"""
function _subset_common(dat, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting channels and samples"
    # Get subset selected channels, samples, and layout
    selected_channels = get_selected_channels(dat, channel_selection, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    layout_subset = subset_layout(dat.layout, channel_selection = channel_selection)

    return selected_channels, selected_samples, layout_subset
end

"""
    _subset_common(dat::EpochData, epoch_selection, channel_selection, sample_selection, include_extra)

Internal helper for subsetting EpochData with epoch selection.

# Arguments
- `dat::EpochData`: Epoch data object
- `epoch_selection::Function`: Epoch selection predicate
- `channel_selection::Function`: Channel selection predicate
- `sample_selection::Function`: Sample selection predicate
- `include_extra::Bool`: Whether to include extra columns

# Returns
- `Tuple`: (selected_epochs, selected_channels, selected_samples, layout_subset)
"""
function _subset_common(dat::MultiDataFrameEeg, epoch_selection, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting epochs, channels and samples"
    # Get subset selected epochs, channels, samples, and layout
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    return selected_epochs, selected_channels, selected_samples, layout_subset
end

"""
    _interval_to_samples(interval::Interval)

Internal helper to convert a Interval (nothing, tuple, or range)
into a samples() predicate function for use in subset().

# Arguments
- `interval::Interval`: Time interval specification

# Returns
- `Function`: samples() predicate that selects the specified time range
"""
function _interval_to_samples(sel::TimeSelection)
    start_time, stop_time = sel.start, sel.stop
    return x -> begin
        t = x[!, :time]
        mask = (t .>= start_time) .& (t .<= stop_time)
        selected_times = t[mask]
        if !isempty(selected_times)
            actual_start = round(first(selected_times); digits = 4)
            actual_stop = round(last(selected_times); digits = 4)
            n = length(selected_times)
            @debug "times(): Requested $(start_time)s–$(stop_time)s, actual sample range $(actual_start)s–$(actual_stop)s ($(n) sample$(n > 1 ? "s" : ""))"
        elseif start_time == stop_time
            # Single time point with no exact match: snap to nearest sample
            _, idx = findmin(abs.(t .- start_time))
            @debug "times(): Requested $(start_time)s, snapped to nearest sample at $(round(t[idx]; digits=4))s"
            mask = fill(false, nrow(x))
            mask[idx] = true
        else
            @minimal_warning "times(): No samples found in range $(start_time)s–$(stop_time)s"
        end
        mask
    end
end

_interval_to_samples(::AllSelection) = samples()
_interval_to_samples(::Nothing) = samples()
_interval_to_samples(t::Tuple{Real,Real}) = _interval_to_samples(TimeSelection(t[1], t[2]))
function _interval_to_samples(v::AbstractVector{<:Real})
    if length(v) != 2
        @minimal_error("Time interval vector must have exactly 2 elements [start, stop], got $(length(v))")
    end
    return _interval_to_samples(TimeSelection(v[1], v[2]))
end
_interval_to_samples(f::Function) = f

"""
    _combine_interval_sample(interval::Interval, sample_selection::Function)

Combine an interval selection with a sample selection predicate.
"""
function _combine_interval_sample(interval::Interval, sample_selection)
    sample_func = _to_sample_func(sample_selection)
    if interval isa Nothing || interval isa AllSelection
        return sample_func
    else
        interval_sel = _interval_to_samples(interval)
        return x -> interval_sel(x) .& sample_func(x)
    end
end

"""Construct a new data object of the same type with subsetted data and layout."""
_create_subset(dat::ContinuousData, ds, ls) = ContinuousData(dat.file, ds, ls, dat.sample_rate, copy(dat.analysis_info))
_create_subset(dat::ErpData, ds, ls) =
    ErpData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info), dat.n_epochs)
_create_subset(dat::EpochData, ds, ls) =
    EpochData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info))
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

# === SUBSET IMPLEMENTATIONS ===

"""
    subset(dat::SingleDataFrameEeg; channel_selection, sample_selection, interval_selection, include_extra)

Create a subset of single-DataFrame EEG data (ContinuousData, ErpData).

# Arguments
- `channel_selection::Function`: Channel predicate (default: `channels()` for all)
- `sample_selection::Function`: Sample predicate (default: `samples()` for all)
- `interval_selection::Interval`: Time interval (default: `times()` for all)
- `include_extra::Bool`: Include extra columns (default: `false`)

# Examples
```julia
# Select specific channels and time window
sub = subset(erp, channel_selection = channels(:Fz, :Cz), interval_selection = times(-0.2, 0.5))
```
"""
function subset(
    dat::T;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::T where {T<:SingleDataFrameEeg}
    # Combine interval and sample selection
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
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, combined_sel, include_extra)
    # Subset BOTH power and phase DataFrames
    power_subset = _subset_dataframe(dat.data_power, selected_channels, selected_samples)
    phase_subset = _subset_dataframe(dat.data_phase, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

"""
    subset(dat::MultiDataFrameEeg; channel_selection, sample_selection, interval_selection, epoch_selection, include_extra)

Create a subset of epoched EEG data (e.g., `EpochData`). 

# Arguments
- `epoch_selection`: Epoch predicate or indices (default: `epochs()` for all).
- `channel_selection`: Channel predicate, Symbol, or indices (default: `channels()` for all)
- `sample_selection`: Sample predicate (default: `samples()` for all)
- `interval_selection::Interval`: Time interval (default: `times()` for all)
- `include_extra::Bool`: Include extra columns (default: `false`)

# Examples
```julia
# Drop trial outliers by examining the interval channel inside the epoch
cleaned_epochs = subset(epochs_data, epoch_selection = epochs(df -> 0.2 < df.interval[1] < 1.0))
```
"""
function subset(
    dat::T;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::T where {T<:MultiDataFrameEeg}
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    dat_subset = _subset_dataframes(dat.data, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, dat_subset, layout_subset)
end

"""
    subset(dat::TimeFreqEpochData; channel_selection, sample_selection, interval_selection, epoch_selection, include_extra)

Create a subset of time-frequency epoched data. Accepts the same arguments as epoched data, including functional `epoch_selection`.
"""
function subset(
    dat::TimeFreqEpochData;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::TimeFreqEpochData
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    # Subset BOTH power and phase Vector{DataFrame}
    power_subset = _subset_dataframes(dat.data_power, selected_epochs, selected_channels, selected_samples)
    phase_subset = _subset_dataframes(dat.data_phase, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

function subset(
    datasets::AbstractVector{<:ErpData};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Vector{ErpData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, and interval selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::AbstractVector{<:EpochData};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::Vector{EpochData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, interval, and epoch selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::AbstractVector{<:TimeFreqData};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Vector{TimeFreqData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, and interval selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::AbstractVector{<:TimeFreqEpochData};
    condition_selection = conditions(),
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    epoch_selection = epochs(),
    include_extra::Bool = false,
)::Vector{TimeFreqEpochData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, interval, and epoch selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
end


# === SELECTION HELPERS ===

"""
    channels()
    channels(name::Symbol)
    channels(names::Symbol...)
    channels(names::Vector{Symbol})
    channels(index::Int)
    channels(indices::Union{Vector{Int}, UnitRange})
    channels(ranges::Vector{UnitRange{Int}})
    channels(predicate::Function)

Predicate generators for channel selection in subsetting and plotting operations.

Returns a function that, given a vector of channel labels, produces a boolean mask
selecting the matching channels. Used with the `channel_selection` keyword argument.

# Examples
```julia
channels()                          # Select all channels (default)
channels(:Cz)                       # Select a single channel by name
channels(:Fz, :Cz, :Pz)            # Select multiple channels by name
channels([:Fz, :Cz, :Pz])          # Same, from a Vector
channels(1:10)                      # Select first 10 channels by index
channels(x -> startswith.(string.(x), "F"))  # Custom predicate
```
"""
channels() = x -> fill(true, length(x))  # Default: select all channels given
channels(channel_names::AbstractVector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_names::AbstractVector{<:AbstractString}) = channels(Symbol.(channel_names))
channels(channel_name::Symbol) = channels([channel_name])
channels(channel_name::AbstractString) = channels([Symbol(channel_name)])
channels(channel_names::Symbol...) = channels(collect(channel_names))
channels(channel_names::AbstractString...) = channels(Symbol.(channel_names)...)
channels(channel_number::Int) = channels([channel_number])
channels(channel_numbers::AbstractVector{Int}) = x -> [i in channel_numbers for i = 1:length(x)]
channels(channel_numbers::UnitRange{Int}) = channels(collect(Int, channel_numbers))
channels(channel_ranges::Vector{UnitRange{Int}}) = channels(collect(Int, union(channel_ranges...)))
channels(predicate::Function) = predicate  # Allow custom function predicates
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

Predicate generators that **exclude** the specified channels.
Inverse of `channels(...)`.

# Examples
```julia
channels_not(:Fp1, :Fp2)           # Exclude two frontal channels
channels_not(1:4)                   # Exclude first 4 channels by index
```
"""
channels_not(args...) = _negate_selector(channels(args...))

"""
    components()
    components(index::Int)
    components(indices::Int...)
    components(indices::Union{Vector{Int}, UnitRange})
    components(predicate::Function)

Predicate generators for ICA component selection.

# Examples
```julia
components()                        # Select all components (default)
components(1, 2, 3)                 # Select components 1, 2, and 3
components(1:5)                     # Select first 5 components
```
"""
components() = x -> fill(true, length(x))  # Default: select all components given
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> [i in component_numbers for i = 1:length(x)]
components(component_number::Int) = components([component_number])
components(component_ranges::Vector{UnitRange{Int}}) = components(union(component_ranges...))
components(component_numbers::Int...) = components(collect(component_numbers))
components(predicate::Function) = predicate  # Allow custom function predicates

"""
    components_not(index::Int)
    components_not(indices::Int...)
    components_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified ICA components.

# Examples
```julia
components_not(1)                   # Exclude component 1
components_not(1:3)                 # Exclude first 3 components
```
"""
components_not(args...) = _negate_selector(components(args...))

"""
    samples()
    samples(column::Symbol)
    samples(predicate::Function)

Predicate generators for sample (row) selection. Returns a function that,
given a DataFrame, produces a boolean mask over rows.

# Examples
```julia
samples()                           # Select all samples (default)
samples(:is_clean)                  # Select rows where :is_clean is true
samples(df -> df.time .> 0)         # Custom predicate on the DataFrame
```
"""
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
# Allow custom function predicates
samples(predicate::Function) = predicate
"""
    samples_or(columns::Vector{Symbol})
    samples_and(columns::Vector{Symbol})
    samples_not(column::Symbol)
    samples_or_not(columns::Vector{Symbol})
    samples_and_not(columns::Vector{Symbol})

Combine or negate sample predicates across boolean columns.

# Examples
```julia
samples_or([:is_extreme, :is_drift])    # Rows where ANY flag is true
samples_and([:is_clean, :is_valid])     # Rows where ALL flags are true
samples_not(:is_bad)                    # Rows where :is_bad is false
samples_or_not([:is_extreme, :is_bad])  # Rows where NO flag is true
```
"""
samples_or(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(false, nrow(x)) : reduce((a, b) -> a .| b, (x[!, col] for col in columns))
"""See `samples_or`."""
samples_and(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : reduce((a, b) -> a .& b, (x[!, col] for col in columns))

"""See `samples_or`."""
samples_not(column::Symbol) = x -> .!(x[!, column])

"""See `samples_or`."""
samples_or_not(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a, b) -> a .| b, (x[!, col] for col in columns)))

"""See `samples_or`."""
samples_and_not(columns::Vector{Symbol}) =
    x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a, b) -> a .& b, (x[!, col] for col in columns)))

# Helper functions for time interval selection (Interval type)
# These return Interval values (tuples, ranges, or nothing), not predicates
"""
    times()
    times(t::Real)
    times(start::Real, stop::Real)
    times(interval::Tuple{Real,Real})

Create time interval values for the `interval_selection` keyword argument.
Unlike other predicate generators, `times()` returns an `AbstractSelection` value
(a `TimeSelection` or `AllSelection`), not a predicate function.

# Examples
```julia
times()                             # Select all time points (returns AllSelection)
times(0.3)                          # Select single time point (snaps to nearest sample)
times(-0.2, 0.5)                    # Select time window from -200ms to 500ms
times((-0.2, 0.5))                  # Same, from a Tuple
```
"""
times() = AllSelection()
times(t::Real) = TimeSelection(Float64(t), Float64(t))
times(start::Real, stop::Real) = TimeSelection(Float64(start), Float64(stop))
times(interval::Tuple{Real,Real}) = TimeSelection(Float64(interval[1]), Float64(interval[2]))


"""
    epochs()
    epochs(index::Int)
    epochs(indices::Int...)
    epochs(indices::Union{Vector{Int}, UnitRange})
    epochs(predicate::Function)

Predicate generators for epoch selection in subsetting operations.
The `predicate` variant receives individual epoch DataFrames and should return `Bool`.

# Examples
```julia
epochs()                              # Select all epochs (default)
epochs(1:5)                           # Select first 5 epochs by numerical index
epochs(2, 4, 6)                       # Select epochs 2, 4, and 6
epochs(df -> nrow(df) > 100)          # Keep epochs with >100 samples
epochs(df -> df.interval[1] < 1.0)    # Keep epochs by metadata value
```
"""
epochs() = x -> fill(true, length(x))  # Default: select all epochs given
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> [i in epoch_numbers for i = 1:length(x)]
epochs(epoch_number::Int) = epochs([epoch_number])
epochs(epoch_ranges::Vector{UnitRange{Int}}) = epochs(union(epoch_ranges...))
epochs(epoch_numbers::Int...) = epochs(collect(epoch_numbers))
epochs(predicate::Function) = x -> map(predicate, x)  # Allow custom function predicates over dataframes

"""
    epochs_not(index::Int)
    epochs_not(indices::Int...)
    epochs_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified epochs. Operates inversely to `epochs`.

# Examples
```julia
epochs_not(1)                         # Exclude epoch 1
epochs_not(1:10)                      # Exclude first 10 epochs
```
"""
epochs_not(args...) = _negate_selector(epochs(args...))

"""
    participants()
    participants(id::Int)
    participants(ids::Int...)
    participants(ids::Union{Vector{Int}, UnitRange})
    participants(predicate::Function)

Predicate generators for participant ID selection in batch operations.

# Examples
```julia
participants()                        # Select all participants (default)
participants(1, 2, 3)                 # Select participants 1, 2, and 3
participants(1:10)                    # Select participants 1 through 10
```
"""
participants() = x -> fill(true, length(x))  # Default: select all participants given
participants(participant_ids::AbstractVector{Int}) = x -> [id in participant_ids for id in x]
participants(participant_ids::UnitRange{Int}) = participants(collect(Int, participant_ids))
participants(participant_id::Int) = participants([participant_id])
participants(participant_ranges::Vector{UnitRange{Int}}) = participants(collect(Int, union(participant_ranges...)))
participants(participant_ids::Int...) = participants(collect(Int, participant_ids))
participants(predicate::Function) = predicate  # Allow custom function predicates

"""
    participants_not(id::Int)
    participants_not(ids::Int...)
    participants_not(ids::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified participant IDs.

# Examples
```julia
participants_not(5)                   # Exclude participant 5
participants_not(1:3)                 # Exclude participants 1, 2, and 3
```
"""
participants_not(args...) = _negate_selector(participants(args...))

"""
    conditions()
    conditions(index::Int)
    conditions(indices::Int...)
    conditions(indices::Union{Vector{Int}, UnitRange})
    conditions(name::String)
    conditions(names::Vector{String})
    conditions(predicate::Function)

Predicate generators for condition selection by index or name.
Used with the `condition_selection` keyword argument on `Vector{ErpData}`,
`Vector{EpochData}`, and similar multi-condition containers.

# Examples
```julia
conditions()                          # Select all conditions (default)
conditions(1, 3)                      # Select conditions 1 and 3 by index
conditions("congruent")               # Select condition by name
conditions(["congruent", "neutral"])   # Select multiple conditions by name
```
"""
conditions() = x -> fill(true, length(x))  # Default: select all conditions given
conditions(condition_indices::AbstractVector{Int}) = x -> intersect(condition_indices, 1:length(x))
conditions(condition_indices::UnitRange{Int}) = conditions(collect(Int, condition_indices))
conditions(condition_index::Int) = conditions([condition_index])
conditions(condition_ranges::Vector{UnitRange{Int}}) = conditions(collect(Int, union(condition_ranges...)))
conditions(condition_indices::Int...) = conditions(collect(Int, condition_indices))
conditions(condition_names::Vector{String}) =
    x -> Vector{Int}(filter(!isnothing, [findfirst(dat -> condition_name(dat) == n, x) for n in condition_names]))
conditions(cond_name::String) = conditions([cond_name])
conditions(cond_names::String...) = conditions(collect(cond_names))
conditions(predicate::Function) = predicate  # Allow custom function predicates

"""
    conditions_not(index::Int)
    conditions_not(indices::Int...)
    conditions_not(indices::Union{Vector{Int}, UnitRange})
    conditions_not(name::String)
    conditions_not(names::String...)
    conditions_not(names::Vector{String})

Predicate generators that **exclude** the specified conditions.

# Examples
```julia
conditions_not(1)                     # Exclude condition 1
conditions_not("incongruent")         # Exclude condition by name
```
"""
conditions_not(args...) = _negate_selector(conditions(args...))

# === SELECTOR CONVERSION HELPERS ===
# Convert symbols, strings, numbers, vectors, ranges, or colons directly to predicates
_to_sample_func(::Nothing) = samples()
_to_sample_func(::Colon) = samples()
_to_sample_func(f::Function) = f === Colon() ? samples() : f
_to_sample_func(col::Symbol) = samples(col)
_to_sample_func(idx::Integer) = x -> [i == idx for i = 1:nrow(x)]
_to_sample_func(idx) = x -> [i in idx for i = 1:nrow(x)]

_to_channel_func(::Nothing) = channels()
_to_channel_func(::Colon) = channels()
_to_channel_func(f::Function) = f === Colon() ? channels() : f
_to_channel_func(ch::Symbol) = channels(ch)
_to_channel_func(chans::AbstractVector{Symbol}) = channels(collect(chans))
_to_channel_func(idx::AbstractString) = channels(Symbol(idx))
_to_channel_func(idx::AbstractVector{<:AbstractString}) = channels(Symbol.(idx))
_to_channel_func(idx::Integer) = channels(Int(idx))
_to_channel_func(idx::AbstractVector{<:Integer}) = channels(collect(Int, idx))
_to_channel_func(idx::AbstractRange) = channels(collect(Int, idx))
_to_channel_func(idx) = channels(idx)

_to_epoch_func(::Nothing) = epochs()
_to_epoch_func(::Colon) = epochs()
_to_epoch_func(f::Function) = f === Colon() ? epochs() : f
_to_epoch_func(idx::Integer) = epochs(Int(idx))
_to_epoch_func(idx::AbstractVector{<:Integer}) = epochs(collect(Int, idx))
_to_epoch_func(idx::AbstractRange) = epochs(collect(Int, idx))
_to_epoch_func(idx) = epochs(idx)

_to_condition_func(::Nothing) = conditions()
_to_condition_func(::Colon) = conditions()
_to_condition_func(f::Function) = f === Colon() ? conditions() : f
_to_condition_func(idx::Integer) = conditions(Int(idx))
_to_condition_func(idx::AbstractVector{<:Integer}) = conditions(collect(Int, idx))
_to_condition_func(idx::AbstractRange) = conditions(collect(Int, idx))
_to_condition_func(name::AbstractString) = conditions(String(name))
_to_condition_func(names::AbstractVector{<:AbstractString}) = conditions(collect(String, names))
_to_condition_func(idx) = conditions(idx)

_to_participant_func(::Nothing) = participants()
_to_participant_func(::Colon) = participants()
_to_participant_func(f::Function) = f === Colon() ? participants() : f
_to_participant_func(idx::Integer) = participants(Int(idx))
_to_participant_func(idx::AbstractVector{<:Integer}) = participants(collect(Int, idx))
_to_participant_func(idx::AbstractRange) = participants(collect(Int, idx))
_to_participant_func(idx) = participants(idx)

"""Preserve user-specified channel name ordering if all names are in the selection."""
function _handle_channel_names_order(user_order::Vector{Symbol}, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    # Validate: check for missing channels and duplicates
    seen = Set{Symbol}()
    for ch in user_order
        if ch ∉ selectable_cols
            @minimal_warning "Channel $(ch) not found in data!"
        elseif ch ∈ seen
            @minimal_warning "Channel $(ch) already specified!"
        else
            push!(seen, ch)
        end
    end

    # Check if this is channels([...]) or channels_not([...])
    existing_in_order = [ch for ch in user_order if ch in selectable_cols]
    existing_in_selected = [ch for ch in existing_in_order if ch in selected]

    # If all existing channels are in selected, it's channels([...]) - preserve order
    if !isempty(existing_in_order) && length(existing_in_selected) == length(existing_in_order)
        result = Symbol[]
        for ch in user_order
            if ch ∉ selectable_cols || ch ∈ result
                continue
            end
            push!(result, ch)
        end
        return result
    end
    return nothing  # Use default order
end

"""Preserve user-specified channel number ordering if all indices are in the selection."""
function _handle_channel_numbers_order(user_order_numbers, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
    # Validate: check for invalid indices and duplicates
    seen = Set{Int}()
    for i in user_order_numbers
        if i < 1 || i > length(selectable_cols)
            @minimal_warning "Channel index $(i) out of range (valid: 1:$(length(selectable_cols)))!"
        elseif i ∈ seen
            @minimal_warning "Channel index $(i) already specified!"
        else
            push!(seen, i)
        end
    end

    # Check if this is channels([...]) or channels_not([...])
    valid_indices = [i for i in user_order_numbers if 1 <= i <= length(selectable_cols)]
    selected_indices = [i for i in valid_indices if selectable_cols[i] in selected]

    # If all valid indices are in selected, it's channels([...]) - preserve order
    if !isempty(valid_indices) && length(selected_indices) == length(valid_indices)
        result = Symbol[]
        for i in user_order_numbers
            if i < 1 || i > length(selectable_cols) || selectable_cols[i] ∈ result
                continue
            end
            push!(result, selectable_cols[i])
        end
        return result
    end
    return nothing  # Use default order
end

get_selected_channels(dat, ::Colon; include_meta::Bool = true, include_extra::Bool = true) =
    get_selected_channels(dat, channels(); include_meta = include_meta, include_extra = include_extra)

get_selected_channels(dat, channel_selection; include_meta::Bool = true, include_extra::Bool = true) =
    get_selected_channels(dat, _to_channel_func(channel_selection); include_meta = include_meta, include_extra = include_extra)

function get_selected_channels(dat, channel_selection::Function; include_meta::Bool = true, include_extra::Bool = true)
    channel_selection === Colon() && (channel_selection = channels())
    # Columns/channels in dataframe to include
    metadata_cols = include_meta ? meta_labels(dat) : Symbol[]
    selectable_cols = include_extra ? vcat(channel_labels(dat), extra_labels(dat)) : channel_labels(dat)

    # Apply channel selection to non-metadata columns
    selection_mask = channel_selection(selectable_cols)
    selected = selectable_cols[selection_mask]

    # Preserve user-specified order if available
    selection_type = typeof(channel_selection)
    if hasfield(selection_type, :channel_names)
        user_order_names = getfield(channel_selection, :channel_names)
        ordered = _handle_channel_names_order(user_order_names, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    elseif hasfield(selection_type, :channel_numbers)
        user_order_numbers = getfield(channel_selection, :channel_numbers)
        ordered = _handle_channel_numbers_order(user_order_numbers, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    end

    return vcat(metadata_cols, selected)
end


"""Apply a component predicate and return matching indices."""
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    component_selection === Colon() && return collect(1:length(ica_result.ica_label))
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end
get_selected_components(ica_result::InfoIca, ::Union{Nothing,Colon}) = collect(1:length(ica_result.ica_label))

"""Apply a sample predicate and return matching row indices."""
function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    sample_selection === Colon() && return 1:nrow(dat.data)
    return findall(sample_selection(dat.data))
end
function get_selected_samples(dat::SingleDataFrameEeg, ::Union{Nothing,AllSelection,Colon})
    return 1:nrow(dat.data)
end

# Helper to select samples based on a predicate
function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    sample_selection === Colon() && return 1:nrow(dat.data[1])
    return findall(sample_selection(dat.data[1])) # assume all data have the same samples
end
function get_selected_samples(dat::MultiDataFrameEeg, ::Union{Nothing,AllSelection,Colon})
    return 1:nrow(dat.data[1])
end

# Helper to select samples from a DataFrame
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    sample_selection === Colon() && return 1:nrow(dat)
    return findall(sample_selection(dat))
end
function get_selected_samples(dat::DataFrame, ::Union{Nothing,AllSelection,Colon})
    return 1:nrow(dat)
end

# Helper to convert TimeSelection and Tuple into sample functions
function get_selected_samples(dat::Union{SingleDataFrameEeg,MultiDataFrameEeg,DataFrame}, sel::Union{TimeSelection,Tuple{Real,Real}})
    func = _interval_to_samples(sel)
    return get_selected_samples(dat, func)
end

get_selected_samples(dat::Union{SingleDataFrameEeg,MultiDataFrameEeg,DataFrame}, sample_selection) =
    get_selected_samples(dat, _to_sample_func(sample_selection))

"""
    get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection) -> Vector{Int}

Get the indices of epochs that match the epoch selection predicate.

The predicate receives `dat.data` (a `Vector{DataFrame}`) and can return a boolean vector
or indices. This enables both index-based and content-based filtering.

# Examples
```julia
# Get all epochs
selected = get_selected_epochs(dat, epochs())

# Get specific epochs by index
selected = get_selected_epochs(dat, epochs(1:10))

# Content-based: keep only epochs where interval < 1 second
selected = get_selected_epochs(dat, epochs(df -> df.interval[1] < 1.0))
```
"""
get_selected_epochs(dat::MultiDataFrameEeg, ::Union{Nothing,Colon}) = collect(1:length(dat.data))
get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection) = get_selected_epochs(dat, _to_epoch_func(epoch_selection))

function get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection::Function)
    epoch_selection === Colon() && return collect(1:length(dat.data))
    # The dat.data accessor maps to data_power for TimeFreqEpochData
    res = epoch_selection(dat.data)
    return res isa AbstractVector{Bool} ? findall(res) : res
end

"""Apply a condition predicate and return matching dataset indices."""
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

Extract a single epoch from multi-epoch EEG data and return it as ContinuousData.

# Arguments
- `dat::MultiDataFrameEeg`: The multi-DataFrame EEG data
- `epoch_idx::Int`: Index of the epoch to extract (1-based)

# Returns
- `ContinuousData`: Single DataFrame containing only the specified epoch

# Examples
```julia
# Extract epoch 3 as continuous data
single_dat = epoch_to_continuous(dat, 3)

# Now you can use single DataFrame functions
summary = channel_summary(single_dat)
```
"""
function epoch_to_continuous(dat::MultiDataFrameEeg, epoch_idx::Int)
    # Validate epoch index
    if epoch_idx < 1 || epoch_idx > length(dat.data)
        @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    end
    return ContinuousData(dat.file, dat.data[epoch_idx], dat.layout, dat.sample_rate, dat.analysis_info)
end
