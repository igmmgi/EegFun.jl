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
_interval_to_samples(f::Function) = f

"""
    _combine_interval_sample(interval::Interval, sample_selection::Function)

Combine an interval selection with a sample selection predicate.
"""
function _combine_interval_sample(interval::Interval, sample_selection::Function)
    if interval isa Nothing || interval isa AllSelection
        return sample_selection
    else
        interval_sel = _interval_to_samples(interval)
        return x -> interval_sel(x) .& sample_selection(x)
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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
- `epoch_selection::Function`: Epoch predicate. You can supply numerical indices using `epochs(1:5)` or a custom DataFrame predicate via `epochs(df -> nrow(df) > 100)` to efficiently filter trials by their contents (default: `epochs()` for all).
- `channel_selection::Function`: Channel predicate (default: `channels()` for all)
- `sample_selection::Function`: Sample predicate (default: `samples()` for all)
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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
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
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
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
    datasets::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
    datasets::Vector{EpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
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
    datasets::Vector{TimeFreqData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
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
    datasets::Vector{TimeFreqEpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
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
