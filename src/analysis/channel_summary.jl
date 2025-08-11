
# ============================================================================ #
#                           CHANNEL SUMMARY FUNCTIONS                         #
# ============================================================================ #

"""
    _channel_summary_impl(data::DataFrame, sample_selection::Vector{Int}, channel_selection::Vector{Symbol})::DataFrame

Internal implementation for computing channel summary statistics.

# Arguments
- `data::DataFrame`: The data frame containing EEG data
- `sample_selection::Vector{Int}`: Indices of samples to include
- `channel_selection::Vector{Symbol}`: Names of channels to include

# Returns
- `DataFrame`: Summary statistics with columns: channel, min, max, std, range, var, zvar

# Statistics Computed
- `min`: Minimum value per channel
- `max`: Maximum value per channel  
- `std`: Standard deviation per channel
- `range`: Range (max - min) per channel
- `var`: Variance per channel
- `zvar`: Z-scored variance (relative to other channels)
"""
function _channel_summary_impl(
    data::DataFrame,
    sample_selection::Vector{Int},
    channel_selection::Vector{Symbol},
)::DataFrame
    # Input validation
    isempty(sample_selection) && @minimal_error_throw("No samples selected for channel summary")
    isempty(channel_selection) && @minimal_error_throw("No channels selected for channel summary")
    
    # Check that all selected channels exist in data
    missing_channels = setdiff(channel_selection, propertynames(data))
    !isempty(missing_channels) && @minimal_error_throw("Channels not found in data: $(missing_channels)")
    
    # Check that all sample indices are valid
    invalid_samples = sample_selection[sample_selection .< 1 .|| sample_selection .> nrow(data)]
    !isempty(invalid_samples) && @minimal_error_throw("Invalid sample indices: $(invalid_samples)")
    
    selected_data = @view data[sample_selection, channel_selection]

    # Get base statistics from describe
    stats = describe(selected_data, :min, :max, :std)

    # Add our custom columns
    stats.range = stats.max .- stats.min
    stats.var = var.(eachcol(selected_data))
    
    # Handle case where all channels have zero variance (avoid NaN in zscore)
    if all(stats.var .== 0.0)
        stats.zvar = zeros(length(stats.var))
    else
        stats.zvar = zscore(stats.var)
    end

    # Rename the variable column to channel
    rename!(stats, :variable => :channel)

    return stats
end

# ============================================================================ #
#                      SINGLE DATAFRAME EEG CHANNEL SUMMARY                   #
# ============================================================================ #

"""
    channel_summary(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_extra::Bool = false)::DataFrame

Computes summary statistics for EEG channels.

# Arguments
- `dat::ContinuousData`: The ContinuousData object containing EEG data.
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples).
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels).
- `include_extra::Bool`: Whether to include additional channels (default: false).

# Returns
A DataFrame containing summary statistics for each channel.

# Examples

## Basic Usage
```julia
# Channel summary for layout channels only (default)
summary = channel_summary(dat)

# Channel summary for specific layout channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Channel summary excluding reference channels from layout
summary = channel_summary(dat, channel_selection = channels_not([:M1, :M2]))
```

## Including Additional Channels
```julia
# Channel summary for additional channels (EOG, extreme value flags, etc.)
# The function automatically detects when you specify additional channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```

## Channel Selection
```julia
# Summary for specific channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :F3, :F4]))

# Summary excluding reference channels
summary = channel_summary(dat, channel_selection = channels_not([:M1, :M2]))

# Summary for frontal channels only (channels 1-10)
summary = channel_summary(dat, channel_selection = channels(1:10))

# Summary for channels starting with "F" (frontal)
summary = channel_summary(dat, channel_selection = x -> startswith.(string.(x), "F"))
```

## Sample Selection
```julia
# Exclude extreme values
summary = channel_summary(dat, sample_selection = samples_not(:is_extreme_value_100))

# Exclude multiple types of bad samples
summary = channel_summary(dat, sample_selection = samples_or_not([:is_extreme_value_100, :is_vEOG, :is_hEOG]))

# Only include samples within epoch windows
summary = channel_summary(dat, sample_selection = samples(:epoch_window))

# Include samples that are both in epoch window AND not extreme
summary = channel_summary(dat, sample_selection = samples_and([:epoch_window, samples_not(:is_extreme_value_100)]))
```

## Combined Selection
```julia
# Exclude reference channels and extreme values
summary = channel_summary(dat, 
    channel_selection = channels_not([:M1, :M2]),
    sample_selection = samples_not(:is_extreme_value_100)
)

# Only frontal channels, exclude bad samples
summary = channel_summary(dat, 
    channel_selection = channels(1:10),
    sample_selection = samples_or_not([:is_extreme_value_100, :is_vEOG])
)

# Complex filtering: frontal channels, good samples, within epochs
summary = channel_summary(dat, 
    channel_selection = channels([:Fp1, :Fp2, :F3, :F4, :F5, :F6, :F7, :F8]),
    sample_selection = samples_and([
        :epoch_window, 
        samples_not(:is_extreme_value_100),
        samples_not(:is_vEOG),
        samples_not(:is_hEOG)
    ])
)
```

## Additional Channels (not in layout)
```julia
# Include derived channels like EOG
# The function automatically switches to all available channels when needed
summary = channel_summary(dat, channel_selection = channels([:vEOG, :hEOG]))

# Mix layout channels and additional channels
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2, :vEOG, :hEOG]))
```
"""
function channel_summary(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Input validation
    nrow(dat.data) == 0 && @minimal_error_throw("Cannot compute channel summary: data is empty")
    
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    
    return _channel_summary_impl(dat.data, selected_samples, selected_channels)
end

# ============================================================================ #
#                       MULTI DATAFRAME EEG CHANNEL SUMMARY                   #
# ============================================================================ #

"""
    channel_summary(dat::MultiDataFrameEeg; sample_selection::Function = samples(), channel_selection::Function = channels(), include_meta::Bool = false, include_extra::Bool = false)::DataFrame

Computes summary statistics for EEG channels across multiple epochs.

# Arguments
- `dat::MultiDataFrameEeg`: The MultiDataFrameEeg object containing epoch data (e.g., EpochData)
- `sample_selection::Function`: Function that returns boolean vector for sample filtering (default: include all samples)
- `channel_selection::Function`: Function that returns boolean vector for channel filtering (default: include all channels)
- `include_meta::Bool`: Whether to include metadata columns (default: false)
- `include_extra::Bool`: Whether to include additional channels (default: false)

# Returns
A DataFrame containing summary statistics for each channel in each epoch, with an additional `epoch` column.

# Examples
```julia
# Basic epoch-wise channel summary
summary = channel_summary(epoch_data)

# Summary for specific channels across epochs
summary = channel_summary(epoch_data, channel_selection = channels([:Fp1, :Fp2]))

# Summary excluding bad samples
summary = channel_summary(epoch_data, sample_selection = samples_not(:is_bad))
```
"""
function channel_summary(
    dat::MultiDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Input validation
    isempty(dat.data) && @minimal_error_throw("Cannot compute channel summary: no epochs in data")
    
    # Process each epoch and collect results
    results = DataFrame[]

    for (epoch_idx, epoch_df) in enumerate(dat.data)
        # Input validation for this epoch
        if nrow(epoch_df) == 0
            @minimal_warning("Skipping empty epoch $(epoch_idx)")
            continue
        end
        
        # Get the original epoch number from the data
        original_epoch_number = epoch_df.epoch[1]  # All rows in an epoch have the same epoch number

        # Create ContinuousData from this epoch DataFrame  
        single_dat = ContinuousData(epoch_df, dat.layout, dat.sample_rate, dat.analysis_info)

        # Get summary for this epoch
        epoch_summary = channel_summary(
            single_dat;
            sample_selection = sample_selection,
            channel_selection = channel_selection,
            include_meta = include_meta,
            include_extra = include_extra,
        )

        # Add epoch column as first column with original epoch number
        insertcols!(epoch_summary, 1, :epoch => fill(original_epoch_number, nrow(epoch_summary)))

        push!(results, epoch_summary)
    end

    # Check if we have any results
    isempty(results) && @minimal_error_throw("No valid epochs found for channel summary")
    
    # Combine all results
    return vcat(results...)
end
