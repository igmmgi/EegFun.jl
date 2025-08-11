

# channel_summary
function _channel_summary_impl(
    data::DataFrame,
    sample_selection::Vector{Int},
    channel_selection::Vector{Symbol},
)::DataFrame
    selected_data = @view data[sample_selection, channel_selection]

    # Get base statistics from describe
    stats = describe(selected_data, :min, :max, :std)

    # Add our custom columns
    stats.range = stats.max .- stats.min
    stats.var = var.(eachcol(selected_data))
    stats.zvar = zscore(stats.var)

    # Rename the variable column to channel
    rename!(stats, :variable => :channel)

    return stats
end


"""
    channel_summary(dat::ContinuousData; sample_selection::Function = samples(), channel_selection::Function = channels(), include_extr ::Bool = false)::DataFrame

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
    selected_channels =
        get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    return _channel_summary_impl(dat.data, selected_samples, selected_channels)
end



function channel_summary(
    dat::MultiDataFrameEeg;
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Process each epoch and collect results
    results = DataFrame[]

    for epoch_df in dat.data
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

    # Combine all results
    return vcat(results...)
end
