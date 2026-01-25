
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
function _channel_summary_impl(data::DataFrame, sample_selection::Vector{Int}, channel_selection::Vector{Symbol})::DataFrame
    # Input validation
    isempty(sample_selection) && @minimal_error_throw("No samples selected for channel summary")
    isempty(channel_selection) && @minimal_error_throw("No channels selected for channel summary")

    # Check that all selected channels exist in data
    missing_channels = setdiff(channel_selection, propertynames(data))
    !isempty(missing_channels) && @minimal_error_throw("Channels not found in data: $(missing_channels)")

    # Check that all sample indices are valid
    invalid_samples = sample_selection[(sample_selection.<1).|(sample_selection.>nrow(data))]
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

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)
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

        # Epoch number is now derived from vector index
        original_epoch_number = epoch_idx

        # Create ContinuousData from this epoch DataFrame  
        single_dat = ContinuousData(dat.file, epoch_df, dat.layout, dat.sample_rate, dat.analysis_info)

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


"""
Batch channel summary statistics for EEG/ERP data.
"""

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC HELPERS
=============================================================================#

"""Generate default output directory name for channel summary."""
function _default_channel_summary_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "channel_summary_$(pattern)")
end

#=============================================================================
    CHANNEL-SUMMARY-SPECIFIC PROCESSING
=============================================================================#

"""
Process a single file through channel summary pipeline.
Returns tuple of (BatchResult, Vector{DataFrame}) with all condition results.
"""
function _process_channel_summary_file(
    filepath::String,
    condition_selection::Function,
    sample_selection::Function,
    channel_selection::Function,
    include_extra::Bool,
)
    filename = basename(filepath)

    # Load data
    data_var = load_data(filepath)
    if isnothing(data_var)
        return (BatchResult(false, filename, "No recognized data variable"), DataFrame[])
    end

    # Get condition count before selection
    n_conditions = length(data_var)

    # Select conditions
    data_var = _condition_select(data_var, condition_selection)

    # Determine actual condition numbers for tracking
    # After _condition_select, data_var is filtered but we need original condition numbers
    condition_numbers = 1:length(data_var)

    # Process each condition and collect results
    summary_dfs = DataFrame[]
    for (cond_idx, data) in enumerate(data_var)
        condition = condition_numbers[cond_idx]

        # Compute channel summary
        summary_df =
            channel_summary(data; sample_selection = sample_selection, channel_selection = channel_selection, include_extra = include_extra)

        # Add metadata columns
        insertcols!(summary_df, 1, :file => splitext(filename)[1])
        insertcols!(summary_df, 2, :condition => condition)

        push!(summary_dfs, summary_df)
    end

    n_conditions = length(summary_dfs)
    return (BatchResult(true, filename, "Processed $n_conditions condition(s)"), summary_dfs)
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    channel_summary(file_pattern::String; 
                    input_dir::String = pwd(), 
                    participant_selection::Function = participants(),
                    condition_selection::Function = conditions(),
                    sample_selection::Function = samples(),
                    channel_selection::Function = channels(),
                    include_extra::Bool = false,
                    output_dir::Union{String, Nothing} = nothing,
                    output_file::String = "channel_summary")

Batch process EEG/ERP data files to compute channel summary statistics.

This function loads JLD2 files, computes channel summary statistics for each file,
and saves the results to CSV files.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "epochs", "erps", "cleaned")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `sample_selection::Function`: Function for sample filtering (default: all samples)
- `channel_selection::Function`: Function for channel filtering (default: all channels)
- `include_extra::Bool`: Whether to include extra channels (default: false)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `output_file::String`: Base name for output files (default: "channel_summary")

# Examples
```julia
# Compute channel summary for all epoch files
channel_summary("epochs")

# Process specific participants and conditions
channel_summary("erps_cleaned", participants=[1, 2, 3], conditions=[1, 2])

# Compute summary for specific channels only
channel_summary("epochs", channel_selection=channels([:Fp1, :Fp2, :F3, :F4]))

# Include extra channels (EOG, etc.)
channel_summary("epochs", include_extra=true)
```
"""
function channel_summary(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    sample_selection::Function = samples(),
    channel_selection::Function = channels(),
    include_extra::Bool = false,
    output_dir::Union{String,Nothing} = nothing,
    output_file::String = "channel_summary",
)

    # Setup logging
    log_file = "$(output_file).log"
    setup_global_logging(log_file)

    try
        @info "Batch channel summary started at $(now())"
        @log_call "channel_summary"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_channel_summary_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files to process"

        # Process all files and collect DataFrames
        all_summaries = DataFrame[]
        n_success = 0
        n_error = 0

        for (i, file) in enumerate(files)
            @info "Channel summary: $file ($i/$(length(files)))"

            input_path = joinpath(input_dir, file)

            result, summary_dfs = try
                _process_channel_summary_file(input_path, condition_selection, sample_selection, channel_selection, include_extra)
            catch e
                @error "Error processing $file" exception = (e, catch_backtrace())
                (BatchResult(false, file, "Exception: $(sprint(showerror, e))"), DataFrame[])
            end

            # Log result
            if result.success
                @info "  ✓ $(result.message)"
                append!(all_summaries, summary_dfs)
                n_success += 1
            else
                @minimal_warning "  ✗ $(result.message)"
                n_error += 1
            end
        end

        # Save combined results to single CSV
        if !isempty(all_summaries)
            combined_df = vcat(all_summaries...)
            output_path = joinpath(output_dir, "$(output_file).csv")
            CSV.write(output_path, combined_df)
            @info "Combined results saved to: $output_path"
        else
            @minimal_warning "No results to save"
        end

        @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
        @info "Output saved to: $output_dir"

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
