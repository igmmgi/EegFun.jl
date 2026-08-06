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
    isempty(sample_selection) && @minimal_error("No samples selected for channel summary")
    isempty(channel_selection) && @minimal_error("No channels selected for channel summary")

    # Check that all selected channels exist in data
    missing_channels = setdiff(channel_selection, propertynames(data))
    !isempty(missing_channels) && @minimal_error("Channels not found in data: $(missing_channels)")

    # Check that all sample indices are valid
    invalid_samples = sample_selection[(sample_selection .< 1) .| (sample_selection .> nrow(data))]
    !isempty(invalid_samples) && @minimal_error("Invalid sample indices: $(invalid_samples)")

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
    channel_summary(dat::SingleDataFrameEeg; sample_selection::Function = samples(),
    interval_selection::Interval = times(), channel_selection::Function = channels(), include_extra::Bool = false)::DataFrame
    channel_summary(dat::MultiDataFrameEeg; sample_selection::Function = samples(),
    interval_selection::Interval = times(), channel_selection::Function = channels(), include_meta::Bool = false, include_extra::Bool = false)::DataFrame
    channel_summary(file_pattern::String; input_dir::String = pwd(),
    participant_selection::Function = participants(), condition_selection::Function = conditions(),
    sample_selection::Function = samples(), channel_selection::Function = channels(),
    include_extra::Bool = false, output_dir = nothing, output_file::String = "channel_summary")

Compute summary statistics (min, max, std, range, var, zvar) for EEG channels.

The `SingleDataFrameEeg` method returns one row per channel; the `MultiDataFrameEeg` method
(e.g. `EpochData`) returns one row per channel per epoch with an `:epoch` column.
The `file_pattern` method batch-processes JLD2 files and writes results to CSV.

# Examples
```julia
# Basic usage
summary = channel_summary(dat)

# Specific channel selection
summary = channel_summary(dat, channel_selection = channels([:Fp1, :Fp2]))

# Exclude bad samples
summary = channel_summary(dat, sample_selection = samples_not(:is_extreme_value_100))

# Epoch-wise summary
summary = channel_summary(epoch_data, channel_selection = channels([:Fp1, :Fp2]))

# Batch processing
channel_summary("epochs", channel_selection = channels([:Fp1, :Fp2]))
```
"""
function channel_summary(
    dat::SingleDataFrameEeg;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Input validation
    nrow(dat.data) == 0 && @minimal_error("Cannot compute channel summary: data is empty")

    selected_channels = get_selected_channels(dat, channel_selection; include_meta = include_meta, include_extra = include_extra)

    # Combine interval and sample selection (consistent with subset() pattern)
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_samples = get_selected_samples(dat, combined_sel)

    return _channel_summary_impl(dat.data, selected_samples, selected_channels)
end

# ============================================================================ #
#                       MULTI DATAFRAME EEG CHANNEL SUMMARY                   #
# ============================================================================ #


function channel_summary(
    dat::MultiDataFrameEeg;
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    channel_selection::Function = channels(),
    include_meta::Bool = false,
    include_extra::Bool = false,
)::DataFrame
    # Input validation
    isempty(dat.data) && @minimal_error("Cannot compute channel summary: no epochs in data")

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
            interval_selection = interval_selection,
            channel_selection = channel_selection,
            include_meta = include_meta,
            include_extra = include_extra,
        )

        # Add epoch column as first column with original epoch number
        insertcols!(epoch_summary, 1, :epoch => fill(original_epoch_number, nrow(epoch_summary)))

        push!(results, epoch_summary)
    end

    # Check if we have any results
    isempty(results) && @minimal_error("No valid epochs found for channel summary")

    # Combine all results
    return vcat(results...)
end

channel_summary(dat::Vector{<:MultiDataFrameEeg}; kwargs...) = channel_summary.(dat; kwargs...)


"""
Batch channel summary statistics for EEG/ERP data.
"""

# === CHANNEL-SUMMARY-SPECIFIC HELPERS ===

"""Generate default output directory name for channel summary."""
function _default_channel_summary_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "channel_summary_$(clean_pattern(pattern))")
end

# === CHANNEL-SUMMARY-SPECIFIC PROCESSING ===

"""
Process a single file through channel summary pipeline.
Returns tuple of (BatchResult, Vector{DataFrame}) with all condition results.
"""
function _process_channel_summary_file(
    filepath::String,
    condition_selection::Function,
    sample_selection::Function,
    interval_selection::Interval,
    channel_selection::Function,
    include_extra::Bool,
)
    filename = basename(filepath)

    # Read data
    data_var = read_data(filepath)
    if isnothing(data_var)
        return (BatchResult(false, filename, "No recognized data variable"), DataFrame[])
    end

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
        summary_df = channel_summary(
            data;
            sample_selection = sample_selection,
            interval_selection = interval_selection,
            channel_selection = channel_selection,
            include_extra = include_extra,
        )

        # Add metadata columns
        insertcols!(summary_df, 1, :file => splitext(filename)[1])
        insertcols!(summary_df, 2, :condition => condition)

        push!(summary_dfs, summary_df)
    end

    n_conditions = length(summary_dfs)
    return (BatchResult(true, filename, "Processed $n_conditions condition(s)"), summary_dfs)
end

# === MAIN API FUNCTION ===


function channel_summary(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
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
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
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
                _process_channel_summary_file(
                    input_path,
                    condition_selection,
                    sample_selection,
                    interval_selection,
                    channel_selection,
                    include_extra,
                )
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
