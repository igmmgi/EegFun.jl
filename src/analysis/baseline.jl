"""
    _apply_baseline!(dat::DataFrame, channels, baseline_interval)

Internal function that applies baseline correction to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to baseline correct
- `channels::Vector{Symbol}`: Names of channels to correct
- `baseline_interval::Union{IntervalIndex,IntervalTime}`: Time interval for baseline calculation

# Notes
- Subtracts the mean of the baseline interval from each specified channel
- Modifies the input DataFrame in-place
"""
function _apply_baseline!(
    dat::DataFrame,
    channels::Vector{Symbol},
    baseline_interval::Union{IntervalIndex,IntervalTime},
)
    # Compute mean baseline interval and apply to each channel
    baseline_means = mean.(eachcol(dat[baseline_interval.start:baseline_interval.stop, channels]))
    for (channel, mean_val) in zip(channels, baseline_means)
        dat[!, channel] .-= mean_val
    end
end

"""
    _apply_baseline!(dat::Vector{DataFrame}, channels, baseline_interval)

Internal function that applies baseline correction to each DataFrame in a vector using broadcasting.
"""
function _apply_baseline!(
    dat::Vector{DataFrame},
    channels::Vector{Symbol},
    baseline_interval::Union{IntervalIndex,IntervalTime},
)
    _apply_baseline!.(dat, Ref(channels), Ref(baseline_interval))
end

"""
    baseline!(dat::EegData, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to EEG data.

# Arguments
- `dat::EegData`: The data to baseline correct
- `baseline_interval::Union{IntervalIndex,IntervalTime}`: Time/index interval for baseline calculation
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies the input data in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
"""
function baseline!(
    dat::EegData,
    baseline_interval::Union{AbstractInterval,Tuple{Real,Real}};
    channel_selection::Function = channels(),
)
    # Validate baseline interval
    baseline_interval = validate_baseline_interval(dat, baseline_interval)

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for baseline correction"
        return
    end

    # Apply baseline correction (dispatch handles DataFrame vs Vector{DataFrame})
    @info "Applying baseline correction to $(length(selected_channels)) channels over interval: $(baseline_interval.start) to $(baseline_interval.stop)"
    _apply_baseline!(dat.data, selected_channels, baseline_interval)
    return nothing
end

"""
    baseline!(dat::EegData; channel_selection=channels())

Apply baseline correction in-place to EEG data using the entire time range.

# Arguments
- `dat::EegData`: The data to baseline correct
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies the input data in-place by subtracting the baseline mean
- Uses the entire time range for baseline calculation
"""
function baseline!(dat::EegData; channel_selection::Function = channels())
    baseline!(dat, IntervalIndex(start = 1, stop = n_samples(dat)); channel_selection = channel_selection)
    return nothing
end

"""
    baseline!(dat::Vector{EpochData}, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to a vector of EpochData objects.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data to baseline correct
- `baseline_interval::Union{IntervalIndex,IntervalTime}`: Time/index interval for baseline calculation
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
"""
function baseline!(
    dat::Vector{EpochData},
    baseline_interval::Union{AbstractInterval,Tuple{Real,Real}};
    channel_selection::Function = channels(),
)
    baseline!.(dat, Ref(baseline_interval); channel_selection = channel_selection)
    return nothing
end

"""
    baseline!(dat::Vector{EpochData}; channel_selection=channels())

Apply baseline correction in-place to a vector of EpochData objects using the entire time range.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data to baseline correct
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by subtracting the baseline mean
- Uses the entire time range for baseline calculation
"""
function baseline!(dat::Vector{EpochData}; channel_selection::Function = channels())
    baseline!.(dat; channel_selection = channel_selection)
    return nothing
end

# generates all non-mutating versions
@add_nonmutating baseline!


#=============================================================================
    BASELINE-SPECIFIC HELPERS
=============================================================================#

# TODO: actually, these are probably not needed as baseline can just be applied inplace whenever it is needed
# For now, juse keep for a little bit of consistency with other functions

"""Generate default output directory name for baseline operation."""
function _default_baseline_output_dir(input_dir::String, pattern::String, baseline_interval::Union{IntervalIndex,IntervalTime})
    interval_str = "$(baseline_interval.start)_to_$(baseline_interval.stop)"
    joinpath(input_dir, "baseline_$(pattern)_$(interval_str)")
end


"""
Process a single file through baseline correction pipeline.
Returns BatchResult with success/failure info.
"""
function _process_baseline_file(
    filepath::String,
    output_path::String,
    baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real}},
    conditions,
)
    filename = basename(filepath)

    # Load data
    data = load_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, conditions)

    # Apply baseline correction (mutates data in-place)
    baseline!.(data, Ref(baseline_interval))

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = data)

    return BatchResult(true, filename, "Baseline corrected successfully")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    baseline(file_pattern::String, baseline_interval; 
             input_dir::String = pwd(), 
             participants::Union{Int, Vector{Int}, Nothing} = nothing,
             conditions::Union{Int, Vector{Int}, Nothing} = nothing,
             output_dir::Union{String, Nothing} = nothing)

Apply baseline correction to EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `baseline_interval`: Baseline interval as:
  - `Tuple{Real,Real}`: Time interval in seconds (e.g., `(-0.2, 0.0)`)
  - `IntervalTime`: Time interval struct (e.g., `IntervalTime(start=-0.2, stop=0.0)`)
  - `IntervalIndex`: Index interval struct (e.g., `IntervalIndex(start=1, stop=51)`)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on baseline settings)

# Example
```julia
# Baseline correct all epochs from -0.2 to 0.0 seconds
baseline("epochs", (-0.2, 0.0))

# Baseline correct specific participant
baseline("epochs", (-0.2, 0.0), participants=3)

# Baseline correct specific participants and conditions
baseline("epochs", (-0.2, 0.0), participants=[3, 4], conditions=[1, 2])

# Use IntervalTime explicitly
baseline("epochs", IntervalTime(start=-0.2, stop=0.0))
```
"""
function baseline(
    file_pattern::String,
    baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real}};
    input_dir::String = pwd(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    conditions::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "baseline.log"
    setup_global_logging(log_file)

    try
        @info ""
        @info "Batch baseline correction started at $(now())"
        @log_call "baseline"

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Validate and normalize baseline interval (converts tuple to IntervalTime, validates structure)
        baseline_interval = validate_baseline_interval(baseline_interval)

        # Setup directories
        output_dir = something(
            output_dir,
            _default_baseline_output_dir(input_dir, file_pattern, baseline_interval),
        )
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        @info "Baseline interval: $baseline_interval"
        process_fn = (input_path, output_path) -> _process_baseline_file(
            input_path,
            output_path,
            baseline_interval,
            conditions,
        )

        # Execute batch operation
        results = _run_batch_operation(
            process_fn,
            files,
            input_dir,
            output_dir;
            operation_name = "Baseline correction",
        )

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
