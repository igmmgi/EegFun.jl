"""
    _apply_baseline!(dat::DataFrame, channels, baseline_interval)

Internal function that applies baseline correction to specified channels in a DataFrame.

# Arguments
- `dat::DataFrame`: The data to baseline correct
- `channels::Vector{Symbol}`: Names of channels to correct
- `baseline_interval::Tuple{Int,Int}`: Time interval for baseline calculation

# Notes
- Subtracts the mean of the baseline interval from each specified channel
- Modifies the input DataFrame in-place
"""
function _apply_baseline!(dat::DataFrame, channels::Vector{Symbol}, baseline_interval::Tuple{Int,Int})
    # Compute mean baseline interval and apply to each channel
    baseline_means = mean.(eachcol(dat[baseline_interval[1]:baseline_interval[2], channels]))
    for (channel, mean_val) in zip(channels, baseline_means)
        dat[!, channel] .-= mean_val
    end
end

"""
    _apply_baseline!(dat::Vector{DataFrame}, channels, baseline_interval)

Internal function that applies baseline correction to each DataFrame in a vector using broadcasting.
"""
function _apply_baseline!(dat::Vector{DataFrame}, channels::Vector{Symbol}, baseline_interval::Tuple{Int,Int})
    _apply_baseline!.(dat, Ref(channels), Ref(baseline_interval))
end

"""
    baseline!(dat::EegData, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to EEG data.

# Arguments
- `dat::EegData`: The data to baseline correct
- `baseline_interval::Interval`: Time/index interval, tuple, range, or predicate function
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies the input data in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
- If a predicate function is provided, it's applied to the data to determine the baseline window
"""
function baseline!(dat::EegData, baseline_interval::Interval; channel_selection::Function = channels())
    # Validate baseline interval
    baseline_interval = _validate_baseline_interval(dat, baseline_interval)

    # Get selected channels
    selected_channels = get_selected_channels(dat, channel_selection, include_meta = false, include_extra = false)
    if isempty(selected_channels)
        @minimal_warning "No channels selected for baseline correction"
        return
    end

    # Apply baseline correction (dispatch handles DataFrame vs Vector{DataFrame})
    @info "Applying baseline correction to $(length(selected_channels)) channels over interval: $(baseline_interval[1]) to $(baseline_interval[2])"
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
    # Use entire time range for baseline
    time_vec = dat.data.time
    baseline!(dat, (first(time_vec), last(time_vec)); channel_selection = channel_selection)
    return nothing
end

"""
    baseline!(dat::Vector{EpochData}, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to a vector of EpochData objects.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data to baseline correct
- `baseline_interval::Interval`: Time/index interval for baseline calculation, or a predicate function
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
- If a predicate function is provided, it's applied to each epoch's data to determine the baseline window
"""
function baseline!(dat::Vector{EpochData}, baseline_interval::Interval; channel_selection::Function = channels())
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

"""
    baseline!(dat::Vector{ErpData}, baseline_interval; channel_selection=channels())

Apply baseline correction in-place to a vector of ErpData objects.

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data to baseline correct
- `baseline_interval::Interval`: Time/index interval for baseline calculation, or a predicate function
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each ErpData in the vector in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
- If a predicate function is provided, it's applied to each ERP's data to determine the baseline window
"""
function baseline!(dat::Vector{ErpData}, baseline_interval::Interval; channel_selection::Function = channels())
    baseline!.(dat, Ref(baseline_interval); channel_selection = channel_selection)
    return nothing
end



"""
    baseline!(dat::Vector{ErpData}; channel_selection=channels())

Apply baseline correction in-place to a vector of ErpData objects using the entire time range.

# Arguments
- `dat::Vector{ErpData}`: Vector of ERP data to baseline correct
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each ErpData in the vector in-place by subtracting the baseline mean
- Uses the entire time range for baseline calculation
"""
function baseline!(dat::Vector{ErpData}; channel_selection::Function = channels())
    baseline!.(dat; channel_selection = channel_selection)
    return nothing
end

# TODO: actually, these are probably not needed as baseline can just be applied inplace whenever it is needed
# For now, just keep for a little bit of consistency with other functions?
@add_nonmutating baseline!


#=============================================================================
    BASELINE-SPECIFIC HELPERS
=============================================================================#
"""Generate default output directory name for baseline operation."""
function _default_baseline_output_dir(input_dir::String, pattern::String, baseline_interval::Tuple{Real,Real})
    interval_str = "$(baseline_interval[1])_to_$(baseline_interval[2])"
    joinpath(input_dir, "baseline_$(pattern)_$(interval_str)")
end


# === BASELINE INTERVAL CONVERSION (Multiple Dispatch) ===

"""Helper to convert and validate baseline intervals to (start, stop) tuples"""
function _to_interval(interval::Tuple{Real,Real})
    length(interval) == 2 || @minimal_error_throw("Baseline interval tuple must have 2 elements (start, stop), got: $interval")
    start, stop = Float64(interval[1]), Float64(interval[2])
    start > stop && @minimal_error_throw("Baseline start ($start) must be <= stop ($stop)")
    return (start, stop)
end


"""
    _validate_baseline_interval(baseline_interval::Interval) -> Tuple{Real,Real}

Validate baseline interval structure (without data). Returns normalized (start, stop) tuple.
Converts ranges to tuples.

# Arguments
- `baseline_interval::Interval`: Interval to validate

# Returns
- `Tuple{Real,Real}`: Normalized (start, stop) tuple

# Throws
- `ArgumentError`: If interval structure is invalid
"""
function _validate_baseline_interval(baseline_interval::Interval)::Tuple{Real,Real}
    return _to_interval(baseline_interval)
end

"""
    _validate_baseline_interval(time::AbstractVector, baseline_interval::Interval) -> Tuple{Int,Int}

Validate and convert baseline interval to (start_idx, stop_idx) tuple with bounds checking.

# Arguments
- `time::AbstractVector`: Time points vector
- `baseline_interval::Interval`: Interval to validate

# Returns
- `Tuple{Int,Int}`: Validated (start_idx, stop_idx) tuple

# Throws
- `ArgumentError`: If interval is invalid
"""
function _validate_baseline_interval(time::AbstractVector, baseline_interval::Interval)::Tuple{Int,Int}
    # First normalize (convert tuple/range, validate structure) using multiple dispatch
    baseline_tuple = _to_interval(baseline_interval)

    # Validate that time values are within the data's time range
    time_min, time_max = first(time), last(time)
    start_time, stop_time = baseline_tuple[1], baseline_tuple[2]

    if start_time < time_min || stop_time > time_max
        @minimal_error_throw "Invalid baseline_interval: $(baseline_tuple) - time values must be within data range ($(time_min) to $(time_max))"
    end

    # Convert time tuple to index tuple
    result = find_idx_start_end(time, start_time, stop_time)

    # Check if interval was found
    if isnothing(result)
        @minimal_error_throw "Invalid baseline_interval: $(baseline_tuple) - interval not found in time vector (range: $(time_min) to $(time_max))"
    end

    start_idx, stop_idx = result

    # Validate bounds (should always pass if range check passed, but keep for safety)
    if !(1 <= start_idx <= length(time)) || !(1 <= stop_idx <= length(time)) || !(start_idx <= stop_idx)
        @minimal_error_throw "Invalid baseline_interval: ($start_idx, $stop_idx) (time vector length: $(length(time)))"
    end

    return (start_idx, stop_idx)
end

function _validate_baseline_interval(dat::MultiDataFrameEeg, baseline_interval::Interval)::Tuple{Int,Int}
    if baseline_interval isa Function
        # assume all data have the same time, so we just use the first epoch to get indices
        mask = baseline_interval(dat.data[1])
        indices = findall(mask)
        isempty(indices) && @minimal_error_throw "Baseline window selection returned no samples!"
        return (first(indices), last(indices))
    end
    return _validate_baseline_interval(dat.data[1].time, baseline_interval) # assume all data have the same time
end

function _validate_baseline_interval(dat::SingleDataFrameEeg, baseline_interval::Interval)::Tuple{Int,Int}
    if baseline_interval isa Function
        mask = baseline_interval(dat.data)
        indices = findall(mask)
        isempty(indices) && @minimal_error_throw "Baseline window selection returned no samples!"
        # Check for contiguity (warn if not, but proceed)
        if !all(diff(indices) .== 1)
            @minimal_warning "Selected baseline window is not contiguous!"
        end
        return (first(indices), last(indices))
    end
    return _validate_baseline_interval(dat.data.time, baseline_interval)
end





"""
Process a single file through baseline correction pipeline.
Returns BatchResult with success/failure info.
"""
function _process_baseline_file(filepath::String, output_path::String, baseline_interval::Interval, condition_selection::Function)
    filename = basename(filepath)

    # Read data
    data = read_data(filepath)
    if isnothing(data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is valid EEG data (Vector of ErpData or EpochData)
    if !(data isa Vector{<:Union{ErpData,EpochData}})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData} or Vector{EpochData}")
    end

    # Select conditions
    data = _condition_select(data, condition_selection)

    # Apply baseline correction (mutates data in-place)
    baseline!.(data, Ref(baseline_interval))

    # Save (always use "data" as variable name since read_data finds by type)
    jldsave(output_path; data = data)

    return BatchResult(true, filename, "Baseline corrected successfully")
end


"""
    baseline(file_pattern::String, baseline_interval; 
             input_dir::String = pwd(), 
             participant_selection::Function = participants(),
             condition_selection::Function = conditions(),
             output_dir::Union{String, Nothing} = nothing)

Apply baseline correction to EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `baseline_interval`: Baseline interval as:
  - Tuple: `(-0.2, 0.0)` - time range in seconds
  - Range: `-0.2:0` - time range using range syntax
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on baseline settings)

# Notes
- For in-memory data with function-based baseline selection, use `baseline!(data, samples((-0.2, 0.0)))` instead

# Example
```julia
# Baseline correct all epochs from -0.2 to 0.0 seconds
baseline("epochs", (-0.2, 0.0))

# Baseline correct specific participant
baseline("epochs", (-0.2, 0.0), participant_selection=participants(3))

# Baseline correct specific participants and conditions
baseline("epochs", (-0.2, 0.0), participant_selection=participants([3, 4]), condition_selection=conditions([1, 2]))

# Use tuple directly
```
"""
function baseline(
    file_pattern::String,
    baseline_interval::Interval;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
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

        # Validate and normalize baseline interval (validates tuple structure)
        baseline_interval = _validate_baseline_interval(baseline_interval)

        # Setup directories
        output_dir = something(output_dir, _default_baseline_output_dir(input_dir, file_pattern, baseline_interval))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # Create processing function with captured parameters
        @info "Baseline interval: $baseline_interval"
        process_fn = (input_path, output_path) -> _process_baseline_file(input_path, output_path, baseline_interval, condition_selection)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Baseline correction")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end




