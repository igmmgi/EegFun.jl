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
- `baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real},Function}`: Time/index interval for baseline calculation, or a predicate function (e.g., `samples((start, end))`)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies the input data in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
- If a predicate function is provided, it's applied to the data to determine the baseline window
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
    baseline!(dat::EegData, baseline_selection::Function; channel_selection=channels())

Apply baseline correction in-place to EEG data using a predicate function.

# Arguments
- `dat::EegData`: The data to baseline correct
- `baseline_selection::Function`: Sample selection predicate (e.g., `samples((start, end))`). Applied to the data DataFrame to determine baseline window.
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies the input data in-place by subtracting the baseline mean
- The predicate function is applied to `dat.data` to get a boolean mask
- If the predicate selects all samples (default), baseline correction is skipped
- If the predicate selects no samples, a warning is issued and baseline is skipped
"""
function baseline!(dat::EegData, baseline_selection::Function; channel_selection::Function = channels())
    # Apply predicate to get baseline mask
    baseline_mask = baseline_selection(dat.data)
    baseline_indices = findall(baseline_mask)

    # Skip baseline if window matches all samples (default) or no samples
    if isempty(baseline_indices)
        @minimal_warning "Baseline selection predicate matched no samples. Skipping baseline correction."
        return nothing
    end

    if length(baseline_indices) >= n_samples(dat)
        # Predicate selected all samples (default case) - skip baseline
        return nothing
    end

    # Convert to interval
    baseline_interval = IntervalIndex(start = first(baseline_indices), stop = last(baseline_indices))

    # Call the interval-based method
    baseline!(dat, baseline_interval; channel_selection = channel_selection)
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
- `baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real},Function}`: Time/index interval for baseline calculation, or a predicate function
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by subtracting the baseline mean
- Uses the specified time/index interval for baseline calculation
- If a predicate function is provided, it's applied to each epoch's data to determine the baseline window
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
    baseline!(dat::Vector{EpochData}, baseline_selection::Function; channel_selection=channels())

Apply baseline correction in-place to a vector of EpochData objects using a predicate function.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data to baseline correct
- `baseline_selection::Function`: Sample selection predicate (e.g., `samples((start, end))`)
- `channel_selection::Function`: Channel selection predicate (default: channels() - all channels)

# Notes
- Modifies each EpochData in the vector in-place by subtracting the baseline mean
- The predicate function is applied to each epoch's data DataFrame to determine the baseline window
"""
function baseline!(dat::Vector{EpochData}, baseline_selection::Function; channel_selection::Function = channels())
    baseline!.(dat, Ref(baseline_selection); channel_selection = channel_selection)
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
function _default_baseline_output_dir(
    input_dir::String,
    pattern::String,
    baseline_interval::Union{IntervalIndex,IntervalTime},
)
    interval_str = "$(baseline_interval.start)_to_$(baseline_interval.stop)"
    joinpath(input_dir, "baseline_$(pattern)_$(interval_str)")
end

function _default_baseline_output_dir(input_dir::String, pattern::String, baseline_selection::Function)
    # For predicates, use a generic name (could be improved to extract time range from predicate)
    joinpath(input_dir, "baseline_$(pattern)_predicate")
end


"""
Process a single file through baseline correction pipeline.
Returns BatchResult with success/failure info.
"""
function _process_baseline_file(
    filepath::String,
    output_path::String,
    baseline_interval::Union{IntervalIndex,IntervalTime,Tuple{Real,Real}},
    condition_selection::Function,
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
    data = _condition_select(data, condition_selection)

    # Apply baseline correction (mutates data in-place)
    baseline!.(data, Ref(baseline_interval))

    # Save (always use "data" as variable name since load_data finds by type)
    jldsave(output_path; data = data)

    return BatchResult(true, filename, "Baseline corrected successfully")
end

function _process_baseline_file(
    filepath::String,
    output_path::String,
    baseline_selection::Function,
    condition_selection::Function,
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
    data = _condition_select(data, condition_selection)

    # Apply baseline correction using predicate (mutates data in-place)
    baseline!.(data, Ref(baseline_selection))

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
             participant_selection::Function = participants(),
             condition_selection::Function = conditions(),
             output_dir::Union{String, Nothing} = nothing)

Apply baseline correction to EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `baseline_interval`: Baseline interval as:
  - `Tuple{Real,Real}`: Time interval in seconds (e.g., `(-0.2, 0.0)`)
  - `IntervalTime`: Time interval struct (e.g., `IntervalTime(start=-0.2, stop=0.0)`)
  - `IntervalIndex`: Index interval struct (e.g., `IntervalIndex(start=1, stop=51)`)
  - `Function`: Sample selection predicate (e.g., `samples((start, end))`)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `condition_selection::Function`: Condition selection predicate (default: `conditions()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on baseline settings)

# Example
```julia
# Baseline correct all epochs from -0.2 to 0.0 seconds
baseline("epochs", (-0.2, 0.0))

# Baseline correct using predicate
baseline("epochs", samples((-0.2, 0.0)))

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

        # Validate and normalize baseline interval (converts tuple to IntervalTime, validates structure)
        baseline_interval = validate_baseline_interval(baseline_interval)

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
        process_fn =
            (input_path, output_path) ->
                _process_baseline_file(input_path, output_path, baseline_interval, condition_selection)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Baseline correction")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

function baseline(
    file_pattern::String,
    baseline_selection::Function;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    condition_selection::Function = conditions(),
    output_dir::Union{String,Nothing} = nothing,
)
    # Convert predicate to interval by loading a sample file
    # (The interval-based version will handle all file finding, validation, etc.)
    files = _find_batch_files(file_pattern, input_dir, participant_selection)
    sample_data = load_data(joinpath(input_dir, files[1]))
    first_item = sample_data isa Vector{<:Union{ErpData,EpochData}} ? sample_data[1] : sample_data
    sample_df = first_item.data isa Vector{DataFrame} ? first_item.data[1] : first_item.data

    baseline_indices = findall(baseline_selection(sample_df))
    isempty(baseline_indices) && error("Baseline selection matched no samples")
    length(baseline_indices) >= nrow(sample_df) && return nothing

    baseline_interval = IntervalIndex(start = first(baseline_indices), stop = last(baseline_indices))
    baseline(
        file_pattern,
        baseline_interval;
        input_dir = input_dir,
        participant_selection = participant_selection,
        condition_selection = condition_selection,
        output_dir = output_dir,
    )
end

# ============================================================================
# TF Baseline for TimeFreqData structs
# ============================================================================

"""
    tf_baseline!(tf_data::TimeFreqData, baseline_window; method=:db)

Apply baseline correction to TimeFreqData in-place.

# Arguments
- `tf_data::TimeFreqData`: Time-frequency data to baseline correct
- `baseline_window::Tuple{Real,Real}`: Time window for baseline (start, stop) in seconds

# Keyword Arguments
- `method::Symbol=:db`: Baseline method
  - `:db`: Decibel change (10 * log10(power/baseline))
  - `:percent`: Percent change (100 * (power - baseline) / baseline)
  - `:relchange`: Relative change (power / baseline)

# Example
```julia
tf_baseline!(tf_data, (-0.3, 0.0); method=:db)
```
"""
function tf_baseline!(tf_data::TimeFreqData, baseline_window::Tuple{Real,Real}; method::Symbol=:db)
    # Check if baseline has already been applied
    if tf_data.baseline !== nothing
        error("Baseline correction has already been applied to this data (method: $(tf_data.baseline.method), window: $(tf_data.baseline.window)). " *
              "Baseline corrections are non-linear and cannot be chained. Use the original data to apply a different baseline.")
    end
    
    times = unique(tf_data.data_power.time)
    freqs_unique = unique(tf_data.data_power.freq)
    n_freqs = length(freqs_unique)
    n_times = length(times)
    
    # Find baseline time indices
    base_mask = baseline_window[1] .<= times .<= baseline_window[2]
    if !any(base_mask)
        error("Baseline window $(baseline_window) does not overlap with data times")
    end
    
    # Get channel columns
    ch_labels = channel_labels(tf_data)
    
    # Process each channel (baseline correction only applies to power, not phase)
    for ch in ch_labels
        # Reshape to freq × time matrix for baseline calculation
        power_mat = zeros(n_freqs, n_times)
        for ti in 1:n_times
            for fi in 1:n_freqs
                row_idx = (ti - 1) * n_freqs + fi
                power_mat[fi, ti] = tf_data.data_power[row_idx, ch]
            end
        end
        
        # Compute baseline per frequency
        baseline_power = mean(power_mat[:, base_mask], dims=2)
        
        # Apply baseline correction
        if method == :db
            # Avoid log(0) by ensuring positive values
            min_power = max.(baseline_power, 1e-10)
            power_mat .= 10 .* log10.(max.(power_mat, 1e-10) ./ min_power)
        elseif method == :percent
            # Avoid division by zero - use minimum threshold
            min_baseline = max.(baseline_power, 1e-10)
            power_mat .= 100 .* (power_mat .- baseline_power) ./ min_baseline
        elseif method == :relchange
            # Avoid division by zero
            min_baseline = max.(baseline_power, 1e-10)
            power_mat .= power_mat ./ min_baseline
        else
            error("Unknown baseline method: $method. Use :db, :percent, or :relchange")
        end
        
        # Write back to DataFrame
        for ti in 1:n_times
            for fi in 1:n_freqs
                row_idx = (ti - 1) * n_freqs + fi
                tf_data.data_power[row_idx, ch] = power_mat[fi, ti]
            end
        end
    end
    
    # Store baseline information
    tf_data.baseline = BaselineInfo(method=method, window=(Float64(baseline_window[1]), Float64(baseline_window[2])))
    
    return nothing
end

# Non-mutating version
function tf_baseline(tf_data::TimeFreqData, baseline_window::Tuple{Real,Real}; method::Symbol=:db)
    tf_copy = copy(tf_data)  # Use custom copy method instead of deepcopy
    tf_baseline!(tf_copy, baseline_window; method=method)
    return tf_copy
end

# Vector version
function tf_baseline!(tf_data::Vector{TimeFreqData}, baseline_window::Tuple{Real,Real}; method::Symbol=:db)
    tf_baseline!.(tf_data, Ref(baseline_window); method=method)
    return nothing
end

function tf_baseline(tf_data::Vector{TimeFreqData}, baseline_window::Tuple{Real,Real}; method::Symbol=:db)
    return [tf_baseline(tf, baseline_window; method=method) for tf in tf_data]
end

"""
    tf_baseline(file_pattern::String, baseline_window;
                method=:db, input_dir=pwd(), output_dir=nothing,
                participant_selection=participants(), condition_selection=conditions())

Apply baseline correction to TimeFreqData files.

# Example
```julia
tf_baseline("tf_epochs_wavelet", (-0.3, 0.0); method=:db)
```
"""
function tf_baseline(
    file_pattern::String,
    baseline_window::Tuple{Real,Real};
    method::Symbol=:db,
    input_dir::String=pwd(),
    output_dir::Union{String,Nothing}=nothing,
    participant_selection::Function=participants(),
    condition_selection::Function=conditions()
)
    log_file = "tf_baseline.log"
    setup_global_logging(log_file)

    try
        @info "Batch TF baseline started at $(now())"
        @log_call "tf_baseline"

        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_baseline_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files, baseline: $baseline_window, method: $method"

        process_fn = (input_path, output_path) -> begin
            filename = basename(input_path)
            data = load_data(input_path)
            if isnothing(data) || !(data isa Vector{TimeFreqData})
                return BatchResult(false, filename, "Invalid data type")
            end
            data = _condition_select(data, condition_selection)
            tf_baseline!(data, baseline_window; method=method)
            jldsave(output_path; data=data)
            return BatchResult(true, filename, "Baseline corrected")
        end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                       operation_name="TF baseline")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
