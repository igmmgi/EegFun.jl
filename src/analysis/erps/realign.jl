"""
Realignment of epoched EEG data to different time points.

Realign stimulus-locked epoched data to a different
event time point (e.g., response time, saccade onset, etc.). This is useful for
creating response-locked waveforms from stimulus-locked epochs.
"""

# === CALCULATING TRIGGER INTERVALS ===

"""
    calculate_trigger_interval!(dat::EpochData, start_triggers::Vector{Int}, end_triggers::Vector{Int}; column_name::Symbol = :interval)
    calculate_trigger_interval!(dat::Vector{EpochData}, start_triggers::Vector{Int}, end_triggers::Vector{Int}; column_name::Symbol = :interval)

Calculate the time interval between two triggers within each epoch. Specifically, it searches for 
the first occurrence of any trigger in `start_triggers`, and the first occurrence of any trigger 
in `end_triggers`. The time difference (end - start) in seconds is appended as a constant column 
to each epoch's DataFrame.

If either the start trigger or end trigger is missing in an epoch, the interval value is set to `NaN`.

# Arguments
- `dat`: `EpochData` or `Vector{EpochData}`
- `start_triggers::Vector{Int}`: A list of trigger codes for the start interval.
- `end_triggers::Vector{Int}`: A list of trigger codes for the end interval.
- `column_name::Symbol`: The name of the column to append (default: `:interval`).

# Examples
```julia
calculate_trigger_interval!(epochs, [101, 103], [111, 112])
```
"""
function calculate_trigger_interval!(
    dat::EpochData,
    start_triggers::Vector{Int},
    end_triggers::Vector{Int};
    column_name::Symbol = :interval,
)::Nothing
    @info "Calculating trigger intervals for EpochData (condition: $(dat.condition_name))"

    n_missing = 0

    for epoch in dat.data
        start_idx = findfirst(x -> x in start_triggers, epoch.trigger)
        end_idx = findfirst(x -> x in end_triggers, epoch.trigger)

        if isnothing(start_idx) || isnothing(end_idx)
            epoch[!, column_name] .= NaN
            n_missing += 1
        else
            interval_seconds = epoch.time[end_idx] - epoch.time[start_idx]
            epoch[!, column_name] .= interval_seconds
        end
    end

    if n_missing > 0
        @minimal_warning "Missing start or end trigger pair in $n_missing / $(length(dat.data)) epochs (set to NaN)."
    end

    return nothing
end

"""
    calculate_trigger_interval!(dat_vec::Vector{EpochData}, start_triggers, end_triggers; column_name = :interval)

Apply `calculate_trigger_interval!` to every condition in a `Vector{EpochData}` in-place.
"""
function calculate_trigger_interval!(
    dat_vec::Vector{EpochData},
    start_triggers::Vector{Int},
    end_triggers::Vector{Int};
    column_name::Symbol = :interval,
)::Nothing
    for dat in dat_vec
        calculate_trigger_interval!(dat, start_triggers, end_triggers; column_name = column_name)
    end
    return nothing
end

"""
    calculate_trigger_interval(dat, start_triggers, end_triggers; column_name = :interval)

Non-mutating version of `calculate_trigger_interval!`.
"""
function calculate_trigger_interval(
    dat::EpochData,
    start_triggers::Vector{Int},
    end_triggers::Vector{Int};
    column_name::Symbol = :interval,
)::EpochData
    dat_copy = EpochData(
        dat.file,
        dat.condition,
        dat.condition_name,
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )
    calculate_trigger_interval!(dat_copy, start_triggers, end_triggers; column_name = column_name)
    return dat_copy
end

"""
    calculate_trigger_interval(dat_vec::Vector{EpochData}, start_triggers, end_triggers; column_name = :interval)

Non-mutating `Vector{EpochData}` overload. Returns a new `Vector{EpochData}` with the interval column added.
"""
function calculate_trigger_interval(
    dat_vec::Vector{EpochData},
    start_triggers::Vector{Int},
    end_triggers::Vector{Int};
    column_name::Symbol = :interval,
)::Vector{EpochData}
    return [calculate_trigger_interval(dat, start_triggers, end_triggers; column_name = column_name) for dat in dat_vec]
end


# === CORE REALIGNMENT FUNCTIONS ===


"""
    realign!(dat::EpochData, realignment_triggers::Vector{Int})::Nothing

Realign epoched data in-place to an event specified by trigger codes.

Realign stimulus-locked epoched data so
that the time of the first matching `realignment_trigger` becomes the new time zero. 
After realignment, all epochs are cropped to a common time interval determined by 
the latest start time and earliest end time across all valid epochs. 

If an epoch does NOT contain any of the `realignment_triggers`, it is automatically 
dropped from the dataset to prevent silent mismatches.

# Arguments
- `dat::EpochData`: Epoched EEG data to realign
- `realignment_triggers::Vector{Int}`: A list of trigger codes to align to.

# Effects
- Modifies the input data in-place (drops invalid epochs)
- Updates the `:time` column in each epoch
- Crops all epochs to a common time interval

# Examples
```julia
using EegFun

# Load stimulus-locked epoched data
epochs = read_data("participant_1_epochs.jld2")

# Realign to response trigger codes
realign!(epochs, [127, 191, 223, 239])
```

# Notes
- After realignment, epoch length may be shorter due to cropping to common interval.
"""
function realign!(dat::EpochData, realignment_triggers::Vector{Int})::Nothing

    @info "Realigning epoched data to triggers: $realignment_triggers"

    # Step 1: Filter and realign each epoch individually
    @info "Step 1/3: Realigning individual epochs"
    _realign_epochs!(dat, realignment_triggers)

    if isempty(dat.data)
        @minimal_warning "All epochs were dropped due to missing realignment triggers!"
        return nothing
    end

    # Step 2: Find common time interval across all realigned epochs
    @info "Step 2/3: Finding common time interval"
    common_interval = _find_common_time_window(dat)

    # Calculate n_samples to strictly avoid off-by-one errors within this single condition
    n_samples = minimum(dat.data) do epoch
        s_idx = find_closest_time_index(epoch.time, common_interval[1])
        e_idx = find_closest_time_index(epoch.time, common_interval[2])
        e_idx - s_idx + 1
    end

    @info "Common time interval: $(common_interval[1]) s to $(common_interval[2]) s (strict: $n_samples samples)"

    # Step 3: Crop all epochs to common interval
    @info "Step 3/3: Cropping epochs to common interval"
    _crop_epochs_to_window!(dat, common_interval, n_samples)

    @info "Realignment complete. Final epoch length: $(nrow(dat.data[1])) samples"

    return nothing
end


"""
    realign(dat::EpochData, realignment_triggers::Vector{Int})::EpochData

Non-mutating version of realign!. Returns a new EpochData object with realigned epochs.
"""
function realign(dat::EpochData, realignment_triggers::Vector{Int})::EpochData
    # Create a deep copy of the data (preserve struct fields)
    dat_copy = EpochData(
        dat.file,
        dat.condition,
        dat.condition_name,
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )

    realign!(dat_copy, realignment_triggers)

    return dat_copy
end

"""
    realign!(dat_vec::Vector{EpochData}, realignment_triggers; global_interval = nothing, global_n_samples = nothing)

In-place realignment for a `Vector{EpochData}`. All conditions are realigned to a **shared
global time window** so that every condition ends up with an identical sample count — a
strict requirement for LRP averaging across conditions.
"""
function realign!(
    dat_vec::Vector{EpochData},
    realignment_triggers::Vector{Int};
    global_interval::Union{Tuple{Float64,Float64},Nothing} = nothing,
    global_n_samples::Union{Int,Nothing} = nothing,
)::Nothing
    @info "Realigning $(length(dat_vec)) conditions in Vector{EpochData}..."

    # Step 1: Validate and realign all epochs individually
    for dat in dat_vec
        if !isempty(dat.data)
            _realign_epochs!(dat, realignment_triggers)
        end
    end

    # Step 2: Find global common time interval across all conditions if not provided
    if isnothing(global_interval)
        common_starts = Float64[]
        common_ends = Float64[]
        for dat in dat_vec
            if !isempty(dat.data)
                start_time, end_time = _find_common_time_window(dat)
                push!(common_starts, start_time)
                push!(common_ends, end_time)
            end
        end

        if isempty(common_starts)
            @info "No data to realign."
            return nothing
        end

        global_interval = (maximum(common_starts), minimum(common_ends))

        # Pre-calculate global strict n_samples
        all_n_samples = Int[]
        for dat in dat_vec
            for epoch in dat.data
                s_idx = find_closest_time_index(epoch.time, global_interval[1])
                e_idx = find_closest_time_index(epoch.time, global_interval[2])
                push!(all_n_samples, e_idx - s_idx + 1)
            end
        end
        global_n_samples = minimum(all_n_samples)
    end

    @info "Global common time interval: $(global_interval[1]) s to $(global_interval[2]) s (strict: $global_n_samples samples)"

    # Step 3: Crop all epochs in all conditions to this global interval and global length
    for dat in dat_vec
        if !isempty(dat.data)
            _crop_epochs_to_window!(dat, global_interval, global_n_samples)
        end
    end

    # Just grab an epoch length to print (assume first cond with data)
    for dat in dat_vec
        if !isempty(dat.data)
            @info "Realignment complete. Final epoch length: $(nrow(dat.data[1])) samples"
            break
        end
    end

    return nothing
end

"""
    realign(dat_vec::Vector{EpochData}, realignment_triggers; global_interval = nothing, global_n_samples = nothing)

Non-mutating `Vector{EpochData}` overload of `realign!`. Returns a deep-copied, globally realigned vector.
"""
function realign(
    dat_vec::Vector{EpochData},
    realignment_triggers::Vector{Int};
    global_interval::Union{Tuple{Float64,Float64},Nothing} = nothing,
    global_n_samples::Union{Int,Nothing} = nothing,
)::Vector{EpochData}
    dat_vec_copy = [
        EpochData(
            dat.file,
            dat.condition,
            dat.condition_name,
            [copy(epoch, copycols = true) for epoch in dat.data],
            copy(dat.layout),
            dat.sample_rate,
            copy(dat.analysis_info),
        ) for dat in dat_vec
    ]
    realign!(dat_vec_copy, realignment_triggers; global_interval = global_interval, global_n_samples = global_n_samples)
    return dat_vec_copy
end



# === INTERNAL HELPER FUNCTIONS ===

"""
Realign each epoch by searching for the first occurrence of the realignment triggers 
and shifting the time vector so that this event becomes time zero.
Episodes missing the trigger are removed from `dat.data`.
"""
function _realign_epochs!(dat::EpochData, realignment_triggers::Vector{Int})
    # Filter out epochs missing the trigger and realign in a single pass
    valid_data = DataFrame[]
    for (i, epoch) in enumerate(dat.data)
        resp_idx = findfirst(x -> x in realignment_triggers, epoch.trigger)
        if !isnothing(resp_idx)
            # Find trigger index and shift time vector so realignment_time becomes 0
            realignment_time = epoch.time[resp_idx]
            epoch[!, :time] .-= realignment_time
            push!(valid_data, epoch)
        else
            @minimal_warning "Epoch $i dropped: missing realignment trigger."
        end
    end
    dat.data = valid_data
end


"""
Find the common time interval that is valid for all realigned epochs.

This finds the latest start time (maximum of all minimum times) and the earliest
end time (minimum of all maximum times) across all epochs.
"""
function _find_common_time_window(dat::EpochData)::Tuple{Float64,Float64}
    # Common interval is the overlap of all individual intervals
    # Since time is sorted, [1] is minimum and [end] is maximum (O(1) vs O(T))

    # Latest start time
    common_start = maximum(epoch.time[1] for epoch in dat.data)
    # Earliest end time
    common_end = minimum(epoch.time[end] for epoch in dat.data)

    if common_start >= common_end
        @minimal_error(
            "No common time interval found after realignment. " *
            "Common start: $(common_start) s, Common end: $(common_end) s. " *
            "The original epochs may not be long enough to accommodate all realignment times."
        )
    end

    return (common_start, common_end)
end


"""
Crop all epochs to the specified time interval.

Uses sample-count based cropping to ensure all epochs have exactly the same length,
avoiding floating-point precision issues that can result in different epoch lengths.
After cropping, regenerates uniform time vectors for all epochs.
"""
function _crop_epochs_to_window!(dat::EpochData, window::Tuple{Float64,Float64}, global_n_samples::Union{Int,Nothing} = nothing)
    start_time, end_time = window

    # First pass: find start and end indices for each epoch
    indices = []
    for epoch in dat.data
        start_idx = find_closest_time_index(epoch.time, start_time)
        end_idx = find_closest_time_index(epoch.time, end_time)
        push!(indices, (start_idx, end_idx))
    end

    # Find the minimum valid range across all epochs
    # (some epochs might have fewer samples due to varying RTs)
    n_samples = minimum([end_idx - start_idx + 1 for (start_idx, end_idx) in indices])

    # Enforce global length mathematically if requested
    if !isnothing(global_n_samples)
        n_samples = global_n_samples
    end

    # Second pass: crop all epochs to this exact length
    for (i, epoch) in enumerate(dat.data)
        start_idx, _ = indices[i]
        end_idx = start_idx + n_samples - 1

        # Validate that we have enough samples
        if end_idx > nrow(epoch)
            @minimal_error(
                "Epoch $i does not have enough samples for common interval. " *
                "Need samples up to index $end_idx, but epoch only has $(nrow(epoch)) samples. " *
                "Window: [$start_time, $end_time]"
            )
        end

        # Crop to exact sample range
        dat.data[i] = epoch[start_idx:end_idx, :]
    end

    # Third pass: Regenerate uniform time vector for all epochs
    # This ensures all epochs have identical time vectors
    uniform_time = range(start_time, stop = end_time, length = n_samples)

    for epoch in dat.data
        epoch[!, :time] = collect(uniform_time)
    end
end


# === BATCH PROCESSING FUNCTIONS ===

"""Generate default output directory name for realignment operation."""
function _default_realign_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "realigned_$(clean_pattern(pattern))")
end


"""
Process a single epoch file through realignment pipeline.
Returns BatchResult with success/failure info.
"""
function _process_realign_file(
    filepath::String,
    output_path::String,
    realignment_triggers::Vector{Int},
    global_interval::Tuple{Float64,Float64},
    global_n_samples::Int,
)
    filename = basename(filepath)

    # Read data using read_data (handles single variable files automatically)
    epochs_data = read_data(filepath)

    if isnothing(epochs_data)
        return BatchResult(false, filename, "No data found in file")
    end

    if !(epochs_data isa EpochData || epochs_data isa Vector{<:EpochData})
        return BatchResult(false, filename, "Data is not EpochData type")
    end

    # Realign
    try
        if epochs_data isa EpochData
            # For a single EpochData object, we cannot enforce multi-condition constraints directly through the simple API
            # But in this batch mode we are passing global interval/n_samples.
            # So we manually apply them here to ensure batch adherence:
            _realign_epochs!(epochs_data, realignment_triggers)
            if isempty(epochs_data.data)
                return BatchResult(true, filename, "All epochs dropped (missing trigger).")
            end
            _crop_epochs_to_window!(epochs_data, global_interval, global_n_samples)
            realigned_epochs = epochs_data
        else
            realigned_epochs =
                realign(epochs_data, realignment_triggers; global_interval = global_interval, global_n_samples = global_n_samples)
        end

        # Save results
        jldsave(output_path; data = realigned_epochs)

        if realigned_epochs isa Vector
            n_epochs = sum(length(cond.data) for cond in realigned_epochs)
            n_samples = isempty(realigned_epochs) || isempty(realigned_epochs[1].data) ? 0 : nrow(realigned_epochs[1].data[1])
        else
            n_epochs = length(realigned_epochs.data)
            n_samples = isempty(realigned_epochs.data) ? 0 : nrow(realigned_epochs.data[1])
        end
        return BatchResult(true, filename, "Realigned $n_epochs total epochs to $n_samples samples each")
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end

function realign(
    file_pattern::String,
    realignment_triggers::Vector{Int};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "realign.log"
    setup_global_logging(log_file)

    try
        @info "Batch realignment started at $(now())"
        @log_call "realign"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_realign_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"

        # --- PASS 1: SCAN GLOBAL BOUNDS ---
        @info "Pass 1/3: Global timeline scanning across ($(length(files))) files..."
        grand_starts = Float64[]
        grand_ends = Float64[]

        for file in files
            input_path = joinpath(input_dir, file)
            epochs_data = read_data(input_path)
            if !isnothing(epochs_data)
                vec_data = epochs_data isa Vector ? epochs_data : [epochs_data]
                for dat in vec_data
                    if !isempty(dat.data)
                        _realign_epochs!(dat, realignment_triggers)
                        if !isempty(dat.data)
                            s, e = _find_common_time_window(dat)
                            push!(grand_starts, s)
                            push!(grand_ends, e)
                        end
                    end
                end
            end
        end

        if isempty(grand_starts)
            @minimal_error("No valid realigned data found across any file! (All trials likely missing the trigger)")
        end

        global_interval = (maximum(grand_starts), minimum(grand_ends))

        # --- PASS 2: EXACT SAMPLE CROP SCAN ---
        @info "Pass 2/3: Enforcing exact sample constraints..."
        grand_n_samples = Int[]
        for file in files
            input_path = joinpath(input_dir, file)
            epochs_data = read_data(input_path)
            if !isnothing(epochs_data)
                vec_data = epochs_data isa Vector ? epochs_data : [epochs_data]
                for dat in vec_data
                    _realign_epochs!(dat, realignment_triggers)
                    for epoch in dat.data
                        s_idx = find_closest_time_index(epoch.time, global_interval[1])
                        e_idx = find_closest_time_index(epoch.time, global_interval[2])
                        push!(grand_n_samples, e_idx - s_idx + 1)
                    end
                end
            end
        end
        global_n_samples = minimum(grand_n_samples)

        @info "Grand Global Interval: $(global_interval[1]) to $(global_interval[2]) (Strict: $global_n_samples samples)"

        # --- PASS 3: REALIGN & RECORD ---
        @info "Pass 3/3: Realigning to disk..."

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) ->
                _process_realign_file(input_path, output_path, realignment_triggers, global_interval, global_n_samples)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Realigning epochs")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

# Alias for string pattern since users often write realign!("pattern", ...) intuitively wanting to batch process
function realign!(file_pattern::String, realignment_triggers::Vector{Int}; kwargs...)
    realign(file_pattern, realignment_triggers; kwargs...)
end


# === TRIGGER INTERVAL BATCH PROCESSING FUNCTIONS ===

"""Generate default output directory name for calculate_trigger_interval operation."""
function _default_calculate_trigger_interval_output_dir(input_dir::String, pattern::String, column_name::Symbol)
    joinpath(input_dir, "epochs_calculate_trigger_interval")
end

"""
Process a single epoch file through calculate_trigger_interval pipeline.
Returns BatchResult with success/failure info.
"""
function _process_calculate_trigger_interval_file(
    filepath::String,
    output_path::String,
    start_triggers::Vector{Int},
    end_triggers::Vector{Int},
    column_name::Symbol,
)
    filename = basename(filepath)

    # Read data
    epochs_data = read_data(filepath)

    if isnothing(epochs_data)
        return BatchResult(false, filename, "No data found in file")
    end

    if !(epochs_data isa EpochData || epochs_data isa Vector{<:EpochData})
        return BatchResult(false, filename, "Data is not EpochData type")
    end

    # Calculate Interval
    try
        epochs_with_interval = calculate_trigger_interval(epochs_data, start_triggers, end_triggers; column_name = column_name)

        # Save results
        jldsave(output_path; data = epochs_with_interval)

        if epochs_with_interval isa Vector
            n_epochs = sum(length(cond.data) for cond in epochs_with_interval)
        else
            n_epochs = length(epochs_with_interval.data)
        end
        return BatchResult(true, filename, "Appended interval column (:$column_name) to $n_epochs total epochs")
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end


"""
    calculate_trigger_interval!(file_pattern::String, start_triggers::Vector{Int}, end_triggers::Vector{Int}; kwargs...)
    calculate_trigger_interval(file_pattern::String, start_triggers::Vector{Int}, end_triggers::Vector{Int}; kwargs...)

Batch process JLD2 files that match `file_pattern`, automatically computing time intervals 
between `start_triggers` and `end_triggers` and appending an interval column to the epochs.
"""
function calculate_trigger_interval(
    file_pattern::String,
    start_triggers::Vector{Int},
    end_triggers::Vector{Int};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
    column_name::Symbol = :interval,
)
    # Setup logging
    log_file = "calculate_trigger_interval.log"
    setup_global_logging(log_file)

    try
        @info "Batch calculate_trigger_interval started at $(now())"
        @log_call "calculate_trigger_interval"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_calculate_trigger_interval_output_dir(input_dir, file_pattern, column_name))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Looking for start triggers $(start_triggers) and end triggers $(end_triggers)"

        # Create processing function with captured parameters
        process_fn =
            (input_path, output_path) ->
                _process_calculate_trigger_interval_file(input_path, output_path, start_triggers, end_triggers, column_name)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Calculating trigger intervals")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end

# Alias for string pattern
function calculate_trigger_interval!(file_pattern::String, start_triggers::Vector{Int}, end_triggers::Vector{Int}; kwargs...)
    calculate_trigger_interval(file_pattern, start_triggers, end_triggers; kwargs...)
end
