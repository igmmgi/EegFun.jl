"""
Realignment of epoched EEG data to different time points.

This function takes stimulus-locked epoched data and realigns it to a different
event time point (e.g., response time, saccade onset, etc.). This is useful for
creating response-locked waveforms from stimulus-locked epochs.
"""

#=============================================================================
    CORE REALIGNMENT FUNCTIONS
=============================================================================#

"""
    realign!(dat::EpochData, realignment_column::Symbol)::Nothing

Realign epoched data in-place to a different time point specified in a DataFrame column.

This function takes stimulus-locked (or any reference-locked) epoched data and
realigns each epoch so that the time specified in `realignment_column` becomes
the new time zero. After realignment, all epochs are cropped to a common time
window determined by the latest start time and earliest end time across all epochs.

# Arguments
- `dat::EpochData`: Epoched EEG data to realign
- `realignment_column::Symbol`: Name of the column containing realignment times (e.g., `:rt`, `:response_time`)

# Effects
- Modifies the input data in-place
- Updates the `:time` column in each epoch
- Crops all epochs to a common time window
- The realignment column values should be relative to the current time zero

# Examples
```julia
using eegfun

# Load stimulus-locked epoched data
epochs = load("participant_1_epochs.jld2", "epochs")

# Realign to response times (stored in :rt column)
realign!(epochs, :rt)

# Now time=0 corresponds to the response for each trial
# and all epochs have the same time window
```

# Notes
- The realignment column must exist in all epochs
- The realignment value should be constant within each epoch (checked automatically)
- After realignment, epoch length may be shorter due to cropping to common window
- Trials with insufficient data before/after the realignment point will determine the final epoch length

# Use Cases
- Response-locked LRP analysis from stimulus-locked data
- Saccade-locked analysis from fixation-locked data
- Any re-referencing to a different event within the epoch
"""
function realign!(dat::EpochData, realignment_column::Symbol)::Nothing

    @info "Realigning epoched data to column: $realignment_column"

    # Validate that the realignment column exists
    _validate_realignment_column(dat, realignment_column)

    # Step 1: Realign each epoch individually
    @info "Step 1/3: Realigning individual epochs"
    _realign_epochs!(dat, realignment_column)

    # Step 2: Find common time window across all realigned epochs
    @info "Step 2/3: Finding common time window"
    common_window = _find_common_time_window(dat)

    @info "Common time window: $(common_window[1]) s to $(common_window[2]) s"

    # Step 3: Crop all epochs to common window
    @info "Step 3/3: Cropping epochs to common window"
    _crop_epochs_to_window!(dat, common_window)

    @info "Realignment complete. Final epoch length: $(nrow(dat.data[1])) samples"

    return nothing
end


"""
    realign(dat::EpochData, realignment_column::Symbol)::EpochData

Non-mutating version of realign!. Returns a new EpochData object with realigned epochs.
"""
function realign(dat::EpochData, realignment_column::Symbol)::EpochData
    # Create a deep copy of the data
    dat_copy = EpochData(
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )

    # Apply realignment
    realign!(dat_copy, realignment_column)

    return dat_copy
end


#=============================================================================
    INTERNAL HELPER FUNCTIONS
=============================================================================#

"""
Validate that the realignment column exists and has appropriate values.
"""
function _validate_realignment_column(dat::EpochData, realignment_column::Symbol)
    # Check that column exists in all epochs
    for (i, epoch) in enumerate(dat.data)
        if !hasproperty(epoch, realignment_column)
            @minimal_error_throw(
                "Realignment column :$realignment_column not found in epoch $i. " *
                "Available columns: $(propertynames(epoch))"
            )
        end

        # Check that the realignment value is constant within the epoch
        realignment_values = epoch[!, realignment_column]
        if !all(realignment_values .≈ realignment_values[1])
            @minimal_error_throw(
                "Realignment column :$realignment_column has varying values within epoch $i. " *
                "The realignment time should be constant for each epoch."
            )
        end

        # Check that the realignment value is finite
        if !isfinite(realignment_values[1])
            @minimal_error_throw(
                "Realignment column :$realignment_column has non-finite value (NaN or Inf) in epoch $i"
            )
        end
    end

    # Check that :time column exists
    if !hasproperty(dat.data[1], :time)
        @minimal_error_throw("Time column :time not found in epochs")
    end
end


"""
Realign each epoch by subtracting the realignment time from the time vector.
"""
function _realign_epochs!(dat::EpochData, realignment_column::Symbol)
    for (i, epoch) in enumerate(dat.data)
        # Get realignment time (should be constant within epoch)
        realignment_time = epoch[1, realignment_column]

        # Shift time vector so realignment_time becomes 0
        epoch[!, :time] .-= realignment_time

        # Set realignment column to exactly 0 (avoid floating-point residuals)
        epoch[!, realignment_column] .= 0.0
    end
end


"""
Find the common time window that is valid for all realigned epochs.

This finds the latest start time (maximum of all minimum times) and the earliest
end time (minimum of all maximum times) across all epochs.
"""
function _find_common_time_window(dat::EpochData)::Tuple{Float64,Float64}
    # Get time ranges for all epochs
    min_times = [minimum(epoch.time) for epoch in dat.data]
    max_times = [maximum(epoch.time) for epoch in dat.data]

    # Common window is the overlap of all individual windows
    # Latest start time
    common_start = maximum(min_times)
    # Earliest end time
    common_end = minimum(max_times)

    if common_start >= common_end
        @minimal_error_throw(
            "No common time window found after realignment. " *
            "Common start: $(common_start) s, Common end: $(common_end) s. " *
            "The original epochs may not be long enough to accommodate all realignment times."
        )
    end

    return (common_start, common_end)
end


"""
Crop all epochs to the specified time window.

Uses sample-count based cropping to ensure all epochs have exactly the same length,
avoiding floating-point precision issues that can result in different epoch lengths.
After cropping, regenerates uniform time vectors for all epochs.
"""
function _crop_epochs_to_window!(dat::EpochData, window::Tuple{Float64,Float64})
    start_time, end_time = window

    # First pass: find start and end indices for each epoch
    indices = []
    for epoch in dat.data
        start_idx = argmin(abs.(epoch.time .- start_time))
        end_idx = argmin(abs.(epoch.time .- end_time))
        push!(indices, (start_idx, end_idx))
    end

    # Find the minimum valid range across all epochs
    # (some epochs might have fewer samples due to varying RTs)
    n_samples = minimum([end_idx - start_idx + 1 for (start_idx, end_idx) in indices])

    # Second pass: crop all epochs to this exact length
    for (i, epoch) in enumerate(dat.data)
        start_idx, _ = indices[i]
        end_idx = start_idx + n_samples - 1

        # Validate that we have enough samples
        if end_idx > nrow(epoch)
            @minimal_error_throw(
                "Epoch $i does not have enough samples for common window. " *
                "Need samples up to index $end_idx, but epoch only has $(nrow(epoch)) samples. " *
                "Window: [$start_time, $end_time]"
            )
        end

        # Crop to exact sample range
        dat.data[i] = epoch[start_idx:end_idx, :]
    end

    # Third pass: Regenerate uniform time vector for all epochs
    # This ensures all epochs have identical time vectors
    dt = 1.0 / dat.sample_rate
    uniform_time = range(start_time, stop = end_time, length = n_samples)

    for epoch in dat.data
        epoch[!, :time] = collect(uniform_time)
    end
end


#=============================================================================
    BATCH PROCESSING FUNCTIONS
=============================================================================#

"""Generate default output directory name for realignment operation."""
function _default_realign_output_dir(input_dir::String, pattern::String, realignment_column::Symbol)
    joinpath(input_dir, "realigned_$(pattern)_$(realignment_column)")
end


"""
Process a single epoch file through realignment pipeline.
Returns BatchResult with success/failure info.
"""
function _process_realign_file(filepath::String, output_path::String, realignment_column::Symbol)
    filename = basename(filepath)

    # Load data
    file_data = load(filepath)

    # Try common variable names for epoched data
    epoch_var_names = ["epochs", "epoch_data", "data"]
    epochs_data = nothing

    for var_name in epoch_var_names
        if haskey(file_data, var_name)
            epochs_data = file_data[var_name]
            break
        end
    end

    if isnothing(epochs_data)
        return BatchResult(false, filename, "No epoched data variable found (tried: $(epoch_var_names))")
    end

    if !(epochs_data isa EpochData)
        return BatchResult(false, filename, "Data is not EpochData type")
    end

    # Realign
    try
        realigned_epochs = realign(epochs_data, realignment_column)

        # Save results
        save(output_path, "epochs", realigned_epochs)

        n_epochs = length(realigned_epochs.data)
        n_samples = nrow(realigned_epochs.data[1])
        return BatchResult(true, filename, "Realigned $n_epochs epochs to $n_samples samples each")
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end


"""
    realign(file_pattern::String, realignment_column::Symbol;
            input_dir::String = pwd(),
            participants::Union{Int, Vector{Int}, Nothing} = nothing,
            output_dir::Union{String, Nothing} = nothing)

Batch realign epoched data from JLD2 files and save to a new directory.

This function processes multiple participant files at once, realigning epoched
data to a different time point specified by a column in the epoch DataFrames.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "epochs", "epochs_cleaned")
- `realignment_column::Symbol`: Column name containing realignment times (e.g., `:rt`, `:response_time`)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)

# Examples
```julia
# Realign all participants' epochs to response times
realign("epochs_cleaned", :rt)

# Specific participants only
realign("epochs_cleaned", :rt, participants = [1, 2, 3])

# Custom output directory
realign("epochs_cleaned", :rt, output_dir = "/path/to/output")

# Full example workflow
realign("epochs_cleaned", :rt,
        input_dir = "/data/study1",
        participants = 1:20)
```

# Output
- Creates new directory with realigned epoch data files
- Each output file contains "epochs" variable with realigned EpochData
- Log file saved to output directory

# Notes
- Input files should contain EpochData with the specified realignment column
- All trials within each file are realigned and cropped to a common window
- Useful for response-locked LRP analysis or other event-locked analyses
"""
function realign(
    file_pattern::String,
    realignment_column::Symbol;
    input_dir::String = pwd(),
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "realign.log"
    setup_global_logging(log_file)

    try
        @info "Batch realignment started at $(now())"
        @log_call "realign" (file_pattern, realignment_column)

        # Validation
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir = something(output_dir, _default_realign_output_dir(input_dir, file_pattern, realignment_column))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Realigning to column: :$realignment_column"

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_realign_file(input_path, output_path, realignment_column)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Realigning epochs")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
