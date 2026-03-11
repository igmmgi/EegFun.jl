"""
Helper function to resample a single DataFrame by downsampling.
Handles trigger preservation and sample column updates.
"""
function _resample_dataframe!(df::DataFrame, factor::Int, trigger_col::Symbol)
    # Get indices of samples to keep (regular downsampling grid)
    keep_indices = 1:factor:nrow(df)

    # If trigger column exists, preserve triggers by scaling their positions
    if hasproperty(df, trigger_col)
        # Find all triggers in original data
        trigger_indices = findall(df[!, trigger_col] .!= 0)
        trigger_values = df[!, trigger_col][trigger_indices]

        # Downsample the data
        df_new = df[keep_indices, :]

        # Clear all triggers in downsampled data
        df_new[!, trigger_col] .= 0

        # For each original trigger, scale its position by the downsampling factor
        for (orig_idx, trig_val) in zip(trigger_indices, trigger_values)
            new_idx = round(Int, orig_idx / factor)
            new_idx = clamp(new_idx, 1, nrow(df_new))
            df_new[!, trigger_col][new_idx] = trig_val
        end

        df_resampled = df_new
    else # No triggers, just downsample
        df_resampled = df[keep_indices, :]
    end

    return df_resampled
end

"""
    resample!(dat::Union{ContinuousData, ErpData}, factor::Int)
    resample!(dat::EpochData, factor::Int)
    resample!(data_vec::Vector{T}, factor::Int) where {T<:EegData}
    resample(dat, factor)
    resample(data_vec::Vector{T}, factor) where {T<:EegData}
    resample(file_pattern::String, factor::Int; input_dir = pwd(), participant_selection = participants(), output_dir = nothing)

Downsample EEG data by keeping every `factor`-th sample.

Mutating (`!`) variants modify in-place and update `sample_rate`; non-mutating variants return a copy.
The `file_pattern` form batch-processes JLD2 files across participants.

⚠️ Simple decimation is used — **low-pass filter the data first** to avoid aliasing.

# Notes
- `factor` must be a positive integer; `sample_rate` must be divisible by `factor`
- All DataFrame columns (triggers, time, metadata) are preserved
- Trigger positions are scaled to the new sample grid for `ContinuousData`

# Examples
```julia
resample!(dat, 2)                    # 512 Hz → 256 Hz in-place
dat_256 = resample(dat, 2)           # non-mutating copy
resample!(all_epochs, 4)             # Vector: broadcasts to every element
resample("continuous", 2)            # batch: all JLD2 files in current dir
resample("continuous", 2, input_dir = "/data/study1", participant_selection = participants(1:20))
```
"""
function resample!(dat::SingleDataFrameEeg, factor::Int)::Nothing

    # Validation
    if factor < 1
        @minimal_error("Downsampling factor must be positive, got $factor")
    end

    if factor == 1
        @info "Downsampling factor is 1, no resampling needed"
        return nothing
    end

    if dat.sample_rate % factor != 0
        @minimal_error(
            "Sample rate $(dat.sample_rate) Hz is not evenly divisible by factor $factor. " *
            "New sample rate would be $(dat.sample_rate / factor) Hz. " *
            "Choose a factor that results in an integer sample rate."
        )
    end

    @info "Resampling data from $(dat.sample_rate) Hz to $(dat.sample_rate ÷ factor) Hz (factor: $factor)"

    # Resample the DataFrame
    dat.data = _resample_dataframe!(dat.data, factor, :trigger)

    # For continuous/ERP data, renumber sample column to be sequential
    if hasproperty(dat.data, :sample)
        dat.data.sample = 1:nrow(dat.data)
    end

    # Update sample rate
    dat.sample_rate = dat.sample_rate ÷ factor

    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data)) samples"

    return nothing
end




function resample!(dat::EpochData, factor::Int)::Nothing
    # Validation
    if factor < 1
        @minimal_error("Downsampling factor must be positive, got $factor")
    end

    if factor == 1
        @info "Downsampling factor is 1, no resampling needed"
        return nothing
    end

    if dat.sample_rate % factor != 0
        @minimal_error(
            "Sample rate $(dat.sample_rate) Hz is not evenly divisible by factor $factor. " *
            "New sample rate would be $(dat.sample_rate / factor) Hz. " *
            "Choose a factor that results in an integer sample rate."
        )
    end

    @info "Resampling $(length(dat.data)) epochs from $(dat.sample_rate) Hz to $(dat.sample_rate ÷ factor) Hz (factor: $factor)"

    # Downsample each epoch
    for (i, epoch) in enumerate(dat.data)
        dat.data[i] = _resample_dataframe!(epoch, factor, :trigger)

        # For epoch data, scale the first sample number and make rest sequential
        if hasproperty(dat.data[i], :sample)
            first_sample = round(Int, dat.data[i].sample[1] / factor)
            n_samples = nrow(dat.data[i])
            dat.data[i].sample = first_sample:(first_sample+n_samples-1)
        end
    end

    # Update sample rate
    dat.sample_rate = dat.sample_rate ÷ factor

    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"

    return nothing
end




function resample(dat::T, factor::Int)::T where {T<:EegData}
    dat_copy = copy(dat)
    resample!(dat_copy, factor)
    return dat_copy
end


function resample(data_vec::Vector{T}, factor::Int)::Vector{T} where {T<:EegData}
    return [resample(dat, factor) for dat in data_vec]
end


function resample!(data_vec::Vector{T}, factor::Int)::Nothing where {T<:EegData}
    for dat in data_vec
        resample!(dat, factor)
    end
    return nothing
end



# =============================================================================
#     BATCH PROCESSING FUNCTIONS
# =============================================================================

"""Generate default output directory name for resampling operation."""
function _default_resample_output_dir(input_dir::String, pattern::String, factor::Int)
    new_rate_str = "resampled_by_$(factor)"
    joinpath(input_dir, "$(new_rate_str)_$(pattern)")
end


"""
Process a single data file through resampling pipeline.
Returns BatchResult with success/failure info.
"""
function _process_resample_file(filepath::String, output_path::String, factor::Int)
    filename = basename(filepath)

    # Read data using read_data (handles single variable files automatically)
    loaded_data = read_data(filepath)

    if isnothing(loaded_data)
        return BatchResult(false, filename, "No data found in file")
    end

    if !(loaded_data isa Union{ContinuousData,EpochData,ErpData})
        return BatchResult(false, filename, "Data is not a recognized EEG data type")
    end

    # Resample
    try
        old_rate = loaded_data.sample_rate
        resampled_data = resample(loaded_data, factor)
        new_rate = resampled_data.sample_rate

        # Save results (always use "data" as variable name since read_data finds by type)
        jldsave(output_path; data = resampled_data)

        message = "Resampled from $old_rate Hz to $new_rate Hz (factor: $factor)"
        return BatchResult(true, filename, message)
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end



function resample(
    file_pattern::String,
    factor::Int;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "resample.log"
    setup_global_logging(log_file)

    try
        @info "Batch resampling started at $(now())"
        @log_call "resample"

        # Validation
        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        if factor <= 0
            @minimal_error("Downsampling factor must be positive, got $factor")
        end

        # Setup directories
        output_dir = something(output_dir, _default_resample_output_dir(input_dir, file_pattern, factor))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Downsampling factor: $factor"

        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> _process_resample_file(input_path, output_path, factor)

        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Resampling")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
