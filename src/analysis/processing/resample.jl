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
    resample!(dat::Union{ContinuousData, ErpData}, factor::Int)::Nothing

Resample data in-place by downsampling by the specified factor.

This function reduces the sampling rate of the data by keeping every nth sample,
where n is the downsampling factor. All metadata columns (triggers, trial info,
conditions, etc.) are preserved. The sample rate is updated accordingly.

# Arguments
- `dat::Union{ContinuousData, ErpData}`: Data to resample
- `factor::Int`: Downsampling factor (e.g., 2 = keep every 2nd sample, halving the sample rate)

# Effects
- Modifies the input data in-place
- Updates the sample rate: new_rate = old_rate / factor
- Preserves all metadata columns

# Examples
```julia
using EegFun

# Load continuous data at 512 Hz
data = load("participant_1_continuous.jld2", "data")

# Downsample to 256 Hz (factor of 2)
resample!(data, 2)
```

# Notes
- Factor must be a positive integer
- New sample rate must be an integer (old_rate must be divisible by factor)
- Simple decimation is used 
- For proper downsampling, low-pass filter data first to avoid aliasing
- All DataFrame columns are preserved, including time, triggers, and metadata
"""
function resample!(dat::SingleDataFrameEeg, factor::Int)::Nothing

    # Validation
    if factor < 1
        @minimal_error_throw("Downsampling factor must be positive, got $factor")
    end

    if factor == 1
        @info "Downsampling factor is 1, no resampling needed"
        return nothing
    end

    if dat.sample_rate % factor != 0
        @minimal_error_throw(
            "Sample rate $(dat.sample_rate) Hz is not evenly divisible by factor $factor. " *
            "New sample rate would be $(dat.sample_rate / factor) Hz. " *
            "Choose a factor that results in an integer sample rate."
        )
    end

    @info "Resampling data from $(dat.sample_rate) Hz to $(dat.sample_rate ÷ factor) Hz (factor: $factor)"

    # Resample the DataFrame
    dat.data = _resample_dataframe!(dat.data, factor, :triggers)

    # For continuous/ERP data, renumber sample column to be sequential
    if hasproperty(dat.data, :sample)
        dat.data.sample = 1:nrow(dat.data)
    end

    # Update sample rate
    dat.sample_rate = dat.sample_rate ÷ factor

    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data)) samples"

    return nothing
end



"""
    resample!(dat::EpochData, factor::Int)::Nothing

Resample epoched data in-place by downsampling each epoch by the specified factor.

This function reduces the sampling rate of each epoch by keeping every nth sample,
where n is the downsampling factor. All metadata columns (triggers, trial info,
conditions, etc.) are preserved in each epoch. The sample rate is updated accordingly.

# Arguments
- `dat::EpochData`: Epoched data to resample
- `factor::Int`: Downsampling factor (e.g., 2 = keep every 2nd sample, halving the sample rate)

# Effects
- Modifies the input data in-place
- Updates the sample rate: new_rate = old_rate / factor
- Resamples all epochs
- Preserves all metadata columns in each epoch

# Examples
```julia
using EegFun

# Load epoched data at 512 Hz
epochs = load("participant_1_epochs.jld2", "epochs")

# Downsample to 256 Hz (factor of 2)
resample!(epochs, 2)
```

# Notes
- Factor must be a positive integer
- New sample rate must be an integer (old_rate must be divisible by factor)
- Simple decimation is used (no anti-aliasing filter applied)
- For proper downsampling, low-pass filter data first to avoid aliasing
- All DataFrame columns in each epoch are preserved
- Triggers and metadata are maintained for each epoch
"""
function resample!(dat::EpochData, factor::Int)::Nothing
    # Validation
    if factor < 1
        @minimal_error_throw("Downsampling factor must be positive, got $factor")
    end

    if factor == 1
        @info "Downsampling factor is 1, no resampling needed"
        return nothing
    end

    if dat.sample_rate % factor != 0
        @minimal_error_throw(
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



"""
    resample(dat::ContinuousData, factor::Int)

Non-mutating version of resample!. Returns a new data object with resampled data.

# Examples
```julia
# Downsample without modifying original
data_256hz = resample(data_512hz, 2)

# Original data unchanged
@assert data_512hz.sample_rate == 512
@assert data_256hz.sample_rate == 256
```
"""
function resample(dat::T, factor::Int)::T where {T<:EegData}
    dat_copy = copy(dat)
    resample!(dat_copy, factor)
    return dat_copy
end

"""
    resample(data_vec::Vector{T}, factor::Int)::Vector{T} where T <: EegData

Resample a vector of EEG data objects (e.g., multiple participants or conditions).
Returns a new vector with each element resampled.

# Examples
```julia
# Resample multiple participants' data
resampled_data = resample(all_participants, 2)
```
"""
function resample(data_vec::Vector{T}, factor::Int)::Vector{T} where {T<:EegData}
    return [resample(dat, factor) for dat in data_vec]
end

"""
    resample!(data_vec::Vector{T}, factor::Int)::Nothing where {T<:EegData}

Mutating version of resample for vector of EEG data objects.

# Arguments
- `data_vec::Vector{T}`: Vector of EEG data objects to resample (modified in-place)
- `factor::Int`: Downsampling factor (must be positive integer)

# Returns
- `Nothing`: All objects in the vector are modified in-place

# Examples
```julia
# Resample multiple EpochData objects
resample!(epochs_vector, 2)  # Downsample by factor of 2
```
"""
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

    # Load data using load_data (handles single variable files automatically)
    loaded_data = load_data(filepath)

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

        # Save results (always use "data" as variable name since load_data finds by type)
        jldsave(output_path; data = resampled_data)

        message = "Resampled from $old_rate Hz to $new_rate Hz (factor: $factor)"
        return BatchResult(true, filename, message)
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end


"""
    resample(file_pattern::String, factor::Int;
            input_dir::String = pwd(),
            participant_selection::Function = participants(),
            output_dir::Union{String, Nothing} = nothing)

Batch resample data files and save to a new directory.

This function processes multiple participant files at once, downsampling each
by the specified factor.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "continuous", "epochs", "erp")
- `factor::Int`: Downsampling factor
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory)

# Examples
```julia
# Downsample all continuous data files by factor of 2
resample("continuous", 2)

# Specific participants only
resample("epochs", 4, participants = [1, 2, 3])

# Custom output directory
resample("erp", 2, output_dir = "/path/to/output")

# Full example workflow
resample("continuous", 2,
        input_dir = "/data/study1",
        participants = 1:20)
```

# Output
- Creates new directory with resampled data files
- Each output file contains the same variable name as input
- Log file saved to output directory

# Notes
- All input files should have the same sample rate
- Factor must result in integer sample rates for all files
- For proper downsampling, low-pass filter data first to avoid aliasing
- Typical factors: 2, 4, 8 (powers of 2 work best)
"""
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
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if factor <= 0
            @minimal_error_throw("Downsampling factor must be positive, got $factor")
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
