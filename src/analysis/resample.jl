"""
Resampling of EEG data to different sampling rates.

This function downsamples EEG data by a specified factor, reducing the sampling
rate while preserving all metadata including triggers, trial information, and
other columns. Works with continuous, epoched, and ERP data.
"""

#=============================================================================
    CORE RESAMPLING FUNCTIONS
=============================================================================#

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
using eegfun

# Load continuous data at 512 Hz
data = load("participant_1_continuous.jld2", "data")

# Downsample to 256 Hz (factor of 2)
resample!(data, 2)

# Now data.sample_rate == 256
```

# Notes
- Factor must be a positive integer
- New sample rate must be an integer (old_rate must be divisible by factor)
- Simple decimation is used (no anti-aliasing filter applied)
- For proper downsampling, low-pass filter data first to avoid aliasing
- All DataFrame columns are preserved, including time, triggers, and metadata

# Use Cases
- Reduce data size for faster processing
- Match sampling rates across different recordings
- Prepare data for algorithms requiring specific sampling rates
"""
function resample!(dat::Union{ContinuousData,ErpData}, factor::Int)::Nothing
    # Validation
    if factor <= 0
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
   
     # Get indices of samples to keep
     base_indices = 1:factor:nrow(dat.data)
    
     # Find indices of samples with triggers (if trigger column exists)
     if hasproperty(dat.data, :trigger)
         trigger_indices = findall(dat.data.trigger .!= 0)
         all_indices = sort(unique(vcat(base_indices, trigger_indices)))
     else
         all_indices = base_indices
     end
     
     # Downsample by keeping selected indices
     dat.data = dat.data[collect(all_indices), :]

   
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
using eegfun

# Load epoched data at 512 Hz
epochs = load("participant_1_epochs.jld2", "epochs")

# Downsample to 256 Hz (factor of 2)
resample!(epochs, 2)

# Now epochs.sample_rate == 256
```

# Notes
- Factor must be a positive integer
- New sample rate must be an integer (old_rate must be divisible by factor)
- Simple decimation is used (no anti-aliasing filter applied)
- For proper downsampling, low-pass filter data first to avoid aliasing
- All DataFrame columns in each epoch are preserved
- Triggers and metadata are maintained for each epoch

# Use Cases
- Reduce epoch data size for faster processing
- Match sampling rates for cross-study comparisons
- Prepare epoched data for machine learning pipelines
"""
function resample!(dat::EpochData, factor::Int)::Nothing
    # Validation
    if factor <= 0
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
        # Get indices of samples to keep
        base_indices = 1:factor:nrow(epoch)
        
        # Find indices of samples with triggers (if trigger column exists)
        if hasproperty(epoch, :trigger)
            trigger_indices = findall(epoch.trigger .!= 0)
            all_indices = sort(unique(vcat(base_indices, trigger_indices)))
        else
            all_indices = base_indices
        end
        
        # Downsample by keeping selected indices
        dat.data[i] = epoch[collect(all_indices), :]
    end
    
    # Update sample rate
    dat.sample_rate = dat.sample_rate ÷ factor
    
    @info "Resampling complete. New sample rate: $(dat.sample_rate) Hz, $(nrow(dat.data[1])) samples per epoch"
    
    return nothing
end



"""
    resample(dat::Union{ContinuousData, ErpData, EpochData}, factor::Int)

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
function resample(dat::ContinuousData, factor::Int)::ContinuousData
    # Create a deep copy
    dat_copy = ContinuousData(
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info)
    )
    
    # Apply resampling
    resample!(dat_copy, factor)
    
    return dat_copy
end

function resample(dat::ErpData, factor::Int)::ErpData
    # Create a deep copy
    dat_copy = ErpData(
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
        dat.n_epochs
    )
    
    # Apply resampling
    resample!(dat_copy, factor)
    
    return dat_copy
end

function resample(dat::EpochData, factor::Int)::EpochData
    # Create a deep copy
    dat_copy = EpochData(
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info)
    )
    
    # Apply resampling
    resample!(dat_copy, factor)
    
    return dat_copy
end


#=============================================================================
    BATCH PROCESSING FUNCTIONS
=============================================================================#

"""Generate default output directory name for resampling operation."""
function _default_resample_output_dir(input_dir::String, pattern::String, factor::Int)
    new_rate_str = "resampled_by_$(factor)"
    joinpath(input_dir, "$(new_rate_str)_$(pattern)")
end


"""
Process a single data file through resampling pipeline.
Returns BatchResult with success/failure info.
"""
function _process_resample_file(
    filepath::String,
    output_path::String,
    factor::Int,
)
    filename = basename(filepath)
    
    # Load data
    file_data = load(filepath)
    
    # Try common variable names
    var_names = ["data", "epochs", "erp", "continuous", "epoch_data", "erp_data", "continuous_data"]
    loaded_data = nothing
    data_var_name = nothing
    
    for var_name in var_names
        if haskey(file_data, var_name)
            loaded_data = file_data[var_name]
            data_var_name = var_name
            break
        end
    end
    
    if isnothing(loaded_data)
        return BatchResult(false, filename, "No EEG data variable found (tried: $(var_names))")
    end
    
    if !(loaded_data isa Union{ContinuousData,EpochData,ErpData})
        return BatchResult(false, filename, "Data is not a recognized EEG data type")
    end
    
    # Resample
    try
        old_rate = loaded_data.sample_rate
        resampled_data = resample(loaded_data, factor)
        new_rate = resampled_data.sample_rate
        
        # Save results using original variable name
        save(output_path, data_var_name, resampled_data)
        
        message = "Resampled from $old_rate Hz to $new_rate Hz (factor: $factor)"
        return BatchResult(true, filename, message)
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end


"""
    resample(file_pattern::String, factor::Int;
            input_dir::String = pwd(),
            participants::Union{Int, Vector{Int}, Nothing} = nothing,
            output_dir::Union{String, Nothing} = nothing)

Batch resample data files and save to a new directory.

This function processes multiple participant files at once, downsampling each
by the specified factor.

# Arguments
- `file_pattern::String`: Pattern to match files (e.g., "continuous", "epochs", "erp")
- `factor::Int`: Downsampling factor
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
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
    participants::Union{Int,Vector{Int},Nothing} = nothing,
    output_dir::Union{String,Nothing} = nothing,
)
    
    # Setup logging
    log_file = "resample.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch resampling started at $(now())"
        @log_call "resample" (file_pattern, factor)
        
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
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Downsampling factor: $factor"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) ->
            _process_resample_file(input_path, output_path, factor)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Resampling")
        
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end

