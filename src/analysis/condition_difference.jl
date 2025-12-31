"""
Batch computation of condition difference waves for ERP data.
"""

#=============================================================================
    DIFFERENCE-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for ERP data."""
function _validate_erps_pattern(pattern::String)
    !contains(pattern, "erps") &&
        return "condition_difference_data only works with ERP data. File pattern must contain 'erps', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for difference operation."""
function _condition_difference_default_output_dir(
    input_dir::String,
    pattern::String,
    pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)
    pairs_str = join([join(pair, "-") for pair in pairs], "_")
    joinpath(input_dir, "differences_$(pattern)_$(pairs_str)")
end

#=============================================================================
    DIFFERENCE-SPECIFIC PROCESSING
=============================================================================#

"""
Create a difference wave by subtracting ERP2 from ERP1.
"""
function _create_difference_wave(erp1::ErpData, erp2::ErpData, cond1::Int, cond2::Int, diff_cond::Int)
    # Validate that both ERPs have the same structure
    have_same_structure(erp1, erp2) || @minimal_error_throw("ERPs have inconsistent structure")

    # Get EEG channels (exclude metadata columns)
    metadata_cols = meta_labels(erp1)
    eeg_channels = setdiff(propertynames(erp1.data), metadata_cols)

    # Create a copy of erp1's data for the difference
    diff_data = copy(erp1.data)

    # Remove condition/condition_name/n_epochs columns if they exist (they're in struct now)
    cols_to_remove = [:condition, :condition_name, :n_epochs]
    for col in cols_to_remove
        if hasproperty(diff_data, col)
            select!(diff_data, Not(col))
        end
    end

    # Subtract EEG channels
    for ch in eeg_channels
        if hasproperty(erp2.data, ch)
            diff_data[!, ch] = erp1.data[!, ch] .- erp2.data[!, ch]
        else
            @minimal_warning "Channel $ch not found in condition $cond2, keeping original values"
        end
    end

    # Update n_epochs to reflect the minimum (conservative estimate)
    min_epochs = min(erp1.n_epochs, erp2.n_epochs)
    diff_condition_name = "difference_$(cond1)_$(cond2)"

    return ErpData(
        erp1.file,
        diff_cond,
        diff_condition_name,
        diff_data,
        erp1.layout,
        erp1.sample_rate,
        erp1.analysis_info,
        min_epochs,
    )
end

"""
Process a single ERP file through difference wave creation.
Returns BatchResult with success/failure info.
"""
function _condition_difference_process_file(
    filepath::String,
    output_path::String,
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}},
)
    filename = basename(filepath)

    # Load data (using load_data which finds by type)
    erps_data = load_data(filepath)

    if isnothing(erps_data)
        return BatchResult(false, filename, "No data variables found")
    end

    # Validate that data is Vector{ErpData}
    if !(erps_data isa Vector{<:ErpData})
        return BatchResult(false, filename, "Invalid data type: expected Vector{ErpData}, got $(typeof(erps_data))")
    end

    # Create difference waves for each condition pair
    difference_waves = ErpData[]

    for (pair_idx, (cond1, cond2)) in enumerate(condition_pairs)
        # Find the ERP data for each condition
        erp1 = nothing
        erp2 = nothing

        for erp in erps_data
            # Extract condition number from ErpData.data DataFrame
            cond_num = erp.data[1, :condition]
            if cond_num == cond1
                erp1 = erp
            elseif cond_num == cond2
                erp2 = erp
            end
        end

        # Check if both conditions exist
        if isnothing(erp1)
            @minimal_warning "  Condition $cond1 not found in $filename. Skipping pair ($cond1, $cond2)."
            continue
        end
        if isnothing(erp2)
            @minimal_warning "  Condition $cond2 not found in $filename. Skipping pair ($cond1, $cond2)."
            continue
        end

        # Create difference wave with sequential condition number
        diff_wave = _create_difference_wave(erp1, erp2, cond1, cond2, pair_idx)
        push!(difference_waves, diff_wave)
    end

    if isempty(difference_waves)
        return BatchResult(false, filename, "No valid condition pairs found")
    end

    # Save difference waves
    jldsave(output_path; data = difference_waves)

    return BatchResult(true, filename, "Created $(length(difference_waves)) difference wave(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    condition_difference(file_pattern::String, condition_pairs; 
                         input_dir::String = pwd(), 
                         participant_selection::Function = participants(),
                         output_dir::Union{String, Nothing} = nothing)

Batch process ERP data files to create condition difference waves.

This function loads JLD2 files containing ERP data, computes differences between specified condition pairs
by subtracting EEG channel columns, and saves the resulting difference waves to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `condition_pairs`: Pairs of conditions to subtract. Can be:
  - `Vector{Tuple{Int, Int}}`: e.g., `[(1,2), (3,4)]`
  - `Vector{Vector{Int}}`: e.g., `[[1,2], [3,4]]`
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participant_selection::Function`: Participant selection predicate (default: `participants()` for all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Create difference waves for conditions 1-2 and 3-4 (tuples)
condition_difference("erps_cleaned", [(1,2), (3,4)])

# Or with vectors
condition_difference("erps_cleaned", [[1,2], [3,4]])

# Process specific participants
condition_difference("erps_cleaned", [(1,2)], 
                     input_dir = "/path/to/data", 
                     participants = [1, 2, 3])

# Specify custom output directory
condition_difference("erps_cleaned", [(1,2), (3,4)], 
                     input_dir = "/path/to/data", 
                     output_dir = "/path/to/output")
```
"""
function condition_difference(
    file_pattern::String,
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}};
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
)

    # Setup logging
    log_file = "condition_difference.log"
    setup_global_logging(log_file)

    result = (success = 0, errors = 0)  # Default return value
    try
        @info "Batch condition differencing started at $(now())"
        @log_call "condition_difference"

        # Validation (early return on error)
        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_erps_pattern(file_pattern)) !== nothing
            @minimal_error_throw(error_msg)
        end

        if (error_msg = _validate_condition_pairs(condition_pairs)) !== nothing
            @minimal_error_throw(error_msg)
        end

        # Setup directories
        output_dir =
            something(output_dir, _condition_difference_default_output_dir(input_dir, file_pattern, condition_pairs))
        mkpath(output_dir)

        # Find files
        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            result = (success = 0, errors = 0)
        else
            @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
            @info "Condition pairs: $condition_pairs"

            # Create processing function with captured parameters
            process_fn =
                (input_path, output_path) ->
                    _condition_difference_process_file(input_path, output_path, condition_pairs)

            # Execute batch operation
            results =
                _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Creating differences")

            result = _log_batch_summary(results, output_dir)
        end

    finally
        _cleanup_logging(log_file, output_dir)
    end

    return result
end

"""
Batch computation of condition difference waves for ERP data.
"""

# ============================================================================
# TF Difference for TimeFreqData
# ============================================================================

"""
    tf_difference(tf1::TimeFreqData, tf2::TimeFreqData) -> TimeFreqData

Compute the difference between two TimeFreqData objects (tf1 - tf2).

# Arguments
- `tf1::TimeFreqData`: First TF data (minuend)
- `tf2::TimeFreqData`: Second TF data (subtrahend)

# Returns
- `TimeFreqData`: Difference wave (tf1 - tf2)

# Example
```julia
diff_tf = tf_difference(tf_cond1, tf_cond2)
```
"""
function tf_difference(tf1::TimeFreqData, tf2::TimeFreqData)
    # Validate same structure
    times1, times2 = unique(tf1.data.time), unique(tf2.data.time)
    freqs1, freqs2 = unique(tf1.data.freq), unique(tf2.data.freq)
    ch1, ch2 = channel_labels(tf1), channel_labels(tf2)
    
    times1 == times2 || error("TimeFreqData objects have different time vectors")
    freqs1 == freqs2 || error("TimeFreqData objects have different frequency vectors")
    ch1 == ch2 || error("TimeFreqData objects have different channels")
    
    # Create difference DataFrame
    diff_data = copy(tf1.data)
    for ch in ch1
        diff_data[!, ch] = tf1.data[!, ch] .- tf2.data[!, ch]
    end
    
    diff_name = "$(tf1.condition_name)_minus_$(tf2.condition_name)"
    
    return TimeFreqData(
        tf1.file,
        tf1.condition * 100 + tf2.condition,  # Combined condition number
        diff_name,
        diff_data,
        tf1.layout,
        tf1.sample_rate,
        tf1.method,
        nothing,  # Baseline info not preserved for condition differences
        tf1.analysis_info
    )
end

"""
    tf_difference(file_pattern::String, condition_pairs;
                  input_dir=pwd(), output_dir=nothing,
                  participant_selection=participants())

Compute condition difference waves for TimeFreqData files.

# Arguments
- `file_pattern::String`: Pattern to match files
- `condition_pairs`: Vector of (cond1, cond2) pairs to compute cond1 - cond2

# Example
```julia
tf_difference("tf_epochs_wavelet", [(1, 2), (3, 4)])
```
"""
function tf_difference(
    file_pattern::String,
    condition_pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}};
    input_dir::String=pwd(),
    output_dir::Union{String,Nothing}=nothing,
    participant_selection::Function=participants()
)
    log_file = "tf_difference.log"
    setup_global_logging(log_file)

    try
        @info "Batch TF difference started at $(now())"
        @log_call "tf_difference"

        if (error_msg = _validate_input_dir(input_dir)) !== nothing
            @minimal_error_throw(error_msg)
        end

        output_dir = something(output_dir, joinpath(input_dir, "tf_difference_$(file_pattern)"))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern'"
            return nothing
        end

        @info "Found $(length(files)) files, condition pairs: $condition_pairs"

        process_fn = (input_path, output_path) -> begin
            filename = basename(input_path)
            data = load_data(input_path)
            if isnothing(data) || !(data isa Vector{TimeFreqData})
                return BatchResult(false, filename, "Invalid data type")
            end
            
            diff_results = TimeFreqData[]
            for (cond1, cond2) in condition_pairs
                tf1 = findfirst(tf -> tf.condition == cond1, data)
                tf2 = findfirst(tf -> tf.condition == cond2, data)
                
                if isnothing(tf1) || isnothing(tf2)
                    @minimal_warning "Missing condition $cond1 or $cond2 in $filename"
                    continue
                end
                
                diff_tf = tf_difference(data[tf1], data[tf2])
                push!(diff_results, diff_tf)
            end
            
            if isempty(diff_results)
                return BatchResult(false, filename, "No differences computed")
            end
            
            jldsave(output_path; data=diff_results)
            return BatchResult(true, filename, "Created $(length(diff_results)) difference(s)")
        end

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                       operation_name="TF difference")
        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
end
