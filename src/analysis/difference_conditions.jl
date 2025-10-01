"""
Batch computation of condition difference waves for ERP data.
"""

#=============================================================================
    DIFFERENCE-SPECIFIC VALIDATION
=============================================================================#

"""Validate that file pattern is for ERP data."""
function _validate_erps_pattern(pattern::String)
    !contains(pattern, "erps") && 
        return "difference_conditions_data only works with ERP data. File pattern must contain 'erps', got: '$pattern'"
    return nothing
end

"""Generate default output directory name for difference operation."""
function _default_difference_output_dir(input_dir::String, pattern::String, pairs::Union{Vector{Tuple{Int, Int}}, Vector{Vector{Int}}})
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
    if nrow(erp1.data) != nrow(erp2.data)
        @minimal_error_throw("ERPs have different numbers of time points: $(nrow(erp1.data)) vs $(nrow(erp2.data))")
    end
    
    if erp1.sample_rate != erp2.sample_rate
        @minimal_error_throw("ERPs have different sample rates: $(erp1.sample_rate) vs $(erp2.sample_rate)")
    end
    
    # Get EEG channels (exclude metadata columns)
    metadata_cols = meta_labels(erp1)
    eeg_channels = setdiff(propertynames(erp1.data), metadata_cols)
    
    # Create a copy of erp1's data for the difference
    diff_data = copy(erp1.data)
    
    # Update condition information
    diff_data.condition .= diff_cond  # Sequential condition number for differences
    diff_data.condition_name .= "difference_$(cond1)_$(cond2)"
    
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
    
    return ErpData(diff_data, erp1.layout, erp1.sample_rate, erp1.analysis_info, min_epochs)
end

"""
Process a single ERP file through difference wave creation.
Returns BatchResult with success/failure info.
"""
function _process_difference_file(filepath::String, output_path::String,
                                 condition_pairs::Union{Vector{Tuple{Int, Int}}, Vector{Vector{Int}}})
    filename = basename(filepath)
    
    # Load data
    file_data = load(filepath)
    
    if !haskey(file_data, "erps")
        return BatchResult(false, filename, "No 'erps' variable found")
    end
    
    erps_data = file_data["erps"]
    
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
    save(output_path, "differences", difference_waves)
    
    return BatchResult(true, filename, "Created $(length(difference_waves)) difference wave(s)")
end

#=============================================================================
    MAIN API FUNCTION
=============================================================================#

"""
    difference_conditions(file_pattern::String, condition_pairs; 
                         input_dir::String = pwd(), 
                         participants::Union{Int, Vector{Int}, Nothing} = nothing,
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
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Examples
```julia
# Create difference waves for conditions 1-2 and 3-4 (tuples)
difference_conditions("erps_cleaned", [(1,2), (3,4)])

# Or with vectors
difference_conditions("erps_cleaned", [[1,2], [3,4]])

# Process specific participants
difference_conditions("erps_cleaned", [(1,2)], 
                     input_dir = "/path/to/data", 
                     participants = [1, 2, 3])

# Specify custom output directory
difference_conditions("erps_cleaned", [(1,2), (3,4)], 
                     input_dir = "/path/to/data", 
                     output_dir = "/path/to/output")
```
"""
function difference_conditions(file_pattern::String, condition_pairs::Union{Vector{Tuple{Int, Int}}, Vector{Vector{Int}}}; 
                                   input_dir::String = pwd(), 
                                   participants::Union{Int, Vector{Int}, Nothing} = nothing,
                                   output_dir::Union{String, Nothing} = nothing)
    
    # Setup logging
    log_file = "difference_conditions.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch condition differencing started at $(now())"
        @log_call "difference_conditions" (file_pattern, condition_pairs)
        
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
        output_dir = something(output_dir, _default_difference_output_dir(input_dir, file_pattern, condition_pairs))
        mkpath(output_dir)
        
        # Find files
        files = _find_batch_files(file_pattern, input_dir; participants)
        
        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            return nothing
        end
        
        @info "Found $(length(files)) JLD2 files matching pattern '$file_pattern'"
        @info "Condition pairs: $condition_pairs"
        
        # Create processing function with captured parameters
        process_fn = (input_path, output_path) -> 
            _process_difference_file(input_path, output_path, condition_pairs)
        
        # Execute batch operation
        results = _run_batch_operation(process_fn, files, input_dir, output_dir; 
                                      operation_name="Creating differences")
        
        _log_batch_summary(results, output_dir)
        
    finally
        _cleanup_logging(log_file, output_dir)
    end
end

"""
Batch computation of condition difference waves for ERP data.
"""
