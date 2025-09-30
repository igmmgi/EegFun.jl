"""
    difference_conditions_data(file_pattern::String, condition_pairs::Vector{Tuple{Int, Int}}; 
                              input_dir::String = pwd(), 
                              participants::Union{Int, Vector{Int}, Nothing} = nothing,
                              output_dir::Union{String, Nothing} = nothing)

Batch process ERP data files to create condition difference waves.

This function loads JLD2 files containing ERP data, computes differences between specified condition pairs
by subtracting EEG channel columns, and saves the resulting difference waves to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `condition_pairs::Vector{Tuple{Int, Int}}`: Pairs of conditions to subtract (e.g., [(1,2), (3,4)])
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Returns
- `Nothing`: Saves difference wave data to output directory

# Examples
```julia
# Create difference waves for conditions 1-2 and 3-4
difference_conditions_data("erps_cleaned", [(1,2), (3,4)])

# Process specific participants
difference_conditions_data("erps_cleaned", [(1,2)], 
                          input_dir = "/path/to/data", 
                          participants = [1, 2, 3])

# Specify custom output directory
difference_conditions_data("erps_cleaned", [(1,2), (3,4)], 
                          input_dir = "/path/to/data", 
                          output_dir = "/path/to/output")
```
"""
function difference_conditions_data(file_pattern::String, condition_pairs::Vector{Tuple{Int, Int}}; 
                                   input_dir::String = pwd(), 
                                   participants::Union{Int, Vector{Int}, Nothing} = nothing,
                                   output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "difference_conditions_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch condition differencing started at $(now())"
        
        # Log the function call
        @log_call "difference_conditions_data" (file_pattern, condition_pairs)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        @info "Condition pairs: $condition_pairs"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        # Validate that file pattern contains "erps"
        if !contains(file_pattern, "erps")
            @minimal_error_throw("difference_conditions_data only works with ERP data. File pattern must contain 'erps', got: '$file_pattern'")
        end
        
        if isempty(condition_pairs)
            @minimal_error_throw("Condition pairs cannot be empty")
        end
        
        # Validate condition pairs
        for (i, (cond1, cond2)) in enumerate(condition_pairs)
            if cond1 == cond2
                @minimal_warning "Condition pair $i: ($cond1, $cond2) has identical conditions. Difference will be zero."
            end
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            pairs_str = join([join(pair, "-") for pair in condition_pairs], "_")
            output_dir = joinpath(input_dir, "differences_$(file_pattern)_$(pairs_str)")
        end
        
        # Create output directory if it doesn't exist
        mkpath(output_dir)
        
        # Find JLD2 files matching the pattern
        all_files = readdir(input_dir)
        jld2_files = filter(x -> endswith(x, ".jld2") && contains(x, file_pattern), all_files)
        
        # Filter by participant number if specified
        if participants !== nothing
            jld2_files = _filter_files(jld2_files; include = participants)
        end
        
        if isempty(jld2_files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            @info "Available files: $(filter(x -> endswith(x, ".jld2"), all_files))"
            return nothing
        end
        
        @info "Found $(length(jld2_files)) JLD2 files matching pattern '$file_pattern'"
        
        processed_count = 0
        error_count = 0
        
        for (i, file) in enumerate(jld2_files)
            input_path = joinpath(input_dir, file)
            output_path = joinpath(output_dir, file)
            
            @info("Processing: $file ($(i)/$(length(jld2_files)))")
            
            try
                # Load the file and determine what data we have
                file_data = load(input_path)
                
                # Check what type of data we have
                if haskey(file_data, "erps")
                    erps_data = file_data["erps"]
                    var_name = "erps"
                else
                    @minimal_warning "No 'erps' variable found in $file. Available: $(keys(file_data))"
                    continue
                end
                
                # Get available conditions
                available_conditions = [condition_number(erp) for erp in erps_data]
                @info "  Found conditions: $available_conditions"
                
                # Create difference waves for each condition pair
                difference_waves = ErpData[]
                for (pair_idx, (cond1, cond2)) in enumerate(condition_pairs)
                    @info "  Creating difference wave $pair_idx: condition $cond1 - condition $cond2"
                    
                    # Find the ERP data for each condition
                    erp1 = nothing
                    erp2 = nothing
                    
                    for erp in erps_data
                        if condition_number(erp) == cond1
                            erp1 = erp
                        elseif condition_number(erp) == cond2
                            erp2 = erp
                        end
                    end
                    
                    # Check if both conditions exist
                    if erp1 === nothing
                        @minimal_warning "  Condition $cond1 not found in data. Skipping pair ($cond1, $cond2)."
                        continue
                    end
                    if erp2 === nothing
                        @minimal_warning "  Condition $cond2 not found in data. Skipping pair ($cond1, $cond2)."
                        continue
                    end
                    
                    # Create difference wave with sequential condition number
                    diff_wave = _create_difference_wave(erp1, erp2, cond1, cond2, pair_idx)
                    push!(difference_waves, diff_wave)
                    
                    @info "    Created difference wave with $(nrow(diff_wave.data)) time points"
                end
                
                if isempty(difference_waves)
                    @minimal_warning "  No valid condition pairs found in $file. Skipping."
                    continue
                end
                
                # Save difference waves
                save(output_path, "differences", difference_waves)
                
                @info "  Saved: $output_path"
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                error_count += 1
                continue
            end
        end
        
        @info "Condition differencing complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "difference_conditions_data.log"
        log_dest = joinpath(output_dir, "difference_conditions_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end

"""
    _create_difference_wave(erp1::ErpData, erp2::ErpData, cond1::Int, cond2::Int, diff_cond::Int) -> ErpData

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
