"""
    combine_conditions(file_pattern::String, condition_groups::Vector{Vector{Int}}; 
                      input_dir::String = pwd(), 
                      participants::Union{Int, Vector{Int}, Nothing} = nothing,
                      output_dir::Union{String, Nothing} = nothing)

Combine EEG epoch conditions from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files (should contain "epochs")
- `condition_groups::Vector{Vector{Int}}`: Groups of condition numbers to combine (e.g., [[1,2], [3,4]])
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on groups)

# Example
```julia
# Combine conditions 1,2 into first group and 3,4 into second group
combine_conditions("epochs", [[1, 2], [3, 4]])

# Combine specific participant
combine_conditions("epochs_cleaned", [[1, 2], [3, 4]], participants=3)

# Combine with custom output directory
combine_conditions("epochs", [[1, 2], [3, 4]], output_dir="/path/to/output/")
```

# Note
- Only works with epoch data (not ERPs)
- Conditions are combined (concatenated) into new conditions
- Use `average_epochs()` separately to create ERPs from combined epochs
"""
function combine_conditions(file_pattern::String, condition_groups::Vector{Vector{Int}}; 
                           input_dir::String = pwd(), 
                           participants::Union{Int, Vector{Int}, Nothing} = nothing,
                           output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "combine_conditions.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch condition combining started at $(now())"
        
        # Log the function call generically
        args = [file_pattern, condition_groups]
        kwargs = [:input_dir => input_dir, :participants => participants, :output_dir => output_dir]
        _log_function_call("combine_conditions", args, kwargs)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        @info "Condition groups: $condition_groups"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        if isempty(condition_groups)
            @minimal_error_throw("Condition groups cannot be empty")
        end
        
        # Validate that file pattern contains "epochs"
        if !contains(file_pattern, "epochs")
            @minimal_error_throw("combine_conditions only works with epoch data. File pattern must contain 'epochs', got: '$file_pattern'")
        end
        
        # Validate condition groups for duplicates and overlaps
        all_conditions = Int[]
        for (group_idx, group) in enumerate(condition_groups)
            # Check for duplicates within the group
            if length(group) != length(unique(group))
                duplicates = group[findall(x -> count(==(x), group) > 1, group)]
                @minimal_warning "Group $group_idx contains duplicate conditions: $duplicates. Only unique conditions will be used."
            end
            
            # Remove duplicates and update the group
            unique_group = unique(group)
            if length(unique_group) != length(group)
                @info "  Group $group_idx: $(group) → $(unique_group) (removed duplicates)"
            end
            condition_groups[group_idx] = unique_group
            
            # Check for overlaps between groups
            overlap = intersect(all_conditions, unique_group)
            if !isempty(overlap)
                @minimal_warning "Condition(s) $overlap appear in multiple groups. This may not be intended."
            end
            
            append!(all_conditions, unique_group)
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            groups_str = join([join(group, "-") for group in condition_groups], "_")
            output_dir = joinpath(input_dir, "combined_$(file_pattern)_$(groups_str)")
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
                # Load epochs data
                file_data = load(input_path)
                
                if !haskey(file_data, "epochs")
                    @minimal_warning "No 'epochs' variable found in $file. Available: $(keys(file_data))"
                    continue
                end
                
                data = file_data["epochs"]
                @info "  Found $(length(data)) conditions in data"
                
                # Combine conditions for epochs
                combined_data = Vector{Any}()
                
                for (group_idx, original_conditions) in enumerate(condition_groups)
                    @info "  Creating condition group $group_idx from conditions $original_conditions"
                    
                    # Validate that all requested conditions exist
                    max_condition = length(data)
                    missing_conditions = filter(c -> c > max_condition || c < 1, original_conditions)
                    if !isempty(missing_conditions)
                        @error "  Condition(s) $missing_conditions not found in data (only has conditions 1-$max_condition)"
                        continue
                    end
                    
                    # Get data for the specified conditions
                    condition_data = data[original_conditions]
                    
                    # For epochs: combine (concatenate) the conditions
                    # Each condition is an EpochData object, we need to concatenate their data fields
                    combined_data_frames = vcat([epoch_data.data for epoch_data in condition_data]...)
                    
                    # Create new EpochData with concatenated data, using metadata from first condition
                    combined_epochs = EpochData(
                        combined_data_frames,
                        condition_data[1].layout,
                        condition_data[1].sample_rate,
                        condition_data[1].analysis_info
                    )
                    
                    push!(combined_data, combined_epochs)
                    @info "    Combined $(length(condition_data)) conditions into $(length(combined_data_frames)) epochs"
                end
                
                # Save combined data
                save(output_path, "epochs", combined_data)
                
                @info "  Saved: $output_path"
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                error_count += 1
                continue
            end
        end
        
        @info "Condition combining complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "combine_conditions.log"
        log_dest = joinpath(output_dir, "combine_conditions.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end
