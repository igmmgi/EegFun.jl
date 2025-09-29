"""
    combine_channels_data(file_pattern::String, channel_groups::Vector{Vector{Symbol}}; 
                         input_dir::String = pwd(), 
                         participants::Union{Int, Vector{Int}, Nothing} = nothing,
                         conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                         output_dir::Union{String, Nothing} = nothing,
                         reduce::Bool = false)

Batch process EEG data files to combine specified channels.

This function loads JLD2 files containing EEG data, combines specified channel groups by averaging,
and saves the resulting data with new combined channels to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "epochs_cleaned")
- `channel_groups::Vector{Vector{Symbol}}`: Groups of channels to combine (e.g., `[[:Fp1, :Fp2], [:PO7, :PO8]]`)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)
- `reduce::Bool`: Whether to keep only combined channels (true) or append to existing (false, default)

# Returns
- `Nothing`: Saves data with combined channels to output directory

# Examples
```julia
# Combine frontal and parietal channels
combine_channels_data("erps_cleaned", [[:Fp1, :Fp2], [:PO7, :PO8]])

# Process specific participants and conditions
combine_channels_data("erps_cleaned", [[:Fp1, :Fp2]], 
                     input_dir = "/path/to/data", 
                     participants = [1, 2, 3], 
                     conditions = [1, 2])

# Create reduced dataset with only combined channels
combine_channels_data("erps_cleaned", [[:Fp1, :Fp2], [:PO7, :PO8]], 
                     reduce = true)
```
"""
function combine_channels_data(file_pattern::String, channel_groups::Vector{Vector{Symbol}}; 
                              input_dir::String = pwd(), 
                              participants::Union{Int, Vector{Int}, Nothing} = nothing,
                              conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                              output_dir::Union{String, Nothing} = nothing,
                              reduce::Bool = false)
    
    # Set up global logging
    log_file = "combine_channels_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch channel combining started at $(now())"
        
        # Log the function call generically
        args = [file_pattern, channel_groups]
        kwargs = [:input_dir => input_dir, :participants => participants, :conditions => conditions, :output_dir => output_dir, :reduce => reduce]
        _log_function_call("combine_channels_data", args, kwargs)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        @info "Channel groups: $channel_groups"
        @info "Reduce mode: $reduce"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        if isempty(channel_groups)
            @minimal_error_throw("Channel groups cannot be empty")
        end
        
        # Validate channel groups
        for (i, group) in enumerate(channel_groups)
            if isempty(group)
                @minimal_error_throw("Channel group $i is empty")
            end
            if length(group) < 2
                @minimal_warning "Channel group $i has only $(length(group)) channel(s): $group. Consider using more channels for meaningful averaging."
            end
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            groups_str = join([join(group, "_") for group in channel_groups], "_")
            output_dir = joinpath(input_dir, "combined_channels_$(file_pattern)_$(groups_str)")
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
                possible_vars = ["erps", "epochs"]
                data = nothing
                var_name = nothing
                
                for var in possible_vars
                    if haskey(file_data, var)
                        data = file_data[var]
                        var_name = var
                        break
                    end
                end
                
                if data === nothing
                    @minimal_warning "No recognized data variable found in $file. Available: $(keys(file_data))"
                    continue
                end
                
                # Filter by condition if specified
                if conditions !== nothing
                    condition_nums = conditions isa Int ? [conditions] : conditions
                    # Filter data to only include specified conditions
                    data = data[condition_nums]
                    @info "  Filtering conditions: $condition_nums"
                end
                
                # Combine channels for each data item
                for (j, item) in enumerate(data)
                    @info "  Combining channels $j/$(length(data))"
                    
                    # Create channel selection functions from groups
                    channel_selections = [channels(group) for group in channel_groups]
                    
                    # Generate output labels (underscore-joined channel names)
                    output_labels = [Symbol(join(group, "_")) for group in channel_groups]
                    
                    # Apply channel combination
                    channel_average!(item, 
                                   channel_selections = channel_selections,
                                   output_labels = output_labels,
                                   reduce = reduce)
                end
                
                # Save data with combined channels
                save(output_path, var_name, data)
                
                @info "  Saved: $output_path"
                processed_count += 1
                
            catch e
                @error "Error processing $file: $e"
                error_count += 1
                continue
            end
        end
        
        @info "Channel combining complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "combine_channels_data.log"
        log_dest = joinpath(output_dir, "combine_channels_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end
