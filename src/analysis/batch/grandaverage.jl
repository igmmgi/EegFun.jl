"""
    grandaverage(file_pattern::String; 
                 input_dir::String = pwd(), 
                 participants::Union{Int, Vector{Int}, Nothing} = nothing,
                 conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                 output_dir::Union{String, Nothing} = nothing)

Batch process ERP data files to create grand averages across participants.

This function loads JLD2 files containing ERP data from multiple participants, groups by condition,
and creates grand averages by averaging the EEG channel data across participants for each condition.

# Arguments
- `file_pattern::String`: Pattern to match JLD2 files (e.g., "erps_cleaned", "erps_original")
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant numbers to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition numbers to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: auto-generated)

# Returns
- `Nothing`: Saves grand averaged ERP data to output directory

# Examples
```julia
# Grand average all ERP files in current directory
grandaverage("erps_cleaned")

# Process specific participants and conditions
grandaverage_participants_data("erps_cleaned", 
                              input_dir = "/path/to/data", 
                              participants = [1, 2, 3], 
                              conditions = [1, 2])

# Specify custom output directory
grandaverage_participants_data("erps_cleaned", 
                              input_dir = "/path/to/data", 
                              output_dir = "/path/to/output")
```
"""
function grandaverage(file_pattern::String; 
                      input_dir::String = pwd(), 
                      participants::Union{Int, Vector{Int}, Nothing} = nothing,
                      conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                      output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "grandaverage.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch grand averaging started at $(now())"
        
        # Log the function call
        @log_call "grandaverage" (file_pattern,)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        
        # Validate inputs
        if !isdir(input_dir)
            @minimal_error_throw("Input directory does not exist: $input_dir")
        end
        
        # Validate that file pattern contains "erps"
        if !contains(file_pattern, "erps")
            @minimal_error_throw("grandaverage_participants_data only works with ERP data. File pattern must contain 'erps', got: '$file_pattern'")
        end
        
        # Create default output directory if not specified
        if output_dir === nothing
            output_dir = joinpath(input_dir, "grandaverage_$(file_pattern)")
        end
        
        # Create output directory if it doesn't exist
        mkpath(output_dir)
        
        # Find JLD2 files matching the pattern
        all_files = readdir(input_dir)
        jld2_files = Base.filter(x -> endswith(x, ".jld2") && contains(x, file_pattern), all_files)
        
        # Filter by participant number if specified
        if participants !== nothing
            jld2_files = _filter_files(jld2_files; include = participants)
        end
        
        if isempty(jld2_files)
            @minimal_warning "No JLD2 files found matching pattern '$file_pattern' in $input_dir"
            @info "Available files: $(Base.filter(x -> endswith(x, ".jld2"), all_files))"
            return nothing
        end
        
        @info "Found $(length(jld2_files)) JLD2 files matching pattern '$file_pattern'"
        
        # Load all ERP data and group by condition
        all_erps_by_condition = Dict{Int, Vector{ErpData}}()
        participant_info = Dict{Int, String}()  # Track which participant each ERP came from
        
        for (i, file) in enumerate(jld2_files)
            input_path = joinpath(input_dir, file)
            
            @info("Loading: $file ($(i)/$(length(jld2_files)))")
            
            try
                # Load the file and determine what data we have
                file_data = load(input_path)
                
                # Check what type of data we have
                if haskey(file_data, "erps")
                    erps_data = file_data["erps"]
                else
                    @minimal_warning "No 'erps' variable found in $file. Available: $(keys(file_data))"
                    continue
                end
                
                # Filter by condition if specified
                if conditions !== nothing
                    condition_nums = conditions isa Int ? [conditions] : conditions
                    # Filter data to only include specified conditions
                    erps_data = erps_data[condition_nums]
                    @info "  Filtering conditions: $condition_nums"
                end
                
                # Group ERPs by condition
                for erp in erps_data
                    # For ErpData, condition is stored in the DataFrame
                    cond_num = hasproperty(erp.data, :condition) ? erp.data[1, :condition] : 1
                    if !haskey(all_erps_by_condition, cond_num)
                        all_erps_by_condition[cond_num] = ErpData[]
                    end
                    push!(all_erps_by_condition[cond_num], erp)
                    participant_info[length(all_erps_by_condition[cond_num])] = file  # Track participant
                end
                
            catch e
                @error "Error loading $file: $e"
                continue
            end
        end
        
        if isempty(all_erps_by_condition)
            @minimal_warning "No valid ERP data found in any files"
            return nothing
        end
        
        @info "Found conditions: $(sort(collect(keys(all_erps_by_condition))))"
        
        # Create grand averages for each condition
        grand_averages = ErpData[]
        
        for (cond_num, erps) in all_erps_by_condition
            @info "Creating grand average for condition $cond_num (n=$(length(erps)) participants)"
            
            if length(erps) < 2
                @minimal_warning "Only $(length(erps)) participant(s) for condition $cond_num. Skipping grand average."
                continue
            end
            
            # Create grand average
            grand_avg = _create_grand_average(erps, cond_num)
            push!(grand_averages, grand_avg)
            
            @info "  Created grand average with $(nrow(grand_avg.data)) time points"
        end
        
        if isempty(grand_averages)
            @minimal_warning "No grand averages created (insufficient participants for any condition)"
            return nothing
        end
        
        # Save grand averages
        output_file = "grandaverage_$(file_pattern).jld2"
        output_path = joinpath(output_dir, output_file)
        save(output_path, "grand_averages", grand_averages)
        
        @info "Grand averaging complete! Created $(length(grand_averages)) grand averages, saved to: $output_path"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "grandaverage.log"
        log_dest = joinpath(output_dir, "grandaverage.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end

"""
    _create_grand_average(erps::Vector{ErpData}, cond_num::Int) -> ErpData

Create a grand average by averaging ERP data across participants for a specific condition.
"""
function _create_grand_average(erps::Vector{ErpData}, cond_num::Int)
    if isempty(erps)
        @minimal_error_throw("Cannot create grand average from empty ERP list")
    end
    
    # Validate that all ERPs have the same structure
    first_erp = erps[1]
    sample_rate = first_erp.sample_rate
    n_timepoints = nrow(first_erp.data)
    
    for (i, erp) in enumerate(erps[2:end])
        if erp.sample_rate != sample_rate
            @minimal_error_throw("ERP $i has different sample rate: $(erp.sample_rate) vs $(sample_rate)")
        end
        if nrow(erp.data) != n_timepoints
            @minimal_error_throw("ERP $i has different number of time points: $(nrow(erp.data)) vs $(n_timepoints)")
        end
    end
    
    # Get metadata columns and EEG channels
    metadata_cols = meta_labels(first_erp)
    eeg_channels = setdiff(propertynames(first_erp.data), metadata_cols)
    
    # Create a copy of the first ERP's data as the base
    grand_avg_data = copy(first_erp.data)
    
    # Update condition information
    grand_avg_data.condition .= cond_num
    grand_avg_data.condition_name .= "grand_average_condition_$(cond_num)"
    
    # Average EEG channels across participants
    for ch in eeg_channels
        # Collect data from all participants for this channel
        # Stack as columns: n_timepoints x n_participants
        channel_matrix = hcat([erp.data[!, ch] for erp in erps]...)
        
        # Average across participants (mean of each time point)
        grand_avg_data[!, ch] = vec(mean(channel_matrix, dims=2))
    end
    
    # Calculate total number of epochs across all participants
    total_epochs = sum(erp.n_epochs for erp in erps)
    
    # Use the layout and analysis info from the first ERP
    return ErpData(grand_avg_data, first_erp.layout, first_erp.sample_rate, first_erp.analysis_info, total_epochs)
end
