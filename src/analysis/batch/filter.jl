"""
    filter_data(file_pattern::String, cutoff_freq::Real; 
                input_dir::String = pwd(), 
                filter_type::String = "lp", 
                participants::Union{Int, Vector{Int}, Nothing} = nothing,
                conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                output_dir::Union{String, Nothing} = nothing)

Filter EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `cutoff_freq::Real`: Cutoff frequency in Hz
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `filter_type::String`: Type of filter ("lp", "hp", "bp", "bs")
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on filter settings)

# Example
```julia
# Filter all epochs at 30 Hz
filter_data("epochs", 30.0)

# Filter specific participant
filter_data("epochs", 30.0, participants=3)

# Filter specific participants and conditions
filter_data("epochs", 30.0, participants=[3, 4], conditions=[1, 2])

# Filter specific directory with participant and condition
filter_data("epochs", 30.0, input_dir="/path/to/input/", participants=3, conditions=1)
```
"""
function filter_data(file_pattern::String, cutoff_freq::Real; 
                    input_dir::String = pwd(), 
                    filter_type::String = "lp", 
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "filter_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch filtering started at $(now())"
        
        # Log the function call generically
        args = [file_pattern, cutoff_freq]
        kwargs = [:input_dir => input_dir, :filter_type => filter_type, :participants => participants, :conditions => conditions, :output_dir => output_dir]
        _log_function_call("filter_data", args, kwargs)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        @info "Filter type: $filter_type, Cutoff frequency: $cutoff_freq Hz"
        
        # Validate inputs
    if !isdir(input_dir)
        @minimal_error_throw("Input directory does not exist: $input_dir")
    end
    
    if cutoff_freq <= 0
        @minimal_error_throw("Cutoff frequency must be positive, got: $cutoff_freq")
    end
    
    if filter_type ∉ ["lp", "hp"]
        @minimal_error_throw("Filter type must be one of: lp, hp, got: $filter_type")
    end
    
    # Create default output directory if not specified
    if output_dir === nothing
        output_dir = joinpath(input_dir, "filtered_$(file_pattern)_$(filter_type)_$(cutoff_freq)hz")
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
    @info "Filter settings: $filter_type filter, cutoff: $cutoff_freq Hz"
    
    processed_count = 0
    error_count = 0
    
    for (i, file) in enumerate(jld2_files)
        input_path = joinpath(input_dir, file)
        output_path = joinpath(output_dir, file)
        
        @info("Processing: $file ($(i)/$(length(jld2_files)))")
        
        try
            # Define possible variable names to try
            possible_vars = ["erps", "epochs"]
            file_data = load(input_path)

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
            
            # Apply filter to each item
            for (i, item) in enumerate(data)
                @info "  Filtering $var_name $i/$(length(data))"
                filter_data!(item, filter_type, cutoff_freq)
            end
            
            # Save filtered data
            save(output_path, var_name, data)
            
            @info "  Saved: $output_path"
            processed_count += 1
            
        catch e
            @error "Error processing $file: $e"
            error_count += 1
            continue
        end
    end
    
        @info "Filtering complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "filter_data.log"
        log_dest = joinpath(output_dir, "filter_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end
