"""
    rereference_data(file_pattern::String; 
                    input_dir::String = pwd(), 
                    reference_type::String = "average", 
                    participants::Union{Int, Vector{Int}, Nothing} = nothing,
                    conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                    output_dir::Union{String, Nothing} = nothing)

Rereference EEG/ERP data from JLD2 files and save to a new directory.

# Arguments
- `file_pattern::String`: Pattern to match files ("epochs", "erps", "cleaned", "original", or custom)
- `input_dir::String`: Input directory containing JLD2 files (default: current directory)
- `reference_type::String`: Type of reference ("average", "mastoid", "cz", or custom channel name)
- `participants::Union{Int, Vector{Int}, Nothing}`: Participant number(s) to process (default: all)
- `conditions::Union{Int, Vector{Int}, Nothing}`: Condition number(s) to process (default: all)
- `output_dir::Union{String, Nothing}`: Output directory (default: creates subdirectory based on reference settings)

# Example
```julia
# Rereference all epochs to average reference
rereference_data("epochs")

# Rereference specific participant to mastoid reference
rereference_data("epochs", reference_type="mastoid", participants=3)

# Rereference specific participants and conditions
rereference_data("epochs", reference_type="cz", participants=[3, 4], conditions=[1, 2])

# Rereference specific directory with participant and condition
rereference_data("epochs", input_dir="/path/to/input/", participants=3, conditions=1)
```
"""
function rereference_data(file_pattern::String; 
                         input_dir::String = pwd(), 
                         reference_type::String = "average", 
                         participants::Union{Int, Vector{Int}, Nothing} = nothing,
                         conditions::Union{Int, Vector{Int}, Nothing} = nothing,
                         output_dir::Union{String, Nothing} = nothing)
    
    # Set up global logging
    log_file = "rereference_data.log"
    setup_global_logging(log_file)
    
    try
        @info "Batch rereferencing started at $(now())"
        
        # Log the function call
        @log_call "rereference_data" (file_pattern,)
        
        @info "File pattern: $file_pattern"
        @info "Input directory: $input_dir"
        @info "Reference type: $reference_type"
        
        # Validate inputs
    if !isdir(input_dir)
        @minimal_error_throw("Input directory does not exist: $input_dir")
    end
    
    if reference_type ∉ [:avg, :mastoid]
        @minimal_warning "Reference type '$reference_type' not in common types (avg, mastoid). Proceeding with custom reference."
    end
    
    # Create default output directory if not specified
    if output_dir === nothing
        output_dir = joinpath(input_dir, "rereferenced_$(file_pattern)_$(reference_type)")
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
    @info "Reference settings: $reference_type reference"
    
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
            
            # Apply rereferencing to each item
            for (i, item) in enumerate(data)
                @info "  Rereferencing $var_name $i/$(length(data))"
                rereference!(item, reference_type)
            end
            
            # Save rereferenced data
            save(output_path, var_name, data)
            
            @info "  Saved: $output_path"
            processed_count += 1
            
        catch e
            @error "Error processing $file: $e"
            error_count += 1
            continue
        end
    end
    
        @info "Rereferencing complete! Processed $processed_count files successfully, $error_count errors, output saved to: $output_dir"
        
    finally
        # Close global logging and move log file to output directory
        close_global_logging()
        log_source = "rereference_data.log"
        log_dest = joinpath(output_dir, "rereference_data.log")
        if log_source != log_dest
            mv(log_source, log_dest, force = true)
        end
    end
end
