# configuration parameter
"""
    ConfigParameter{T}

A struct to define configuration parameters with type and validation constraints.

# Fields
- `type::Type{T}`: The Julia type of the parameter
- `min::Union{Nothing,T}`: Minimum value for numeric parameters (or nothing)
- `max::Union{Nothing,T}`: Maximum value for numeric parameters (or nothing)
- `description::String`: Human-readable description of the parameter
- `values::Union{Nothing,Vector{T}}`: List of allowed values (or nothing if any value is allowed)
- `default::Union{Nothing,T}`: Default value (or nothing if required)
"""
struct ConfigParameter{T}
    type::Type{T}
    min::Union{Nothing,T}
    max::Union{Nothing,T}
    description::String
    values::Union{Nothing,Vector{T}}
    default::Union{Nothing,T}
end

# Store parameter definitions as a constant
const PARAMETERS = Dict{String,ConfigParameter}(

    # File paths and settings
    "files.input.raw_data_directory" =>
        ConfigParameter{String}(String, nothing, nothing, "Directory containing raw data files", nothing, "."),
    "files.input.raw_data_files" => ConfigParameter{Union{Vector{String},String}}(
        Union{Vector{String},String},
        nothing,
        nothing,
        "Pattern for raw data files to process",
        nothing,
        "\\.bdf",
    ),
    "files.input.layout_file" => ConfigParameter{String}(
        String,
        nothing,
        nothing,
        "Electrode layout file name (\"*.csv\")",
        nothing,
        "biosemi72.csv",
    ),
    "files.output.output_data_directory" => ConfigParameter{String}(
        String,
        nothing,
        nothing,
        "Directory for processed data output",
        nothing,
        "./preprocessed_files",
    ),
    "files.output.save_continuous_data" =>
        ConfigParameter{Bool}(Bool, nothing, nothing, "Whether to save continuous data", nothing, false),
    "files.output.save_ica_data" =>
        ConfigParameter{Bool}(Bool, nothing, nothing, "Whether to save ICA results", nothing, true),
    "files.output.save_epoch_data" =>
        ConfigParameter{Bool}(Bool, nothing, nothing, "Whether to save epoched data", nothing, true),
    "files.output.save_erp_data" =>
        ConfigParameter{Bool}(Bool, nothing, nothing, "Whether to save ERP data", nothing, true),
    "files.output.delete_current_files" => ConfigParameter{Bool}(
        Bool,
        nothing,
        nothing,
        "Whether to delete existing files in output directory",
        nothing,
        true,
    ),
    "files.output.exit_early" => ConfigParameter{Bool}(
        Bool,
        nothing,
        nothing,
        "Whether to exit early from preprocessing pipeline (i.e., quick epoching only)",
        nothing,
        false,
    ),

    # Preprocessing settings
    "preprocessing.reference_electrode" =>
        ConfigParameter{String}(String, nothing, nothing, "Electrode(s) to use as reference", nothing, "avg"),

    # Filtering settings
    "filtering.highpass.type" => ConfigParameter{String}(
        String,
        nothing,
        nothing,
        "Type of highpass filter",
        ["fir", "iir"],
        "fir",
    ),
    "filtering.highpass.cutoff" =>
        ConfigParameter{Float64}(Float64, 0.01, 20.0, "High-pass filter cutoff in Hz", nothing, 0.1),
    "filtering.highpass.order" => ConfigParameter{Int}(Int, 1, 8, "Filter order", nothing, 1),

    # Filtering settings
    "filtering.lowpass.type" => ConfigParameter{String}(
        String,
        nothing,
        nothing,
        "Type of lowpass filter",
        ["fir", "iir"],
        "fir",
    ),
    "filtering.lowpass.cutoff" =>
        ConfigParameter{Float64}(Float64, 5, 200, "High-pass filter cutoff in Hz", nothing, 40),
    "filtering.lowpass.order" => ConfigParameter{Int}(Int, 1, 8, "Filter order", nothing, 3),

    # ICA settings
    "ica.ica_run" => ConfigParameter{Bool}(Bool, nothing, nothing, "Whether to run ICA", nothing, false),

    # Performance settings
    "performance.use_threading" => ConfigParameter{Bool}(
        Bool,
        nothing,
        nothing,
        "Whether to use threading for parallel processing",
        nothing,
        false,
    ),
    "performance.thread_count" => ConfigParameter{Int}(
        Int,
        1,
        nothing,
        "Number of threads to use for parallel processing (set to 0 to use all available threads)",
        nothing,
        1,
    ),
)

# Used for some validation stuff
"""
    ValidationResult

A struct to track validation results for configuration parameters.

# Fields
- `success::Bool`: Whether the validation was successful
- `error::Union{Nothing,String}`: Error message if validation failed, nothing otherwise
- `path::Union{Nothing,String}`: Path to the parameter that failed validation, nothing if validation succeeded

# Constructors
- `ValidationResult(success::Bool)`: Create a successful validation result
- `ValidationResult(success::Bool, error::String, path::String)`: Create a validation result with error details
"""
struct ValidationResult
    success::Bool
    error::Union{Nothing,String}
    path::Union{Nothing,String}
    # Constructors
    ValidationResult(success::Bool) = new(success, nothing, nothing)
    ValidationResult(success::Bool, error::String, path::String) = new(success, error, path)
end


"""
    load_config(config_file::String; validate::Bool=true)

Load and merge configuration from a TOML file with defaults.

# Arguments
- `config_file::String`: Path to the configuration file
- `validate::Bool`: Whether to validate the configuration after loading

# Returns
- `Union{Dict,Nothing}`: The loaded configuration or nothing if loading failed
"""
function load_config(config_file::String)
    # Load default config
    default_config = TOML.parsefile(joinpath(@__DIR__, "default.toml"))
    
    # Check if user config exists
    if !isfile(config_file)
        @minimal_error "Configuration file not found: $config_file"
    end
    
    # Load user config
    user_config = Dict()
    @info "Loading config file: $config_file"
    try
        user_config = TOML.parsefile(config_file)
    catch e
        @minimal_error "Error parsing TOML file: $(e.msg)"
    end
    
    # Merge and validate
    config = _merge_configs(default_config, user_config)
    validation_result = _validate_config(config)
    if !validation_result.success
        @minimal_error validation_result.error
    end
    
    return _extract_values(config)
end


"""
    _merge_configs(default_config::Dict, user_config::Dict)

Merge user config onto defaults, maintaining simple value structure.

# Arguments
- `default_config::Dict`: The default configuration
- `user_config::Dict`: The user configuration to merge

# Returns
- `Dict`: The merged configuration
"""
function _merge_configs(default_config::Dict, user_config::Dict)
    result = deepcopy(default_config)

    function merge_nested!(target::Dict, source::Dict)
        for (key, value) in source
            if !haskey(target, key)
                # New key, add it directly
                target[key] = value
                continue
            end

            if isa(value, Dict) && isa(target[key], Dict)
                # Both are dictionaries, merge recursively
                merge_nested!(target[key], value)
            else
                # Simple value or non-dict, update directly
                target[key] = value
            end
        end
    end

    merge_nested!(result, user_config)
    return result
end

"""
    _extract_values(config::Dict)

Extract final values from config.
"""
function _extract_values(config::Dict)
    result = Dict()
    for (key, value) in config
        if isa(value, Dict) # nested 
            result[key] = _extract_values(value)
        else # simple value
            result[key] = value
        end
    end
    return result
end

"""
    _validate_config(config::Dict, path="")

Validate config values against their metadata definitions.

# Arguments
- `config::Dict`: The configuration to validate
- `path::String`: The current path in the configuration (for nested validation)

# Returns
- `ValidationResult`: Result of the validation
"""
function _validate_config(config::Dict, path = "")
    # Use a stack to track nested dictionaries to validate
    stack = [(config, path)]

    while !isempty(stack)
        current_config, current_path = pop!(stack)

        for (key, value) in current_config
            new_path = isempty(current_path) ? key : "$current_path.$key"

            if isa(value, Dict)
                # Add nested dictionary to stack
                push!(stack, (value, new_path))
            else
                # Check if we have metadata for this parameter
                if haskey(PARAMETERS, new_path)
                    result = _validate_parameter(value, PARAMETERS[new_path], new_path)
                    !result.success && return result
                end
            end
        end
    end

    return ValidationResult(true)
end

"""
    _validate_parameter(value, parameter_spec::ConfigParameter, parameter_name::String)

Helper function to validate a single parameter value against its specification.

# Arguments
- `value`: The value to validate
- `parameter_spec::ConfigParameter`: The parameter specification to validate against
- `parameter_name::String`: The full name of the parameter in the configuration (e.g., "filtering.highpass.cutoff")

# Returns
- `ValidationResult`: Result of the validation
"""
function _validate_parameter(value, parameter_spec::ConfigParameter, parameter_name::String)
    # For numeric types, allow conversion between numeric types
    if parameter_spec.type <: Number
        # Check if value is any numeric type
        isa(value, Number) ||
            return ValidationResult(false, "$parameter_name must be a number, got $(typeof(value))", parameter_name)

        # Convert to the target type
        try
            value = convert(parameter_spec.type, value)
        catch
            return ValidationResult(
                false,
                "$parameter_name must be convertible to $(parameter_spec.type), got $(typeof(value))",
                parameter_name,
            )
        end

        # Check min/max constraints
        !isnothing(parameter_spec.min) &&
            value < parameter_spec.min &&
            return ValidationResult(false, "$parameter_name ($value) must be >= $(parameter_spec.min)", parameter_name)

        !isnothing(parameter_spec.max) &&
            value > parameter_spec.max &&
            return ValidationResult(false, "$parameter_name ($value) must be <= $(parameter_spec.max)", parameter_name)
    else
        # For non-numeric types, require exact type match
        isa(value, parameter_spec.type) || return ValidationResult(
            false,
            "$parameter_name must be of type $(parameter_spec.type), got $(typeof(value))",
            parameter_name,
        )
    end

    # Check allowed values if they exist
    if !isnothing(parameter_spec.values)
        value_in_allowed = false
        
        # Convert value to string for comparison if needed
        for allowed_value in parameter_spec.values
            if isequal(value, allowed_value)
                value_in_allowed = true
                break
            end
        end
        
        value_in_allowed || return ValidationResult(
            false,
            "$parameter_name ($value) must be one of: $(join(parameter_spec.values, ", "))",
            parameter_name,
        )
    end

    return ValidationResult(true)
end


"""
    show_parameter_info(parameter_name::String="")

Display information about configuration parameters. If parameter_name is empty, shows all parameters.
If parameter_name is provided, shows detailed information about that specific parameter.

# Arguments
- `parameter_name::String`: Optional path to a specific parameter (e.g., "filtering.highpass.cutoff")

# Example
```julia
# Show all parameters
show_parameter_info()

# Show specific parameter
show_parameter_info("filtering.highpass.cutoff")
```
"""
function show_parameter_info(parameter_name::String = "")
    if isempty(parameter_name)
        # Show all parameters
        @info "Available Configuration Parameters:"
        @info "==================================="

        # Group parameters by section and subsection
        sections = Dict{String,Dict{String,Vector{Tuple{String,ConfigParameter}}}}()
        
        for (path, parameter_spec) in PARAMETERS
            parts = split(path, ".")
            section = parts[1]
            subsection = length(parts) > 1 ? parts[2] : ""
            
            if !haskey(sections, section)
                sections[section] = Dict{String,Vector{Tuple{String,ConfigParameter}}}()
            end
            if !haskey(sections[section], subsection)
                sections[section][subsection] = Tuple{String,ConfigParameter}[]
            end
            push!(sections[section][subsection], (path, parameter_spec))
        end

        # Sort sections
        sorted_sections = sort(collect(keys(sections)))

        # Display each section
        for section in sorted_sections
            @info "[$section]"
            @info "-"^(length(section) + 2)

            # Sort subsections
            sorted_subsections = sort(collect(keys(sections[section])))
            
            for subsection in sorted_subsections
                if !isempty(subsection)
                    @info "  [$subsection]"
                end
                
                # Sort parameters within subsection
                sorted_params = sort(sections[section][subsection], by=first)
                for (path, parameter_spec) in sorted_params
                    # Get the parameter name (last part of the path)
                    param_name = last(split(path, "."))
                    indent = isempty(subsection) ? "  " : "    "
                    @info "$indent$param_name: $(parameter_spec.description)"
                end
            end
        end

        @info "Use show_parameter_info(\"section.parameter\") for detailed information about a specific parameter."
    else
        # Show specific parameter
        if !haskey(PARAMETERS, parameter_name)
            @warn "Parameter not found: $parameter_name"
            return
        end

        parameter_spec = PARAMETERS[parameter_name]
        @info "Parameter: $parameter_name"
        @info "="^(length(parameter_name) + 11)
        @info "Description: $(parameter_spec.description)"
        @info "Type: $(parameter_spec.type)"

        if !isnothing(parameter_spec.min) || !isnothing(parameter_spec.max)
            range_str = "Range: "
            if !isnothing(parameter_spec.min)
                range_str *= "$(parameter_spec.min) ≤ "
            end
            range_str *= "value"
            if !isnothing(parameter_spec.max)
                range_str *= " ≤ $(parameter_spec.max)"
            end
            @info range_str
        end

        if !isnothing(parameter_spec.values)
            @info "Allowed values: $(join(parameter_spec.values, ", "))"
        end

        if !isnothing(parameter_spec.default)
            @info "Default: $(parameter_spec.default)"
        end
    end
end

"""
    generate_config_template(filename::String="config_template.toml")

Generate and save a template TOML configuration file with all available parameters and their default values.

# Arguments
- `filename::String`: Name of the template file to create (default: "config_template.toml")
"""
function generate_config_template(filename::String="config_template.toml")
    try
        @info "Starting config template generation"
        
        # Open file for writing
        open(filename, "w") do io
            # Write header
            println(io, "# EEG Processing Configuration Template")
            println(io, "# Generated on ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
            println(io)
            println(io, "# This template shows all available configuration options.")
            println(io, "# Required fields are marked with [REQUIRED]")
            println(io, "# Default values are shown where available")
            println(io)

            # Group parameters by section and subsection
            sections = _group_parameters_by_section()
            
            # Sort sections
            sorted_sections = sort(collect(keys(sections)))

            # Write each section
            for section in sorted_sections
                println(io, "\n# $section Settings")
                println(io, "[$section]")
                
                # Sort subsections
                sorted_subsections = sort(collect(keys(sections[section])))
                
                for subsection in sorted_subsections
                    if !isempty(subsection)
                        println(io, "\n# $subsection Settings")
                        println(io, "[$section.$subsection]")
                    end
                    
                    # Sort parameters within subsection
                    sorted_params = sort(sections[section][subsection], by=first)
                    for (path, parameter_spec) in sorted_params
                        # Get the parameter name (last part of the path)
                        param_name = last(split(path, "."))
                        
                        # Write parameter documentation
                        _write_parameter_docs(io, parameter_spec)
                        
                        # Use the default value from parameter_spec
                        value = parameter_spec.default
                        
                        # Write parameter value in appropriate format
                        _write_parameter_value(io, param_name, value)
                    end
                end
            end
        end
        @info "Configuration template saved to: $filename"
    catch e
        @error "Error in generate_config_template" exception=(e, catch_backtrace())
        if e isa MethodError
            @minimal_error "Failed to save configuration template: Method error occurred"
        else
            @minimal_error "Failed to save configuration template: $e"
        end
    end
end

"""
    _group_parameters_by_section()

Group configuration parameters by section and subsection.

# Returns
- `Dict`: A nested dictionary mapping sections and subsections to parameter specs
"""
function _group_parameters_by_section()
    sections = Dict{String,Dict{String,Vector{Tuple{String,ConfigParameter}}}}()
    
    for (path, parameter_spec) in PARAMETERS
        parts = split(path, ".")
        section = parts[1]
        subsection = length(parts) > 2 ? parts[2] : ""  # Only use subsection if there are more than 2 parts
        
        if !haskey(sections, section)
            sections[section] = Dict{String,Vector{Tuple{String,ConfigParameter}}}()
        end
        if !haskey(sections[section], subsection)
            sections[section][subsection] = Tuple{String,ConfigParameter}[]
        end
        push!(sections[section][subsection], (path, parameter_spec))
    end
    
    return sections
end

"""
    _write_parameter_docs(io::IO, parameter_spec::ConfigParameter)

Write parameter documentation to the given IO stream.
"""
function _write_parameter_docs(io::IO, parameter_spec::ConfigParameter)
    println(io, "\n# $(parameter_spec.description)")
    println(io, "# Type: $(parameter_spec.type)")
    
    if !isnothing(parameter_spec.min) || !isnothing(parameter_spec.max)
        range_str = "# Range: "
        if !isnothing(parameter_spec.min)
            range_str *= "$(parameter_spec.min) ≤ "
        end
        range_str *= "value"
        if !isnothing(parameter_spec.max)
            range_str *= " ≤ $(parameter_spec.max)"
        end
        println(io, range_str)
    end
    
    if !isnothing(parameter_spec.values)
        println(io, "# Allowed values: $(join(parameter_spec.values, ", "))")
    end
    
    # Print default value if it exists
    if !isnothing(parameter_spec.default)
        println(io, "# Default: $(parameter_spec.default)")
    end
    
    # Mark required fields
    if isnothing(parameter_spec.default)
        println(io, "# [REQUIRED]")
    end
end

"""
    _write_parameter_value(io::IO, param_name::Union{String,SubString{String}}, value)

Write a parameter value to the given IO stream in the appropriate TOML format.
"""
function _write_parameter_value(io::IO, param_name::Union{String,SubString{String}}, value)
    if isa(value, String)
        # Double any backslashes in the string for TOML
        escaped_value = replace(value, "\\" => "\\\\")
        println(io, "$(String(param_name)) = \"$escaped_value\"")
    elseif isa(value, Vector)
        if isempty(value)
            println(io, "$(String(param_name)) = []")
        else
            println(io, "$(String(param_name)) = [")
            for item in value
                if isa(item, String)
                    # Double any backslashes in the string for TOML
                    escaped_item = replace(item, "\\" => "\\\\")
                    println(io, "    \"$escaped_item\",")
                else
                    println(io, "    $item,")
                end
            end
            println(io, "]")
        end
    else
        println(io, "$(String(param_name)) = $value")
    end
end
