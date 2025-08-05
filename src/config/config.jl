# configuration parameter
"""
    ConfigParameter{T}

A struct to define configuration parameters with type and validation constraints.

# Fields
- `description::String`: Human-readable description of the parameter
- `min::Union{Nothing,T}`: Minimum value for numeric parameters (or nothing)
- `max::Union{Nothing,T}`: Maximum value for numeric parameters (or nothing)
- `allowed_values::Union{Nothing,Vector{T}}`: List of allowed values (or nothing if any value is allowed)
- `default::Union{Nothing,T}`: Default value (or nothing if required)
"""
struct ConfigParameter{T}
    description::String
    min::Union{Nothing,T}
    max::Union{Nothing,T}
    allowed_values::Union{Nothing,Vector{T}}
    default::Union{Nothing,T}
end

# Keyword constructor with defaults
function ConfigParameter{T}(;
    description::String,
    default::Union{Nothing,T} = nothing,
    min::Union{Nothing,T} = nothing,
    max::Union{Nothing,T} = nothing,
    allowed_values::Union{Nothing,Vector{T}} = nothing,
) where {T}
    ConfigParameter{T}(description, min, max, allowed_values, default)
end


# Store parameter definitions as a constant
const PARAMETERS = Dict{String,ConfigParameter}(

    # File paths and settings
    "files.input.directory" =>
        ConfigParameter{String}(description = "Directory containing raw data files.", default = "."),
    "files.input.raw_data_files" => ConfigParameter{Union{Vector{String},String}}(
        description = "Pattern (regex or explicit list) for raw data files to process.",
        default = "\\.bdf",
    ),
    "files.input.layout_file" =>
        ConfigParameter{String}(description = "Electrode layout file name (\"*.csv\")", default = "biosemi72.csv"),
    "files.input.epoch_condition_file" => ConfigParameter{Union{Nothing,String}}(
        description = "TOML file that defines the condition epochs.",
        default = "",
    ),
    "files.output.directory" => ConfigParameter{String}(
        description = "Directory for processed output files",
        default = "./preprocessed_files",
    ),

    # What data should we save?
    "files.output.save_continuous_data" =>
        ConfigParameter{Bool}(description = "Save continuous data?", default = false),
    "files.output.save_ica_data" => ConfigParameter{Bool}(description = "Save ICA results?", default = true),
    "files.output.save_epoch_data_original" =>
        ConfigParameter{Bool}(description = "Save epoched data?", default = true),
    "files.output.save_epoch_data_cleaned" =>
        ConfigParameter{Bool}(description = "Save epoched data after cleaning?", default = true),
    "files.output.save_erp_data_original" => ConfigParameter{Bool}(description = "Save ERP data", default = true),
    "files.output.save_erp_data_cleaned" =>
        ConfigParameter{Bool}(description = "Save ERP data after cleaning?", default = true),
    "files.output.exit_early" => ConfigParameter{Bool}(
        description = "Exit early from preprocessing pipeline (i.e., quick epoching only)",
        default = false,
    ),

    # Preprocessing settings
    "preprocess.epoch_start" => ConfigParameter{Real}(description = "Epoch start (seconds).", default = -1),
    "preprocess.epoch_end" => ConfigParameter{Real}(description = "Epoch end (seconds).", default = 1),
    "preprocess.reference_channel" =>
        ConfigParameter{String}(description = "Channels(s) to use as reference", default = "avg"),
    "preprocess.layout.neighbour_criterion" => ConfigParameter{Real}(
        description = "Distance criterion (in mm) for channel neighbour definition.",
        default = 40,
        min = 0,
    ),
    "preprocess.eog.vEOG_channels" => ConfigParameter{Vector{Vector{String}}}(
        description = "Channels used in the calculation of vertical eye movements (vEOG).",
        default = [["Fp1", "IO1"], ["Fp2", "IO2"], ["vEOG"]],
    ),
    "preprocess.eog.hEOG_channels" => ConfigParameter{Vector{Vector{String}}}(
        description = "Channels used in the calculation of horizontal eye movements (hEOG).",
        default = [["F9"], ["F10"], ["hEOG"]],
    ),
    "preprocess.eog.vEOG_criterion" => ConfigParameter{Real}(
        description = "Distance criterion for vertical EOG channel definition.",
        default = 50,
        min = 0,
    ),
    "preprocess.eog.hEOG_criterion" => ConfigParameter{Real}(
        description = "Distance criterion for horizontal EOG channel definition.",
        default = 30,
        min = 0,
    ),
    "preprocess.eeg.extreme_value_criterion" => ConfigParameter{Real}(
        description = "Value (mV) for defining data section as an extreme value.",
        default = 500,
    ),
    "preprocess.eeg.artifact_value_criterion" => ConfigParameter{Real}(
        description = "Value (mV) for defining data section as an artifact value.",
        default = 100,
    ),

    # Filtering settings
    "filter.highpass.on" => ConfigParameter{Bool}(description = "Apply highpass filter true/false", default = true),
    "filter.highpass.method" => ConfigParameter{String}(
        description = "Type of filter",
        default = "fir",
        allowed_values = ["fir", "iir"],
    ),
    "filter.highpass.filter_func" => ConfigParameter{Function}(
        description = "Filter function",
        default = filtfilt,
        allowed_values = [filt, filtfilt],
    ),
    "filter.highpass.cutoff_freq" => ConfigParameter{Real}(
        description = "High-pass filter cutoff frequency (Hz)",
        default = 0.1,
        min = 0.01,
        max = 20.0,
    ),
    "filter.highpass.order" => ConfigParameter{Int}(description = "Filter order", default = 1, min = 1, max = 4),

    # Lowpass filtering settings
    "filter.lowpass.on" => ConfigParameter{Bool}(description = "Apply lowpass filter true/false", default = true),
    "filter.lowpass.method" => ConfigParameter{String}(
        description = "Type of filter",
        default = "fir",
        allowed_values = ["fir", "iir"],
    ),
    "filter.lowpass.filter_func" => ConfigParameter{Function}(
        description = "Filter function",
        default = filtfilt,
        allowed_values = [filt, filtfilt],
    ),
    "filter.lowpass.cutoff_freq" => ConfigParameter{Real}(
        description = "Low-pass filter cutoff frequency (Hz)",
        default = 30,
        min = 5,
        max = 500,
    ),
    "filter.lowpass.order" => ConfigParameter{Int}(description = "Filter order", default = 3, min = 1, max = 8),

    # ICA settings
    "ica.run" => ConfigParameter{Bool}(
        description = "Run Independent Component Analysis (ICA) true/false.",
        default = false,
    ),

    "filter.ica_highpass.on" => ConfigParameter{Bool}(description = "Apply highpass filter ICA data true/false", default = true),
    "filter.ica_highpass.method" => ConfigParameter{String}(
        description = "Type of filter",
        default = "fir",
        allowed_values = ["fir", "iir"],
    ),
    "filter.ica_highpass.filter_func" => ConfigParameter{Function}(
        description = "Filter function",
        default = filtfilt,
        allowed_values = [filt, filtfilt],
    ),
    "filter.ica_highpass.cutoff_freq" => ConfigParameter{Real}(
        description = "High-pass filter cutoff frequency (Hz)",
        default = 1,
        min = 1,
        max = 20.0,
    ),
    "filter.ica_highpass.order" => ConfigParameter{Int}(description = "Filter order", default = 1, min = 1, max = 4),

    "filter.ica_lowpass.on" => ConfigParameter{Bool}(description = "Apply lowpass filter ICA data true/false", default = true),
    "filter.ica_lowpass.method" => ConfigParameter{String}(
        description = "Type of filter",
        default = "fir",
        allowed_values = ["fir", "iir"],
    ),
    "filter.ica_lowpass.filter_func" => ConfigParameter{Function}(
        description = "Filter function",
        default = filtfilt,
        allowed_values = [filt, filtfilt],
    ),
    "filter.ica_lowpass.cutoff_freq" => ConfigParameter{Real}(
        description = "Low-pass filter cutoff frequency (Hz)",
        default = 30,
        min = 5,
        max = 500,
    ),
    "filter.ica_highpass.order" => ConfigParameter{Int}(description = "Filter order", default = 3, min = 1, max = 4),
   
)

"""
    ValidationResult

Result of parameter validation.

# Fields
- `success::Bool`: Whether validation succeeded
- `error::Union{Nothing,String}`: Error message if validation failed, nothing if validation succeeded (default: nothing)
- `path::Union{Nothing,String}`: Path to the parameter that failed validation, nothing if validation succeeded (default: nothing)
"""
@kwdef struct ValidationResult
    success::Bool
    error::Union{Nothing,String} = nothing
    path::Union{Nothing,String} = nothing
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
        return nothing
    end

    # Load user config
    user_config = Dict()
    @info "Loading config file: $config_file"
    try
        user_config = TOML.parsefile(config_file)
    catch e
        @minimal_error "Error parsing TOML file: $e"
        return nothing
    end

    # Merge and validate
    config = _merge_configs(default_config, user_config)
    validation_result = _validate_config(config)
    if !validation_result.success
        @minimal_error validation_result.error
        return nothing
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
    result = copy(default_config)

    function merge_nested!(target::Dict, source::Dict)
        for (key, value) in source
            if !haskey(target, key)
                target[key] = value
                continue
            end
            if isa(value, Dict) && isa(target[key], Dict)
                merge_nested!(target[key], value)
            else
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

    return ValidationResult(success = true)
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
    # Get the type from the type parameter
    param_type = typeof(parameter_spec).parameters[1]

    # Check if value is the right type
    if param_type <: Number
        # For numeric types, accept any number
        isa(value, Number) ||
            return ValidationResult(success = false, error = "$parameter_name must be a number, got $(typeof(value))", path = parameter_name)
    else
        # For non-numeric types, require exact type match
        isa(value, param_type) ||
            return ValidationResult(success = false, error = "$parameter_name must be of type $param_type, got $(typeof(value))", path = parameter_name)
    end

    # Check min/max constraints
    if !isnothing(parameter_spec.min) && value < parameter_spec.min
        return ValidationResult(success = false, error = "$parameter_name ($value) must be >= $(parameter_spec.min)", path = parameter_name)
    end
    if !isnothing(parameter_spec.max) && value > parameter_spec.max
        return ValidationResult(success = false, error = "$parameter_name ($value) must be <= $(parameter_spec.max)", path = parameter_name)
    end

    # Check allowed values if they exist
    if !isnothing(parameter_spec.allowed_values) && !(value in parameter_spec.allowed_values)
        return ValidationResult(
            success = false,
            error = "$parameter_name ($value) must be one of: $(join(parameter_spec.allowed_values, ", "))",
            path = parameter_name
        )
    end

    return ValidationResult(success = true)
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

Show information about configuration parameters.

# Arguments
- `parameter_name::String`: Parameter name or section name (default: "" for overview)

# Examples
```julia
# Show overview of all parameters
show_parameter_info()

# Show all parameters in a section
show_parameter_info("preprocess")

# Show all parameters in a subsection
show_parameter_info("files.input")

# Show specific parameter details
show_parameter_info("preprocess.epoch_start")
```
"""
function show_parameter_info(parameter_name::String = "")
    if isempty(parameter_name)
        # Show all parameters
        @info "Available Configuration Parameters:"
        @info "==================================="

        # Group parameters by section and subsection
        sections = _group_parameters_by_section()

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
                sorted_params = sort(sections[section][subsection], by = first)
                for (path, parameter_spec) in sorted_params
                    # Get the parameter name (last part of the path)
                    param_name = last(split(path, "."))
                    indent = isempty(subsection) ? "  " : "    "
                    @info "$indent$param_name: $(parameter_spec.description)"
                end
            end
        end

        @info "Use show_parameter_info(\"section\") for section overview"
        @info "Use show_parameter_info(\"section.parameter\") for specific parameter details"
    else
        # Check if it's an exact parameter match first
        if haskey(PARAMETERS, parameter_name)
            _show_parameter_details(parameter_name)
        else
            # Check if it's a section or partial match
            matching_params = collect(filter(keys(PARAMETERS)) do key
                startswith(key, parameter_name)
            end)
            
            if !isempty(matching_params)
                _show_section_overview(parameter_name, matching_params)
            else
                @warn "Parameter or section not found: $parameter_name"
                @info "Use show_parameter_info() to see all available parameters and sections"
            end
        end
    end
end

"""
    _show_parameter_details(parameter_name::String)

Show detailed information about a specific parameter.
"""
function _show_parameter_details(parameter_name::String)
    parameter_spec = PARAMETERS[parameter_name]
    @info "Parameter: $parameter_name"
    @info "="^(length(parameter_name) + 11)
    @info "Description: $(parameter_spec.description)"
    @info "Type: $(typeof(parameter_spec).parameters[1])"

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

    if !isnothing(parameter_spec.allowed_values)
        @info "Allowed values: $(join(parameter_spec.allowed_values, ", "))"
    end

    if !isnothing(parameter_spec.default)
        @info "Default: $(parameter_spec.default)"
    end
end

"""
    _show_section_overview(section_name::String, matching_params::Vector{String})

Show overview of all parameters in a section.
"""
function _show_section_overview(section_name::String, matching_params::Vector{String})
    @info "Section: $section_name"
    @info "="^(length(section_name) + 9)
    
    # Group parameters by subsection
    sections = Dict{String,Vector{Tuple{String,ConfigParameter}}}()
    
    for param_path in matching_params
        parts = split(param_path, ".")
        if length(parts) > 1
            # For nested paths like "files.input.directory", 
            # subsection is everything after the section name
            section_prefix = section_name * "."
            if startswith(param_path, section_prefix)
                subsection_path = param_path[length(section_prefix)+1:end]
                subsection_parts = split(subsection_path, ".")
                if length(subsection_parts) > 1
                    subsection = join(subsection_parts[1:end-1], ".")
                else
                    subsection = ""
                end
            else
                subsection = ""
            end
        else
            subsection = ""
        end
        
        if !haskey(sections, subsection)
            sections[subsection] = Tuple{String,ConfigParameter}[]
        end
        push!(sections[subsection], (param_path, PARAMETERS[param_path]))
    end
    
    # Sort and display subsections
    sorted_subsections = sort(collect(keys(sections)))
    
    for subsection in sorted_subsections
        if !isempty(subsection)
            @info "  [$subsection]"
        end
        
        # Sort parameters within subsection
        sorted_params = sort(sections[subsection], by = first)
        for (path, parameter_spec) in sorted_params
            # Get the parameter name (last part of the path)
            param_name = last(split(path, "."))
            indent = isempty(subsection) ? "  " : "    "
            @info "$indent$param_name: $(parameter_spec.description)"
        end
    end
    
    @info ""
    @info "Use show_parameter_info(\"$section_name.parameter_name\") for detailed information about a specific parameter"
end

"""
    generate_config_template(filename::String="config_template.toml")

Generate and save a template TOML configuration file with all available parameters and their default values.

# Arguments
- `filename::String`: Name of the template file to create (default: "config_template.toml")
"""
function generate_config_template(filename::String = "config_template.toml")
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
                    sorted_params = sort(sections[section][subsection], by = first)
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

        # For nested paths like "ica.ica_filter.highpass.on", 
        # we want subsection to be "ica_filter.highpass"
        if length(parts) > 2
            subsection = join(parts[2:(end-1)], ".")  # Join all parts except first and last
        else
            subsection = ""
        end

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
    println(io, "# Type: $(typeof(parameter_spec).parameters[1])")

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

    if !isnothing(parameter_spec.allowed_values)
        println(io, "# Allowed values: $(join(parameter_spec.allowed_values, ", "))")
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
