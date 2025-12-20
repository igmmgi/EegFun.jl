# =============================================================================
# TYPES AND STRUCTURES
# =============================================================================

"""
    ConfigParameter{T}

A struct to define configuration parameters with type and validation constraints.

# Fields
- `description::String`: Human-readable description of the parameter
- `default::Union{Nothing,T}`: Default value (or nothing if required)
- `allowed::Union{Nothing,Vector{String}}`: List of allowed values (or nothing if any value is allowed)
- `min::Union{Nothing,T}`: Minimum value for numeric parameters (or nothing)
- `max::Union{Nothing,T}`: Maximum value for numeric parameters (or nothing)
"""
struct ConfigParameter{T}
    description::String
    default::Union{Nothing,T}
    allowed::Union{Nothing,Vector{String}}
    min::Union{Nothing,T}
    max::Union{Nothing,T}
end

# Keyword constructor with defaults + helper for parameter definitions
function ConfigParameter{T}(;
    description::String,
    default::Union{Nothing,T} = nothing,
    allowed::Union{Nothing,Vector{String}} = nothing,
    min::Union{Nothing,T} = nothing,
    max::Union{Nothing,T} = nothing,
) where {T}
    ConfigParameter{T}(description, default, allowed, min, max)
end

# =============================================================================
# PARAMETER CONSTRUCTOR HELPERS
# =============================================================================

# Helper to create ConfigParameter with common defaults
function _param(::Type{T}, desc, default = nothing; allowed = nothing, min = nothing, max = nothing) where {T}
    ConfigParameter{T}(description = desc, default = default, allowed = allowed, min = min, max = max)
end

string_param(desc, default = ""; allowed = nothing) =
    _param(Union{Vector{String},String}, desc, default, allowed = allowed)
simple_string_param(desc, default = ""; allowed = nothing) = _param(String, desc, default, allowed = allowed)
bool_param(desc, default = false) = _param(Bool, desc, default)
number_param(desc, default, min = nothing, max = nothing) = _param(Real, desc, default, min = min, max = max)
channel_groups_param(desc, default) = _param(Vector{Vector{String}}, desc, default)

# Helper function to create filter parameter specifications
function _filter_param_spec(prefix, apply, type, freq, min_freq, max_freq, order, min_order, max_order)
    Dict(
        "$prefix.apply"  => bool_param("Apply: true/false", apply),
        "$prefix.type"   => string_param("Filter type identifier", type, allowed = ["hp", "lp"]),
        "$prefix.method" => string_param("Filter type", "iir", allowed = ["fir", "iir"]),
        "$prefix.func"   => string_param("Filter function", "filtfilt", allowed = ["filt", "filtfilt"]),
        "$prefix.freq"   => number_param("Cutoff frequency (Hz)", freq, min_freq, max_freq),
        "$prefix.order"  => number_param("Filter order", order, min_order, max_order),
    )
end


# =============================================================================
# PARAMETER DEFINITIONS
# =============================================================================

# fmt: off
const PARAMETERS = Dict{String,ConfigParameter}(

    # File paths and settings
    "files.input.directory"            => simple_string_param("Directory containing raw data files.", "."),
    "files.input.raw_data_files"       => string_param("Pattern (regex or explicit list) for raw data files to process.", "\\.bdf"),
    "files.input.layout_file"          => simple_string_param("Electrode layout file name (\"*.csv\")", "biosemi72.csv"),
    "files.input.epoch_condition_file" => simple_string_param("TOML file that defines the condition epochs.", ""),
    "files.output.directory"           => simple_string_param("Directory for processed output files", "./preprocessed_files"),

    # What data should we save?
    "files.output.save_continuous_data_original" => bool_param("Save continuous data original?", true),
    "files.output.save_continuous_data_cleaned"  => bool_param("Save continuous data cleaned?", true),
    "files.output.save_ica_data"                 => bool_param("Save ICA results?", true),
    "files.output.save_epoch_data_original"      => bool_param("Save epoched data original?", true),
    "files.output.save_epoch_data_cleaned"       => bool_param("Save epoched data cleaned?", true),
    "files.output.save_epoch_data_good"          => bool_param("Save epoched data good?", true),
    "files.output.save_erp_data_original"        => bool_param("Save ERP data original?", true),
    "files.output.save_erp_data_cleaned"         => bool_param("Save ERP data cleaned?", true),
    "files.output.save_erp_data_good"            => bool_param("Save ERP data good?", true),

    # Preprocessing settings
    "preprocess.epoch_start"                      => number_param("Epoch start (seconds).", -1),
    "preprocess.epoch_end"                        => number_param("Epoch end (seconds).", 1),
    "preprocess.reference_channel"                => simple_string_param("Channels(s) to use as reference", "avg"),
    "preprocess.layout.neighbour_criterion"       => number_param("Distance criterion (normalized) for channel neighbour definition.", 0.25, 0),
    "preprocess.eog.vEOG_channels"                => channel_groups_param("Channels used in the calculation of vertical eye movements (vEOG).", [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]]),
    "preprocess.eog.hEOG_channels"                => channel_groups_param("Channels used in the calculation of horizontal eye movements (hEOG).", [["F9"], ["F10"], ["hEOG"]]),
    "preprocess.eog.vEOG_criterion"               => number_param("Distance criterion for vertical EOG channel definition.", 50, 0),
    "preprocess.eog.hEOG_criterion"               => number_param("Distance criterion for horizontal EOG channel definition.", 30, 0),
    "preprocess.eeg.extreme_value_abs_criterion"  => number_param("Value (mV) for defining data section as an extreme value.", 500),
    "preprocess.eeg.artifact_value_abs_criterion" => number_param("Value (mV) for defining data section (or epoch) as an artifact value.", 100),
    "preprocess.eeg.artifact_value_z_criterion"   => number_param("Value (z) for defining data section (or epoch) as an artifact value (NB. various statistics with 0 being off!).", 0),

    # ICA settings
    "preprocess.ica.apply"              => bool_param("Independent Component Analysis (ICA) true/false."),
    "preprocess.ica.percentage_of_data" => number_param("Percentage of data to use for ICA (0-100).", 100.0, 0.0, 100.0),

    # Filtering settings - using helper function
    _filter_param_spec("preprocess.filter.highpass", true, "hp", 0.1, 0.01, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.lowpass", false, "lp", 30.0, 5.00, 500.0, 3, 1, 8)...,
    _filter_param_spec("preprocess.filter.ica_highpass", true, "hp", 1.0, 1.00, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.ica_lowpass", false, "lp", 30.0, 5.00, 500.0, 3, 1, 8)...,
)
# fmt: on

# =============================================================================
# VALIDATION TYPES
# =============================================================================

"""
    ValidationResult

Result of parameter validation.

# Fields
- `success::Bool`: Whether validation succeeded
- `error::Union{Nothing,String}`: Error message if validation failed, nothing if validation succeeded (default: nothing)
- `key_path::Union{Nothing,String}`: Path of TOML keys to the parameter that failed validation, nothing if validation succeeded (default: nothing)
"""
@kwdef struct ValidationResult
    success::Bool
    error::Union{Nothing,String} = nothing
    key_path::Union{Nothing,String} = nothing
end


# =============================================================================
# MAIN CONFIGURATION FUNCTIONS
# =============================================================================

"""
    load_config(config_file::String)

Load and merge configuration from a TOML file with defaults.

# Arguments
- `config_file::String`: Path to the configuration file

# Returns
- `Union{Dict,Nothing}`: The loaded configuration or nothing if loading failed
"""
function load_config(config_file::String)

    # Load default config
    default_config = TOML.parsefile(joinpath(@__DIR__, "default.toml"))
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

    # Merge, convert types, and validate
    config = _merge_configs(default_config, user_config)

    # Convert Any arrays to proper types (fixes Julia 1.12 TOML parsing issue)
    _convert_any_arrays!(config)

    validation_result = _validate_config(config)
    if !validation_result.success
        @minimal_error validation_result.error
        return nothing
    end

    return config
end


# =============================================================================
# TYPE CONVERSION FUNCTIONS
# =============================================================================

"""
    _convert_any_arrays!(config::Dict)

Convert Any arrays in the config to their proper types based on parameter specifications.
This fixes the issue where TOML.jl in Julia 1.12 returns Any arrays instead of typed arrays.
"""
function _convert_any_arrays!(config::Dict; path = "")
    for (key, value) in config
        new_path = isempty(path) ? key : "$path.$key"

        if isa(value, Dict)
            # Recursively process nested dictionaries
            _convert_any_arrays!(value; path = new_path)
        elseif haskey(PARAMETERS, new_path)
            # Convert this parameter if we have type information
            param_spec = PARAMETERS[new_path]
            param_type = typeof(param_spec).parameters[1]

            # Handle Vector{Vector{String}} case (like hEOG_channels, vEOG_channels)
            if param_type == Vector{Vector{String}} && isa(value, Vector) && eltype(value) == Any
                try
                    # Convert Any[["F9"], ["F10"]] -> Vector{Vector{String}}
                    converted_value = Vector{Vector{String}}()
                    for item in value
                        if isa(item, Vector)
                            push!(converted_value, String.(item))
                        else
                            push!(converted_value, [String(item)])
                        end
                    end
                    config[key] = converted_value
                catch e
                    @minimal_warning "Failed to convert $new_path from Any array to $param_type: $e"
                end
            elseif param_type <: Vector && isa(value, Vector) && eltype(value) == Any
                # Handle other Vector types
                try
                    inner_type = param_type.parameters[1]
                    config[key] = inner_type.(value)
                catch e
                    @minimal_warning "Failed to convert $new_path from Any array to $param_type: $e"
                end
            end
        end
    end
end

# =============================================================================
# CONFIGURATION MERGING FUNCTIONS
# =============================================================================

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
    _merge_nested!(result, user_config)
    return result
end

function _merge_nested!(target::Dict, source::Dict)
    for (key, value) in source
        if !haskey(target, key)
            target[key] = value
            continue
        end
        if isa(value, Dict) && isa(target[key], Dict)
            _merge_nested!(target[key], value)
        else
            target[key] = value
        end
    end
end


# =============================================================================
# CONFIGURATION VALIDATION FUNCTIONS
# =============================================================================

"""
    _validate_config(config::Dict; path="")

Validate config values against their metadata definitions.

# Arguments
- `config::Dict`: The configuration to validate
- `path::String`: The current path in the configuration (for nested validation)

# Returns
- `ValidationResult`: Result of the validation
"""
function _validate_config(config::Dict; path = "")
    for (key, value) in config
        new_path = isempty(path) ? key : "$path.$key"
        if isa(value, Dict) # Recursively validate nested dictionary
            result = _validate_config(value; path = new_path)
            !result.success && return result
        else # Check if we have data for this parameter
            if haskey(PARAMETERS, new_path)
                result = _validate_parameter(value, PARAMETERS[new_path], new_path)
                !result.success && return result
            else # Unknown parameter?
                return ValidationResult(success = false, error = "Unknown parameter: $new_path", key_path = new_path)
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

    # Helper function for creating validation errors
    function validation_error(msg)
        ValidationResult(success = false, error = msg, key_path = parameter_name)
    end

    # Check if value is the right type
    if param_type <: Number
        value isa Number || return validation_error("$parameter_name must be a number, got $(typeof(value))")
    else
        # Check type compatibility (fixed for Julia 1.12 TOML parsing)
        value isa param_type ||
            return validation_error("$parameter_name must be of type $param_type, got $(typeof(value))")
    end

    # Check min/max constraints
    !isnothing(parameter_spec.min) &&
        value < parameter_spec.min &&
        return validation_error("$parameter_name ($value) must be >= $(parameter_spec.min)")

    !isnothing(parameter_spec.max) &&
        value > parameter_spec.max &&
        return validation_error("$parameter_name ($value) must be <= $(parameter_spec.max)")

    # Check allowed values if they exist
    !isnothing(parameter_spec.allowed) &&
        !(value in parameter_spec.allowed) &&
        return validation_error("$parameter_name ($value) must be one of: $(join(parameter_spec.allowed, ", "))")

    return ValidationResult(success = true)
end




# =============================================================================
# PARAMETER INFORMATION DISPLAY FUNCTIONS
# =============================================================================

"""
    show_parameter_info(; parameter_name::String="")

Display information about configuration parameters. If parameter_name is empty, shows all parameters.
If parameter_name is provided, shows detailed information about that specific parameter.

# Arguments
- `parameter_name::String`: Optional path to a specific parameter (e.g., "filtering.highpass.cutoff")
"""
function show_parameter_info(; parameter_name::String = "")
    isempty(parameter_name) ? _show_all_parameters() : _show_specific_parameter(parameter_name)
end

function _show_all_parameters()
    @info "Available Configuration Parameters:"
    @info "==================================="

    sections = _group_parameters_by_section()
    sorted_sections = sort(collect(keys(sections)))

    for section in sorted_sections
        _display_section(section, sections[section])
    end

    @info "Use show_parameter_info(\"section\") for section overview"
    @info "Use show_parameter_info(\"section.parameter\") for specific parameter details"
end

function _display_section(section, section_data)
    @info "[$section]"
    @info "-"^(length(section) + 2)

    sorted_subsections = sort(collect(keys(section_data)))

    for subsection in sorted_subsections
        _display_subsection(subsection, section_data[subsection])
    end
end

function _display_subsection(subsection, params)
    !isempty(subsection) && @info "  [$subsection]"

    sorted_params = sort(params, by = first)
    for (path, parameter_spec) in sorted_params
        param_name = String(last(split(path, ".")))
        indent = isempty(subsection) ? "  " : "    "
        @info "$indent$param_name: $(parameter_spec.description)"
    end
end

function _show_specific_parameter(parameter_name)
    if haskey(PARAMETERS, parameter_name)
        _show_parameter_details(parameter_name)
    else
        matching_params = collect(filter(keys(PARAMETERS)) do key
            startswith(key, parameter_name)
        end)

        if !isempty(matching_params)
            _show_section_overview(parameter_name, matching_params)
        else
            @minimal_warning "Parameter or section not found: $parameter_name"
            @info "Use show_parameter_info() to see all available parameters and sections"
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
        min_str = isnothing(parameter_spec.min) ? "" : "$(parameter_spec.min) ≤ "
        max_str = isnothing(parameter_spec.max) ? "" : " ≤ $(parameter_spec.max)"
        @info "Range: $(min_str)value$(max_str)"
    end

    if !isnothing(parameter_spec.allowed)
        @info "Allowed values: $(join(parameter_spec.allowed, ", "))"
    end

    # Print default value or mark as required
    if isnothing(parameter_spec.default)
        @info "[REQUIRED]"
    else
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

    grouped_params = _group_params_by_subsection(section_name, matching_params)
    _display_grouped_params(grouped_params)

    @info ""
    @info "Use show_parameter_info(\"$section_name.parameter_name\") for detailed information about a specific parameter"
end

function _group_params_by_subsection(section_name::String, matching_params::Vector{String})
    sections = Dict{String,Vector{Tuple{String,ConfigParameter}}}()

    for param_path in matching_params
        subsection = _extract_subsection(section_name, param_path)
        get!(sections, subsection, Tuple{String,ConfigParameter}[])
        push!(sections[subsection], (param_path, PARAMETERS[param_path]))
    end

    return sections
end

function _extract_subsection(section_name::String, param_path::String)
    section_prefix = section_name * "."
    if !startswith(param_path, section_prefix)
        return ""
    end

    subsection_path = param_path[(length(section_prefix)+1):end]
    subsection_parts = split(subsection_path, ".")

    return length(subsection_parts) > 1 ? join(subsection_parts[1:(end-1)], ".") : ""
end

function _display_grouped_params(grouped_params::Dict{String,Vector{Tuple{String,ConfigParameter}}})
    sorted_subsections = sort(collect(keys(grouped_params)))

    for subsection in sorted_subsections
        !isempty(subsection) && @info "  [$subsection]"

        sorted_params = sort(grouped_params[subsection], by = first)
        for (path, parameter_spec) in sorted_params
            param_name = String(last(split(path, ".")))
            indent = isempty(subsection) ? "  " : "    "
            @info "$indent$param_name: $(parameter_spec.description)"
        end
    end
end

# =============================================================================
# TEMPLATE GENERATION FUNCTIONS
# =============================================================================

"""
    generate_config_template(filename::String="config_template.toml")

Generate and save a template TOML configuration file with all available parameters and their default values.

# Arguments
- `filename::String`: Name of the template file to create (default: "config_template.toml")
"""
function generate_config_template(; filename::String = "config_template.toml")
    try
        @info "Starting config template generation"
        open(filename, "w") do io
            _write_template_header(io)
            _write_template_sections(io)
        end
        @info "Configuration template saved to: $filename"
    catch e
        @minimal_error "Failed to save configuration template: $e"
    end
end

function _write_template_header(io::IO)
    println(io, "# EEG Processing Configuration Template")
    println(io, "# Generated on ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io)
    println(io, "# This template shows all available configuration options.")
    println(io, "# Required fields are marked with [REQUIRED]")
    println(io, "# Default values are shown where available")
    println(io)
end

function _write_template_sections(io::IO)
    sections = _group_parameters_by_section()
    sorted_sections = sort(collect(keys(sections)))

    for section in sorted_sections
        _write_section(io, section, sections[section])
    end
end

function _write_section(io::IO, section::String, section_data::Dict{String,Vector{Tuple{String,ConfigParameter}}})
    println(io, "\n# $section Settings")
    println(io, "[$section]")

    sorted_subsections = sort(collect(keys(section_data)))

    for subsection in sorted_subsections
        _write_subsection(io, section, subsection, section_data[subsection])
    end
end

function _write_subsection(io::IO, section::String, subsection::String, params::Vector{Tuple{String,ConfigParameter}})
    if !isempty(subsection)
        println(io, "\n# $subsection Settings")
        println(io, "[$section.$subsection]")
    end

    sorted_params = sort(params, by = first)
    for (path, parameter_spec) in sorted_params
        param_name = String(last(split(path, ".")))
        _write_parameter_docs(io, parameter_spec)
        _write_parameter_value(io, param_name, parameter_spec.default)
    end
end

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

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
        subsection = length(parts) > 2 ? join(parts[2:(end-1)], ".") : ""

        get!(sections, section, Dict{String,Vector{Tuple{String,ConfigParameter}}}())
        get!(sections[section], subsection, Tuple{String,ConfigParameter}[])
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
        min_str = isnothing(parameter_spec.min) ? "" : "$(parameter_spec.min) ≤ "
        max_str = isnothing(parameter_spec.max) ? "" : " ≤ $(parameter_spec.max)"
        println(io, "# Range: $(min_str)value$(max_str)")
    end

    if !isnothing(parameter_spec.allowed)
        println(io, "# Allowed values: $(join(parameter_spec.allowed, ", "))")
    end

    # Print default value or mark as required
    if isnothing(parameter_spec.default)
        println(io, "# [REQUIRED]")
    else
        println(io, "# Default: $(parameter_spec.default)")
    end
end

"""
    _write_parameter_value(io::IO, param_name::String, value)

Write a parameter value to the given IO stream in the appropriate TOML format.
"""
function _write_parameter_value(io::IO, param_name::String, value)
    if value isa String
        println(io, "$param_name = \"$(replace(value, "\\" => "\\\\"))\"")
    elseif value isa Vector
        println(io, "$param_name = $(isempty(value) ? "[]" : "[$(join(value, ", "))]")")
    else
        println(io, "$param_name = $value")
    end
end
