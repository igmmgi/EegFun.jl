# === TYPES AND STRUCTURES ===

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
@kwdef struct ConfigParameter{T}
    description::String
    default::Union{Nothing,T} = nothing
    allowed::Union{Nothing,Vector} = nothing
    min::Union{Nothing,T} = nothing
    max::Union{Nothing,T} = nothing
end

# === PARAMETER CONSTRUCTOR HELPERS ===

"""Create a `ConfigParameter{T}` with common defaults."""
function _param(::Type{T}, desc, default = nothing; allowed = nothing, min = nothing, max = nothing) where {T}
    ConfigParameter{T}(description = desc, default = default, allowed = allowed, min = min, max = max)
end

"""Create a string (or string-vector) parameter."""
string_param(desc, default = ""; allowed = nothing) = _param(Union{Vector{String},String}, desc, default, allowed = allowed)
"""Create a simple `String` parameter."""
simple_string_param(desc, default = ""; allowed = nothing) = _param(String, desc, default, allowed = allowed)
"""Create a Symbol parameter."""
symbol_param(desc, default; allowed = nothing) = _param(Symbol, desc, default, allowed = allowed)
"""Create a `Bool` parameter."""
bool_param(desc, default = false) = _param(Bool, desc, default)
"""Create a numeric parameter with optional min/max bounds."""
number_param(desc, default, min = nothing, max = nothing) = _param(Real, desc, default, min = min, max = max)
"""Create a numeric vector parameter (`Vector{Real}`)."""
number_vector_param(desc, default) = _param(Vector{Real}, desc, default)
"""Create a channel-groups parameter (`Vector{Vector{String}}`)."""
channel_groups_param(desc, default) = _param(Vector{Vector{String}}, desc, default)

"""Generate the standard set of filter parameters (apply, type, method, func, freq, order) for a given prefix."""
function _filter_param_spec(prefix, apply, type, freq, min_freq, max_freq, order, min_order, max_order)
    Dict(
        "$prefix.apply" => bool_param("Apply: true/false", apply),
        "$prefix.type" => symbol_param("Filter type identifier", type, allowed = [:hp, :lp]),
        "$prefix.method" => symbol_param("Filter method", :iir, allowed = [:fir, :iir]),
        "$prefix.func" => symbol_param("Filter function", :filtfilt, allowed = [:filt, :filtfilt]),
        "$prefix.freq" => number_param("Cutoff frequency (Hz)", freq, min_freq, max_freq),
        "$prefix.order" => number_param("Filter order", order, min_order, max_order),
    )
end


# === PARAMETER DEFINITIONS ===

# fmt: off
const PARAMETERS = Dict{String,ConfigParameter}(

    # File paths and settings
    "files.input.directory" => simple_string_param("Directory containing raw data files.", "."),
    "files.input.raw_data_files" => string_param("Pattern (regex or explicit list) for raw data files to process.", "\\.bdf"),
    "files.input.recursive" => bool_param("Search subdirectories recursively for raw data files (e.g., BIDS).", false),
    "files.input.layout_file" => simple_string_param("Electrode layout file name (\"*.csv\")", "biosemi72.csv"),
    "files.input.epoch_condition_file" => simple_string_param("TOML file that defines the condition epochs.", ""),
    "files.output.directory" => simple_string_param("Directory for processed output files", "./preprocessed_files"),

    # What data should we save?
    "files.output.save_continuous_data_raw" => bool_param("Save continuous data original?", false),
    "files.output.save_continuous_data_corrected" => bool_param("Save continuous data cleaned?", false),
    "files.output.save_ica_data" => bool_param("Save ICA results?", true),
    "files.output.save_epoch_data_raw" => bool_param("Save epoched data original?", false),
    "files.output.save_epoch_data_corrected" => bool_param("Save epoched data cleaned?", false),
    "files.output.save_epoch_data" => bool_param("Save epoched data good?", true),
    "files.output.save_erp_data_raw" => bool_param("Save ERP data original?", false),
    "files.output.save_erp_data_corrected" => bool_param("Save ERP data cleaned?", false),
    "files.output.save_erp_data" => bool_param("Save ERP data good?", true),

    # Preprocessing settings
    "preprocess.interactive_continuous" => bool_param("Pause execution to interactively review continuous data", false),
    "preprocess.interactive_ica" => bool_param("Pause execution to interactively review ICA components", false),
    "preprocess.interactive_epochs" => bool_param("Pause execution to interactively review epoch rejection", false),
    "preprocess.epoch_start" => number_param("Epoch start (seconds).", -1),
    "preprocess.epoch_end" => number_param("Epoch end (seconds).", 1),
    "preprocess.reference_channel" => simple_string_param("Channels(s) to use as reference", "avg"),
    "preprocess.layout.neighbour_criterion" =>
        number_param("Distance criterion (normalized) for channel neighbour definition.", 0.25, 0),
    "preprocess.eog.vEOG_channels" => channel_groups_param(
        "Channels used in the calculation of vertical eye movements (vEOG).",
        [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    ),
    "preprocess.eog.hEOG_channels" =>
        channel_groups_param("Channels used in the calculation of horizontal eye movements (hEOG).", [["F9"], ["F10"], ["hEOG"]]),
    "preprocess.eog.vEOG_criterion" => number_param("Distance criterion for vertical EOG channel definition.", 50, 0),
    "preprocess.eog.hEOG_criterion" => number_param("Distance criterion for horizontal EOG channel definition.", 30, 0),
    "preprocess.eeg.extreme_value_abs_criterion" => number_param("Value (mV) for defining data section as an extreme value.", 500),
    "preprocess.eeg.artifact_value_abs_criterion" =>
        number_param("Value (mV) for defining data section (or epoch) as an artifact value.", 100),
    "preprocess.eeg.artifact_value_z_criterion" => number_param(
        "Value (z) for defining data section (or epoch) as an artifact value (NB. various statistics with 0 being off!).",
        0,
    ),
    "preprocess.eeg.artifact_interval_start" => number_param("Start time (s) for artifact rejection (optional).", nothing),
    "preprocess.channel_repair.method" =>
        simple_string_param("Method for bad channel interpolation (:spherical_spline or :neighbor_interpolation)", "spherical_spline"),

    # ICA settings
    "preprocess.ica.apply" => bool_param("Independent Component Analysis (ICA) true/false.", true),
    "preprocess.ica.percentage_of_data" => number_param("Percentage of data to use for ICA (0-100).", 100.0, 0.0, 100.0),

    # CleanLine settings
    "preprocess.cleanline.apply" => bool_param("Apply CleanLine algorithm to remove line noise?", false),
    "preprocess.cleanline.line_frequencies" => number_vector_param("Line noise frequencies to target (e.g. [50.0])", [50.0]),
    "preprocess.cleanline.bandwidth" => number_param("Bandwidth for scanning frequencies around the line frequency.", 2.0),
    "preprocess.cleanline.sliding_win_length" => number_param("Sliding window length (seconds) for multi-taper regression.", 4.0),
    "preprocess.cleanline.sliding_win_step" => number_param("Sliding window step (seconds) for multi-taper regression.", 2.0),
    "preprocess.cleanline.time_bandwidth" => number_param("Time-bandwidth product (TW) for tapers.", 3.0),
    "preprocess.cleanline.k_tapers" => number_param("Number of tapers (usually 2*TW-1).", 5),
    "preprocess.cleanline.p_value" => number_param("Significance threshold for F-test.", 0.05),
    "preprocess.cleanline.pad" => number_param("Padding factor for FFT.", 2),

    # Resampling settings
    "preprocess.resample.apply" => bool_param("Apply resampling/downsampling?", false),
    "preprocess.resample.target_rate" => number_param("Target sampling rate in Hz (e.g. 512, 256).", 512),

    # Filtering settings - using helper function
    _filter_param_spec("preprocess.filter.highpass", true, :hp, 0.1, 0.01, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.lowpass", false, :lp, 30.0, 5.00, 500.0, 3, 1, 8)...,
    _filter_param_spec("preprocess.filter.ica_highpass", true, :hp, 1.0, 1.00, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.ica_lowpass", false, :lp, 30.0, 5.00, 500.0, 3, 1, 8)...,
)
# fmt: on

# === VALIDATION TYPES ===

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


# === MAIN CONFIGURATION FUNCTIONS ===

"""
    read_config(config_file::String)

Load and merge configuration from a TOML file with defaults.

# Arguments
- `config_file::String`: Path to the configuration file

# Returns
- `Union{Dict,Nothing}`: The loaded configuration or nothing if loading failed
"""
function _generate_default_config()
    config = Dict{String,Any}()
    for (path, param) in PARAMETERS
        isnothing(param.default) && continue

        parts = split(path, ".")
        current = config
        for i = 1:(length(parts)-1)
            key = String(parts[i])
            if !haskey(current, key)
                current[key] = Dict{String,Any}()
            end
            current = current[key]
        end
        current[String(parts[end])] = param.default
    end
    return config
end

"""
    read_config(config_file::String)

Read a TOML configuration file into a dictionary.
"""
function read_config(config_file::String)
    default_config = _generate_default_config()

    if !isfile(config_file)
        @minimal_error "Configuration file not found: $config_file"
    end

    user_config = Dict{String,Any}()
    @info "Loading config file: $config_file"
    try
        user_config = TOML.parsefile(config_file)
    catch e
        @minimal_error "Error parsing TOML file: $e"
    end

    # Merge, convert types, and validate
    config = _merge_configs(default_config, user_config)

    # Convert Any arrays to proper types (fixes Julia 1.12 TOML parsing issue)
    _convert_any_arrays!(config)

    validation_result = _validate_config(config)
    if !validation_result.success
        @minimal_error validation_result.error
    end

    return config
end


# === TYPE CONVERSION FUNCTIONS ===

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
            elseif param_type == Symbol && isa(value, AbstractString)
                config[key] = Symbol(value)
            elseif (param_type <: Vector || param_type == Vector{Real}) && isa(value, Vector)
                # Handle other Vector types (including Vector{Real})
                try
                    inner_type = param_type == Vector{Real} ? Float64 : param_type.parameters[1]
                    config[key] = inner_type.(value)
                catch e
                    @minimal_warning "Failed to convert $new_path from $(typeof(value)) to $param_type: $e"
                end
            end
        end
    end
end

# === CONFIGURATION MERGING FUNCTIONS ===

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

"""Recursively merge `source` into `target`, overwriting leaf values."""
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


# === CONFIGURATION VALIDATION FUNCTIONS ===

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
function _validate_parameter(value, parameter_spec::ConfigParameter{T}, parameter_name::String) where {T}
    param_type = T

    # Helper function for creating validation errors
    """Create a failed `ValidationResult` for this parameter."""
    function validation_error(msg)
        ValidationResult(success = false, error = msg, key_path = parameter_name)
    end

    # Check if value is the right type
    if param_type <: Number
        value isa Number || return validation_error("$parameter_name must be a number, got $(typeof(value))")
    elseif param_type == Vector{Real} || param_type <: AbstractVector{<:Real}
        (value isa AbstractVector && eltype(value) <: Real) ||
            return validation_error("$parameter_name must be a vector of numbers, got $(typeof(value))")
    else
        # Check type compatibility (fixed for Julia 1.12 TOML parsing)
        value isa param_type || return validation_error("$parameter_name must be of type $param_type, got $(typeof(value))")
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




# === PARAMETER INFORMATION DISPLAY FUNCTIONS ===

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

"""Print an overview of every configuration parameter, grouped by section."""
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

"""Print all parameters in one top-level section."""
function _display_section(section, section_data)
    @info "[$section]"
    @info "-"^(length(section) + 2)

    sorted_subsections = sort(collect(keys(section_data)))

    for subsection in sorted_subsections
        _display_subsection(subsection, section_data[subsection])
    end
end

"""Print all parameters in one subsection."""
function _display_subsection(subsection, params)
    !isempty(subsection) && @info "  [$subsection]"

    sorted_params = sort(params, by = first)
    for (path, parameter_spec) in sorted_params
        param_name = String(last(split(path, ".")))
        indent = isempty(subsection) ? "  " : "    "
        @info "$indent$param_name: $(parameter_spec.description)"
    end
end

"""Show details or section overview for a specific parameter name."""
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
    _show_parameter_details(parameter_name, PARAMETERS[parameter_name])
end

function _show_parameter_details(parameter_name::String, parameter_spec::ConfigParameter{T}) where {T}
    @info "Parameter: $parameter_name"
    @info "="^(length(parameter_name) + 11)
    @info "Description: $(parameter_spec.description)"
    @info "Type: $T"

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

"""Group matching parameters into a dict keyed by subsection name."""
function _group_params_by_subsection(section_name::String, matching_params::Vector{String})
    sections = Dict{String,Vector{Tuple{String,ConfigParameter}}}()

    for param_path in matching_params
        subsection = _extract_subsection(section_name, param_path)
        get!(sections, subsection, Tuple{String,ConfigParameter}[])
        push!(sections[subsection], (param_path, PARAMETERS[param_path]))
    end

    return sections
end

"""Extract the subsection portion of a parameter path relative to the section prefix."""
function _extract_subsection(section_name::String, param_path::String)
    section_prefix = section_name * "."
    if !startswith(param_path, section_prefix)
        return ""
    end

    subsection_path = param_path[(length(section_prefix)+1):end]
    subsection_parts = split(subsection_path, ".")

    return length(subsection_parts) > 1 ? join(subsection_parts[1:(end-1)], ".") : ""
end

"""Print parameters grouped by subsection."""
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

# === TEMPLATE GENERATION FUNCTIONS ===

"""
    generate_config_template(; filename::String="config_template.toml")

Generate and save a template TOML configuration file with all available parameters and their default values.

# Arguments
- `filename::String`: Name of the template file to create (keyword argument, default: "config_template.toml")
"""
function generate_config_template(; filename::String = "config_template.toml")
    if !endswith(filename, ".toml")
        filename = filename * ".toml"
    end
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

"""Write the TOML template file header comment block."""
function _write_template_header(io::IO)
    println(io, "# EEG Processing Configuration Template")
    println(io, "# Generated on ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io)
    println(io, "# This template shows all available configuration options.")
    println(io, "# Required fields are marked with [REQUIRED]")
    println(io, "# Default values are shown where available")
    println(io)
end

"""Write all parameter sections to the TOML template file."""
function _write_template_sections(io::IO)
    sections = _group_parameters_by_section()
    sorted_sections = sort(collect(keys(sections)))

    for section in sorted_sections
        _write_section(io, section, sections[section])
    end
end

"""Write a single top-level section and its subsections to the TOML template."""
function _write_section(io::IO, section::String, section_data::Dict{String,Vector{Tuple{String,ConfigParameter}}})
    println(io, "\n# $section Settings")
    println(io, "[$section]")

    sorted_subsections = sort(collect(keys(section_data)))

    for subsection in sorted_subsections
        _write_subsection(io, section, subsection, section_data[subsection])
    end
end

"""Write a subsection header and its parameter entries to the TOML template."""
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

# === UTILITY FUNCTIONS ===

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
function _write_parameter_docs(io::IO, parameter_spec::ConfigParameter{T}) where {T}
    println(io, "\n# $(parameter_spec.description)")
    println(io, "# Type: $T")

    if !isnothing(parameter_spec.min) || !isnothing(parameter_spec.max)
        min_str = isnothing(parameter_spec.min) ? "" : "$(parameter_spec.min) ≤ "
        max_str = isnothing(parameter_spec.max) ? "" : " ≤ $(parameter_spec.max)"
        println(io, "# Range: $(min_str)value$(max_str)")
    end

    if !isnothing(parameter_spec.allowed)
        println(io, "# Allowed values: $(join(parameter_spec.allowed, ", "))")
    end

    # Print default value or mark as required/optional
    if isnothing(parameter_spec.default)
        if occursin("(optional)", lowercase(parameter_spec.description))
            println(io, "# [OPTIONAL]")
        else
            println(io, "# [REQUIRED]")
        end
    else
        println(io, "# Default: $(parameter_spec.default)")
    end
end

"""
    _format_toml_value(value) -> String

Recursively format a value for TOML output.
"""
function _format_toml_value(value)
    if value isa String
        return "\"$(replace(value, "\\" => "\\\\"))\""
    elseif value isa Symbol
        return "\"$(value)\""
    elseif value isa Vector
        return isempty(value) ? "[]" : "[" * join([_format_toml_value(v) for v in value], ", ") * "]"
    elseif value isa Bool
        return value ? "true" : "false"
    else
        return string(value)
    end
end

"""
    _write_parameter_value(io::IO, param_name::String, value)

Write a parameter value to the given IO stream in the appropriate TOML format.
"""
function _write_parameter_value(io::IO, param_name::String, value)
    if isnothing(value)
        println(io, "# $param_name = ")
    else
        println(io, "$param_name = $(_format_toml_value(value))")
    end
end
