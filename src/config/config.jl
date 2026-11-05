"""
    ConfigParameter{T}

Specification for a typed configuration parameter with optional default, allowed values, and numeric bounds.
"""
@kwdef struct ConfigParameter{T}
    description::String
    default::Union{Nothing,T} = nothing
    allowed::Union{Nothing,Vector{Any}} = nothing
    min::Union{Nothing,T} = nothing
    max::Union{Nothing,T} = nothing
end

function _param(::Type{T}, desc, default = nothing; allowed = nothing, min = nothing, max = nothing) where {T}
    ConfigParameter{T}(description = desc, default = default, allowed = allowed, min = min, max = max)
end

string_param(desc, default = ""; allowed = nothing) = _param(Union{Vector{String},String}, desc, default; allowed = allowed)
simple_string_param(desc, default = ""; allowed = nothing) = _param(String, desc, default; allowed = allowed)
symbol_param(desc, default; allowed = nothing) = _param(Symbol, desc, default; allowed = allowed)
bool_param(desc, default = false) = _param(Bool, desc, default)
number_param(desc, default, min = nothing, max = nothing) = _param(Real, desc, default; min = min, max = max)
number_vector_param(desc, default) = _param(Vector{Real}, desc, default)
channel_groups_param(desc, default) = _param(Vector{Vector{String}}, desc, default)

function _filter_param_spec(prefix, apply, type, freq, min_freq, max_freq, order, min_order, max_order)
    Dict(
        "$prefix.apply" => bool_param("Apply: true/false", apply),
        "$prefix.type" => symbol_param("Filter type identifier", type; allowed = [:hp, :lp]),
        "$prefix.method" => symbol_param("Filter method", :iir; allowed = [:fir, :iir]),
        "$prefix.func" => symbol_param("Filter function", :filtfilt; allowed = [:filt, :filtfilt]),
        "$prefix.freq" => number_param("Cutoff frequency (Hz)", freq, min_freq, max_freq),
        "$prefix.order" => number_param("Filter order", order, min_order, max_order),
    )
end

# fmt: off
const PARAMETERS = Dict{String,ConfigParameter}(
    # File paths and settings
    "files.input.directory" => simple_string_param("Directory containing raw data files.", "."),
    "files.input.raw_data_files" => string_param("Pattern (regex or explicit list) for raw data files to process.", "\\.bdf"),
    "files.input.recursive" => bool_param("Search subdirectories recursively for raw data files (e.g., BIDS).", false),
    "files.input.layout_file" => simple_string_param("Electrode layout file name (\"*.csv\")", "biosemi72.csv"),
    "files.input.epoch_condition_file" => simple_string_param("TOML file that defines the condition epochs.", ""),
    "files.output.directory" => simple_string_param("Directory for processed output files", "./preprocessed_files"),

    # Output saving flags
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
    "preprocess.reference_channel" => string_param("Channel(s) to use as reference", "avg"),
    "preprocess.layout.neighbor_criterion" => number_param("Distance criterion (normalized) for channel neighbor definition.", 0.35, 0),
    "preprocess.eog.vEOG_channels" => channel_groups_param(
        "Channels used in the calculation of vertical eye movements (vEOG).",
        [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    ),
    "preprocess.eog.hEOG_channels" =>
        channel_groups_param("Channels used in the calculation of horizontal eye movements (hEOG).", [["F9"], ["F10"], ["hEOG"]]),
    "preprocess.eog.vEOG_criterion" => number_param("Distance criterion for vertical EOG channel definition.", 50, 0),
    "preprocess.eog.hEOG_criterion" => number_param("Distance criterion for horizontal EOG channel definition.", 30, 0),
    "preprocess.eeg.extreme_value_abs_criterion" => number_param("Value (μV) for defining data section as an extreme value.", 500),
    "preprocess.eeg.artifact_value_abs_criterion" =>
        number_param("Value (μV) for defining data section (or epoch) as an artifact value.", 100),
    "preprocess.eeg.artifact_value_z_criterion" => number_param(
        "Value (z) for defining data section (or epoch) as an artifact value (NB. various statistics with 0 being off!).",
        0,
    ),
    "preprocess.eeg.artifact_interval_start" => number_param("Start time (s) for artifact rejection (optional).", nothing),
    "preprocess.eeg.artifact_interval_end" => number_param("End time (s) for artifact rejection (optional).", nothing),
    "preprocess.channel_repair.method" =>
        simple_string_param("Method for bad channel interpolation (:spherical_spline or :neighbor_interpolation)", "spherical_spline"),

    # Autoreject adaptive artifact rejection settings
    "preprocess.autoreject.apply" => bool_param("Apply Autoreject adaptive artifact rejection and repair?", false),
    "preprocess.autoreject.mode" => symbol_param("Autoreject mode (:local or :global).", :local; allowed = [:local, :global]),
    "preprocess.autoreject.n_folds" => number_param("Number of cross-validation folds for Autoreject.", 10, 2, 50),
    "preprocess.autoreject.n_candidates" => number_param("Number of candidate thresholds to evaluate per channel.", 30, 5, 100),
    "preprocess.autoreject.consensus_perc" => number_param("Consensus rejection threshold (fraction of bad channels).", 0.2, 0.0, 1.0),
    "preprocess.autoreject.max_repair_perc" => number_param("Maximum fraction of bad channels to repair in an epoch.", 0.1, 0.0, 1.0),
    "preprocess.autoreject.repair_method" => symbol_param(
        "Channel repair method (:spherical_spline or :neighbor_interpolation).",
        :spherical_spline;
        allowed = [:spherical_spline, :neighbor_interpolation],
    ),
    "preprocess.autoreject.seed" => number_param("Random seed for cross-validation splits.", 42),

    # ICA settings
    "preprocess.ica.apply" => bool_param("Independent Component Analysis (ICA) true/false.", true),
    "preprocess.ica.percentage_of_data" => number_param("Percentage of data to use for ICA (0-100).", 100.0, 0.0, 100.0),
    "preprocess.ica.component_method" =>
        symbol_param("Method for identifying artifact components (:correlation).", :correlation; allowed = [:correlation]),

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

    # Filtering settings
    _filter_param_spec("preprocess.filter.highpass", true, :hp, 0.1, 0.01, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.lowpass", false, :lp, 30.0, 5.00, 500.0, 3, 1, 8)...,
    _filter_param_spec("preprocess.filter.ica_highpass", true, :hp, 1.0, 1.00, 20.0, 1, 1, 4)...,
    _filter_param_spec("preprocess.filter.ica_lowpass", false, :lp, 30.0, 5.00, 500.0, 3, 1, 8)...,
)
# fmt: on

"""
    ValidationResult

Result of configuration parameter validation (`success`, `error`, `key_path`).
"""
@kwdef struct ValidationResult
    success::Bool
    error::Union{Nothing,String} = nothing
    key_path::Union{Nothing,String} = nothing
end

function _generate_default_config()
    config = Dict{String,Any}()
    for (path, param) in PARAMETERS
        def = param.default
        isnothing(def) && continue

        current = config
        start_idx = 1
        while true
            dot_idx = findnext('.', path, start_idx)
            if isnothing(dot_idx)
                leaf_key = SubString(path, start_idx)
                current[String(leaf_key)] = def isa Vector ? deepcopy(def) : def
                break
            end
            seg = String(SubString(path, start_idx, dot_idx - 1))
            current = get!(() -> Dict{String,Any}(), current, seg)::Dict{String,Any}
            start_idx = dot_idx + 1
        end
    end
    return config
end

"""
    read_config(config_file::String) -> Dict{String,Any}

Load, merge with defaults, and validate a TOML configuration file.
"""
function read_config(config_file::String)
    !isfile(config_file) && @minimal_error "Configuration file not found: $config_file"
    default_config = _generate_default_config()

    @info "Loading config file: $config_file"
    user_config = try
        TOML.parsefile(config_file)
    catch e
        @minimal_error "Error parsing TOML file: $(sprint(showerror, e))"
    end

    config = _merge_configs(default_config, user_config)
    _convert_any_arrays!(config)

    validation_result = _validate_config(config)
    !validation_result.success && @minimal_error validation_result.error
    return config
end

function _convert_param_value(::ConfigParameter{T}, value, new_path::String) where {T}
    if T == Vector{Vector{String}} && value isa Vector && eltype(value) == Any
        try
            return [item isa Vector ? String.(item) : [String(item)] for item in value]
        catch e
            @minimal_error "Failed to convert $new_path from Any array to $T: $(sprint(showerror, e))"
        end
    elseif T == Symbol && value isa AbstractString
        return Symbol(value)
    elseif (T <: Vector || T == Vector{Real}) && value isa Vector
        try
            inner_type = T == Vector{Real} ? Float64 : eltype(T)
            return inner_type.(value)
        catch e
            @minimal_error "Failed to convert $new_path from $(typeof(value)) to $T: $(sprint(showerror, e))"
        end
    end
    return value
end

function _convert_any_arrays!(config::Dict; path::String = "")
    for (key, value) in config
        new_path = isempty(path) ? key : "$path.$key"
        if value isa Dict
            _convert_any_arrays!(value; path = new_path)
        elseif haskey(PARAMETERS, new_path)
            config[key] = _convert_param_value(PARAMETERS[new_path], value, new_path)
        end
    end
end

function _merge_configs(default_config::Dict, user_config::Dict)
    result = copy(default_config)
    _merge_nested!(result, user_config)
    return result
end

function _merge_nested!(target::Dict, source::Dict)
    for (key, value) in source
        if haskey(target, key) && value isa Dict && target[key] isa Dict
            _merge_nested!(target[key], value)
        else
            target[key] = value
        end
    end
end

function _validate_config(config::Dict; path::String = "")
    for (key, value) in config
        isempty(path) && key == "metadata" && value isa Dict && continue
        new_path = isempty(path) ? key : "$path.$key"
        if value isa Dict
            result = _validate_config(value; path = new_path)
            !result.success && return result
        elseif haskey(PARAMETERS, new_path)
            result = _validate_parameter(value, PARAMETERS[new_path], new_path)
            !result.success && return result
        else
            return ValidationResult(success = false, error = "Unknown parameter: $new_path", key_path = new_path)
        end
    end
    return ValidationResult(success = true)
end

function _validate_parameter(value, parameter_spec::ConfigParameter{T}, parameter_name::String) where {T}
    fail(msg) = ValidationResult(success = false, error = msg, key_path = parameter_name)

    if T <: Number
        value isa Number || return fail("$parameter_name must be a number, got $(typeof(value))")
    elseif T == Vector{Real} || T <: AbstractVector{<:Real}
        (value isa AbstractVector && eltype(value) <: Real) ||
            return fail("$parameter_name must be a vector of numbers, got $(typeof(value))")
    else
        value isa T || return fail("$parameter_name must be of type $T, got $(typeof(value))")
    end

    !isnothing(parameter_spec.min) && value < parameter_spec.min && return fail("$parameter_name ($value) must be >= $(parameter_spec.min)")

    !isnothing(parameter_spec.max) && value > parameter_spec.max && return fail("$parameter_name ($value) must be <= $(parameter_spec.max)")

    !isnothing(parameter_spec.allowed) &&
        !(value in parameter_spec.allowed) &&
        return fail("$parameter_name ($value) must be one of: $(join(parameter_spec.allowed, ", "))")

    return ValidationResult(success = true)
end

"""
    show_parameter_info(; parameter_name::String = "")

Display documentation for all configuration parameters, a section prefix, or a single parameter.
"""
function show_parameter_info(; parameter_name::String = "")
    isempty(parameter_name) ? _show_all_parameters() : _show_specific_parameter(parameter_name)
end

function _show_all_parameters()
    @info "Available Configuration Parameters:"
    @info "==================================="
    sections = _group_parameters_by_section()
    for section in sort!(collect(keys(sections)))
        _display_section(section, sections[section])
    end
    @info "Use show_parameter_info(\"section\") for section overview"
    @info "Use show_parameter_info(\"section.parameter\") for specific parameter details"
end

function _display_section(section, section_data)
    @info "[$section]"
    @info "-"^(length(section) + 2)
    _display_grouped_params(section_data)
end

@inline function _leaf_param_name(path::String)
    idx = findlast('.', path)
    return isnothing(idx) ? path : path[(idx+1):end]
end

function _display_subsection(subsection, params)
    !isempty(subsection) && @info "  [$subsection]"
    indent = isempty(subsection) ? "  " : "    "
    for (path, parameter_spec) in sort(params; by = first)
        @info "$indent$(_leaf_param_name(path)): $(parameter_spec.description)"
    end
end

function _display_grouped_params(grouped_params::Dict{String,Vector{Tuple{String,ConfigParameter}}})
    for subsection in sort!(collect(keys(grouped_params)))
        _display_subsection(subsection, grouped_params[subsection])
    end
end

function _show_specific_parameter(parameter_name)
    if haskey(PARAMETERS, parameter_name)
        _show_parameter_details(parameter_name)
    else
        matching_params = [k for k in keys(PARAMETERS) if startswith(k, parameter_name)]
        if !isempty(matching_params)
            _show_section_overview(parameter_name, matching_params)
        else
            @minimal_warning "Parameter or section not found: $parameter_name"
            @info "Use show_parameter_info() to see all available parameters and sections"
        end
    end
end

_show_parameter_details(parameter_name::String) = _show_parameter_details(parameter_name, PARAMETERS[parameter_name])

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
    !isnothing(parameter_spec.allowed) && @info "Allowed values: $(join(parameter_spec.allowed, ", "))"
    @info isnothing(parameter_spec.default) ? "[REQUIRED]" : "Default: $(parameter_spec.default)"
end

function _show_section_overview(section_name::String, matching_params::Vector{String})
    @info "Section: $section_name"
    @info "="^(length(section_name) + 9)
    _display_grouped_params(_group_params_by_subsection(section_name, matching_params))
    @info ""
    @info "Use show_parameter_info(\"$section_name.parameter_name\") for detailed information about a specific parameter"
end

function _group_params_by_subsection(section_name::String, matching_params::Vector{String})
    sections = Dict{String,Vector{Tuple{String,ConfigParameter}}}()
    for param_path in matching_params
        subsection = _extract_subsection(section_name, param_path)
        push!(get!(() -> Tuple{String,ConfigParameter}[], sections, subsection), (param_path, PARAMETERS[param_path]))
    end
    return sections
end

function _extract_subsection(section_name::String, param_path::String)
    n_sec = ncodeunits(section_name)
    (ncodeunits(param_path) > n_sec && startswith(param_path, section_name) && codeunit(param_path, n_sec + 1) == UInt8('.')) || return ""
    last_dot = findlast('.', param_path)
    return (isnothing(last_dot) || last_dot <= n_sec + 1) ? "" : param_path[(n_sec+2):(last_dot-1)]
end

"""
    generate_config_template(; filename::String = "config_template.toml")

Generate and save a commented TOML configuration template with all parameters and defaults.
"""
function generate_config_template(; filename::String = "config_template.toml")
    endswith(filename, ".toml") || (filename *= ".toml")
    try
        @info "Starting config template generation"
        buf = IOBuffer(; sizehint = 16384)
        _write_template_header(buf)
        _write_template_sections(buf)
        write(filename, take!(buf))
        @info "Configuration template saved to: $filename"
    catch e
        @minimal_error "Failed to save configuration template: $(sprint(showerror, e))"
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
    for section in sort!(collect(keys(sections)))
        _write_section(io, section, sections[section])
    end
end

function _write_section(io::IO, section::String, section_data::Dict{String,Vector{Tuple{String,ConfigParameter}}})
    println(io, "\n# $section Settings")
    println(io, "[$section]")
    for subsection in sort!(collect(keys(section_data)))
        _write_subsection(io, section, subsection, section_data[subsection])
    end
end

function _write_subsection(io::IO, section::String, subsection::String, params::Vector{Tuple{String,ConfigParameter}})
    if !isempty(subsection)
        println(io, "\n# $subsection Settings")
        println(io, "[$section.$subsection]")
    end
    for (path, parameter_spec) in sort(params; by = first)
        _write_parameter_docs(io, parameter_spec)
        _write_parameter_value(io, _leaf_param_name(path), parameter_spec.default)
    end
end

function _group_parameters_by_section()
    sections = Dict{String,Dict{String,Vector{Tuple{String,ConfigParameter}}}}()
    for (path, parameter_spec) in PARAMETERS
        first_dot = findfirst('.', path)
        last_dot = findlast('.', path)
        section = isnothing(first_dot) ? path : path[1:(first_dot-1)]
        subsection = (isnothing(first_dot) || first_dot == last_dot) ? "" : path[(first_dot+1):(last_dot-1)]

        sec_dict = get!(() -> Dict{String,Vector{Tuple{String,ConfigParameter}}}(), sections, section)
        push!(get!(() -> Tuple{String,ConfigParameter}[], sec_dict, subsection), (path, parameter_spec))
    end
    return sections
end

function _write_parameter_docs(io::IO, parameter_spec::ConfigParameter{T}) where {T}
    println(io, "\n# $(parameter_spec.description)")
    println(io, "# Type: $T")

    if !isnothing(parameter_spec.min) || !isnothing(parameter_spec.max)
        min_str = isnothing(parameter_spec.min) ? "" : "$(parameter_spec.min) ≤ "
        max_str = isnothing(parameter_spec.max) ? "" : " ≤ $(parameter_spec.max)"
        println(io, "# Range: $(min_str)value$(max_str)")
    end

    !isnothing(parameter_spec.allowed) && println(io, "# Allowed values: $(join(parameter_spec.allowed, ", "))")

    if isnothing(parameter_spec.default)
        println(io, occursin("(optional)", lowercase(parameter_spec.description)) ? "# [OPTIONAL]" : "# [REQUIRED]")
    else
        println(io, "# Default: $(parameter_spec.default)")
    end
end

function _format_toml_value(value)
    value isa String && return "\"$(replace(value, "\\" => "\\\\"))\""
    value isa Symbol && return "\"$(value)\""
    if value isa Vector
        return isempty(value) ? "[]" : "[" * join((_format_toml_value(v) for v in value), ", ") * "]"
    end
    value isa Bool && return value ? "true" : "false"
    return string(value)
end

function _write_parameter_value(io::IO, param_name::String, value)
    if isnothing(value)
        println(io, "# $param_name = ")
    else
        println(io, param_name, " = ", _format_toml_value(value))
    end
end
