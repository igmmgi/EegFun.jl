# Custom exception type with cleaner error display
struct ConfigurationError <: Exception
    msg::String
end

# Custom display method to avoid printing a stacktrace
function Base.showerror(io::IO, e::ConfigurationError)
    print(io, "Configuration Error: ", e.msg)
end
struct MinimalError <: Exception
    msg::String
end

# Custom display method to avoid printing a stacktrace
function Base.showerror(io::IO, e::MinimalError)
    print(io, "Configuration Error: ", e.msg)
end

# Override stack trace printing specifically for MinimalError
function Base.show_backtrace(io::IO, bt::Vector, ::MinimalError)
    # Do nothing - suppress backtrace
end

"""
    @minimal_error(msg)

Display an error message without a stacktrace and return nothing.
"""
macro minimal_error(msg)
    quote
        @error "Configuration Error: ", $(esc(msg)) _module=nothing _file=nothing _line=nothing
        throw(ErrorException("")) # Throw a silent exception that will be caught 
    end
end

"""
    validate_config_with_minimal_errors(config::Dict)

Validates the configuration but catches and handles errors with minimal output.
Returns true if validation passed, false if validation failed.
"""
function validate_config_with_minimal_errors(config::Dict)
    try
        validate_config(config)
        return true
    catch e
        if e isa ErrorException && isempty(e.msg)
            # This is from our minimal_error macro - already printed message, no need for stacktrace
            return false
        else
            # Some other unexpected error - rethrow to get full stack trace
            rethrow()
        end
    end
end

"""
    load_config(config_file::String; validate::Bool=true)

Load and merge configuration from a TOML file with defaults.
"""
function load_config(config_file::String; validate::Bool=true)
    # Load default config
    default_config = TOML.parsefile(joinpath(@__DIR__, "default.toml"))
    
    # Check if user config exists
    if !isfile(config_file)
        println(stderr, "Configuration Error: Configuration file not found: $config_file")
        return nothing
    end
    
    # Load user config
    user_config = Dict()
    try
        @info "Loading config file: $config_file"
        user_config = TOML.parsefile(config_file)
    catch e
        println(stderr, "Configuration Error: Error parsing TOML file: $(e.msg)")
        return nothing
    end
    
    # Merge and validate
    config = merge_configs(default_config, user_config)
    
    if validate
        if !validate_config_with_minimal_errors(config)
            return nothing  # Validation failed with minimal error message
        end
    end
    
    return extract_values(config)
end

"""
    merge_configs(default_config::Dict, user_config::Dict)

Merge user config onto defaults, handling both simple values and parameter definitions.
"""
function merge_configs(default_config::Dict, user_config::Dict)
    result = deepcopy(default_config)
    
    for (section, values) in user_config
        if !haskey(result, section)
            result[section] = values
            continue
        end
        
        for (key, value) in values
            if !haskey(result[section], key)
                result[section][key] = value
                continue
            end
            
            # Handle parameter definitions
            if isa(result[section][key], Dict) && haskey(result[section][key], "value")
                result[section][key] = isa(value, Dict) ? 
                    merge(result[section][key], value) : 
                    merge(result[section][key], Dict("value" => value))
            else
                result[section][key] = value
            end
        end
    end
    
    return result
end

"""
    extract_values(config::Dict)

Extract final values from config, handling parameter definitions at any nesting level
with support for both dot notation and nested tables.
"""
function extract_values(config::Dict)
    result = Dict()
    
    # First, identify any dot-notation parameters (like a.value, a.type)
    param_prefixes = Set()
    for key in keys(config)
        if endswith(string(key), ".value")
            prefix = replace(string(key), ".value" => "")
            push!(param_prefixes, prefix)
        end
    end
    
    # Process all keys
    for (key, value) in config
        key_str = string(key)
        
        # Check if this is part of a dot-notation parameter
        is_part_of_param = false
        for prefix in param_prefixes
            if startswith(key_str, "$(prefix).")
                is_part_of_param = true
                break
            end
        end
        
        if is_part_of_param
            # Skip non-value parameter parts - they'll be handled together
            if !endswith(key_str, ".value")
                continue
            end
            
            # Get the parameter name without .value
            param_name = replace(key_str, ".value" => "")
            result[Symbol(param_name)] = value
            
        elseif isa(value, Dict)
            if haskey(value, "value")
                # This is a standard parameter definition
                result[key] = value["value"]
            else
                # This is a nested section
                result[key] = extract_values(value)
            end
        else
            # This is a simple value
            result[key] = value
        end
    end
    return result
end

"""
    validate_config(config::Dict, path="")

Validate config values against their definitions at any nesting level.
Supports both dot notation and nested table parameter definitions.
"""
function validate_config(config::Dict, path="")
    # First, identify any dot-notation parameters (like a.value, a.type)
    param_defs = Dict()
    
    for key in keys(config)
        key_str = string(key)
        if endswith(key_str, ".value")
            # Found a parameter definition using dot notation
            param_name = replace(key_str, ".value" => "")
            param_defs[param_name] = Dict("value" => config[key])
            
            # Look for other parameter properties with the same prefix
            for (prop_key, prop_value) in config
                prop_key_str = string(prop_key)
                if startswith(prop_key_str, "$(param_name).") && prop_key_str != "$(param_name).value"
                    prop_name = replace(prop_key_str, "$(param_name)." => "")
                    param_defs[param_name][prop_name] = prop_value
                end
            end
        end
    end
    
    # Process all entries for validation
    for (key, value) in config
        key_str = string(key)
        current_path = isempty(path) ? key_str : "$path.$key_str"
        
        # Skip keys that are part of dot-notation parameters
        skip = false
        for param_name in keys(param_defs)
            if startswith(key_str, "$(param_name).") 
                skip = true
                break
            end
        end
        if skip
            continue
        end
        
        if isa(value, Dict)
            if haskey(value, "type") && haskey(value, "value")
                # Regular parameter definition, validate it
                _validate_parameter(value, current_path)
            else
                # Nested section
                validate_config(value, current_path)
            end
        end
    end
    
    # Now validate any dot-notation parameters we found
    for (param_name, param_def) in param_defs
        current_path = isempty(path) ? param_name : "$path.$param_name"
        if haskey(param_def, "type")
            _validate_parameter(param_def, current_path)
        end
    end
    
    return true
end

"""
    _validate_parameter(param_def::Dict, path::String)

Helper function to validate a single parameter definition.
"""
function _validate_parameter(param_def::Dict, path::String)
    param_value = param_def["value"]
    param_type = param_def["type"]
   
    if param_type == "number"
        isa(param_value, Number) || 
            @minimal_error "$path must be a number, got $(typeof(param_value))"
        
        haskey(param_def, "min") && param_value < param_def["min"] && 
            @minimal_error "$path ($param_value) must be >= $(param_def["min"])"
        
        haskey(param_def, "max") && param_value > param_def["max"] && 
            @minimal_error "$path ($param_value) must be <= $(param_def["max"])"
    
    elseif param_type == "string"
        isa(param_value, String) || 
            @minimal_error "$path must be a string, got $(typeof(param_value))"
        
        if haskey(param_def, "values")
            param_value in param_def["values"] || 
                @minimal_error "$path ($param_value) must be one of: $(join(param_def["values"], ", "))"
        end
    end
    
    return true
end
