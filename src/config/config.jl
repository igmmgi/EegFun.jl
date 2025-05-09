function load_config(config_file::String; validate::Bool=true)
    # Path to default configuration
    default_config_path = joinpath(@__DIR__, "default.toml")
    
    # Load default configuration first
    if !isfile(default_config_path)
        @warn "Default configuration not found at $default_config_path. Using empty defaults."
        default_config = Dict()
    else
        default_config = TOML.parsefile(default_config_path)
    end
    
    # Check if user config file exists
    if !isfile(config_file)
        throw(SystemError("Configuration file not found: $config_file"))
    end
    
    # Load user configuration
    user_config = Dict()
    try
        println("Loading config file: $config_file")
        user_config = TOML.parsefile(config_file)
    catch e
        throw(ArgumentError("Error parsing TOML file: $(e.msg)"))
    end
    
    # Smart merge user config onto defaults
    config = smart_merge(default_config, user_config)
    
    # Skip validation if not requested
    if !validate
        return finalize_config(config)
    end
    
    # Validate the configuration
    validate_config(config)
    
    # Return finalized configuration
    return finalize_config(config)
end

"""
    smart_merge(default_config::Dict, user_config::Dict)

Merge default and user configs, handling both simple values and complex parameter definitions.
"""
function smart_merge(default_config::Dict, user_config::Dict)
    result = deepcopy(default_config)
    
    for (section, section_values) in user_config
        if !haskey(result, section)
            result[section] = section_values
            continue
        end
        
        for (key, user_value) in section_values
            if !haskey(result[section], key)
                result[section][key] = user_value
                continue
            end
            
            # Handle the case where default is a parameter definition
            if isa(result[section][key], Dict) && haskey(result[section][key], "value")
                # If user provided a simple value, update the "value" field
                if !isa(user_value, Dict)
                    result[section][key]["value"] = user_value
                else
                    # If user provided a dict, do a normal merge
                    for (subkey, subvalue) in user_value
                        result[section][key][subkey] = subvalue
                    end
                end
            else
                # Simple value in default, just replace it
                result[section][key] = user_value
            end
        end
    end
    
    return result
end

"""
    finalize_config(config::Dict)

Extract just the values from a config with parameter definitions.
"""
function finalize_config(config::Dict)
    result = Dict()
    
    for (section, section_values) in config
        result[section] = Dict()
        
        for (key, value) in section_values
            if isa(value, Dict) && haskey(value, "value")
                result[section][key] = value["value"]
            else
                result[section][key] = value
            end
        end
    end
    
    return result
end

"""
    validate_config(config::Dict)

Validate config values against their definitions.
"""
function validate_config(config::Dict)
    for (section, section_values) in config
        for (key, value) in section_values
            if isa(value, Dict) && haskey(value, "type")
                # This is a parameter definition with validation info
                param_def = value
                param_value = param_def["value"]
                
                # Validate based on type
                if param_def["type"] == "float"
                    if !isa(param_value, Real)
                        throw(ArgumentError("$section.$key must be a number, got $(typeof(param_value))"))
                    end
                    
                    if haskey(param_def, "min") && param_value < param_def["min"]
                        throw(ArgumentError("$section.$key ($(param_value)) must be >= $(param_def["min"])"))
                    end
                    
                    if haskey(param_def, "max") && param_value > param_def["max"]
                        throw(ArgumentError("$section.$key ($(param_value)) must be <= $(param_def["max"])"))
                    end
                elseif param_def["type"] == "float_or_null" && haskey(param_def, "values")
                    if !(param_value in param_def["values"])
                        throw(ArgumentError("$section.$key ($(param_value)) must be one of: $(join(param_def["values"], ", "))"))
                    end
                end
                # Add more type validations as needed
            end
        end
    end
    
    return true
end
