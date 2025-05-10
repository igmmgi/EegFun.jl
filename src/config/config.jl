"""
    load_config(config_file::String; validate::Bool=true)

Load and merge configuration from a TOML file with defaults.
"""
function load_config(config_file::String; validate::Bool=true)
    # Load default config
    default_config = TOML.parsefile(joinpath(@__DIR__, "default.toml"))
    
    # Check if user config exists
    if !isfile(config_file)
        throw(SystemError("Configuration file not found: $config_file"))
    end
    
    # Load user config
    user_config = Dict()
    try
        @info "Loading config file: $config_file"
        user_config = TOML.parsefile(config_file)
    catch e
        throw(ArgumentError("Error parsing TOML file: $(e.msg)"))
    end
    
    # Merge and validate
    config = merge_configs(default_config, user_config)
    validate && validate_config(config)
    
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

Extract final values from config, handling parameter definitions.
"""
function extract_values(config::Dict)
    result = Dict()
    for (section, values) in config
        result[section] = Dict()
        for (key, value) in values
            result[section][key] = isa(value, Dict) && haskey(value, "value") ? value["value"] : value
        end
    end
    return result
end

"""
    validate_config(config::Dict)

Validate config values against their definitions.
"""
function validate_config(config::Dict)
    for (section, values) in config
        for (key, value) in values
            isa(value, Dict) && haskey(value, "type") || continue
            
            param_value = value["value"]
            param_type = value["type"]
            
            if param_type == "number"
                isa(param_value, Number) || 
                    throw(ArgumentError("$section.$key must be a number, got $(typeof(param_value))"))
                
                haskey(value, "min") && param_value < value["min"] && 
                    throw(ArgumentError("$section.$key ($param_value) must be >= $(value["min"])"))
                
                haskey(value, "max") && param_value > value["max"] && 
                    throw(ArgumentError("$section.$key ($param_value) must be <= $(value["max"])"))
            
            elseif param_type == "string"
                isa(param_value, String) || 
                    throw(ArgumentError("$section.$key must be a string, got $(typeof(param_value))"))
                
                if haskey(value, "values")
                    param_value in value["values"] || 
                        throw(ArgumentError("$section.$key ($param_value) must be one of: $(join(value["values"], ", "))"))
                end
            end
        end
    end
    
    return true
end
