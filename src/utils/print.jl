function print_vector_(v::Vector; max_length::Int = 10, n_ends::Int = 5)
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join(v, ", ")
end

function print_vector(v::UnitRange; max_length::Int = 10, n_ends::Int = 5)
    print_vector_(collect(v), max_length = max_length, n_ends = n_ends)
end

function print_vector(v::Vector; max_length::Int = 10, n_ends::Int = 5)
    print_vector_(collect(v), max_length = max_length, n_ends = n_ends)
end

"""
    get_git_commit()

Get the current git commit hash.

# Returns
- `String`: The full git commit hash, or "unknown" if not available
"""
function get_git_commit()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

"""
    get_eegfun_version()

Get the EEGfun version from Project.toml.

# Returns
- `String`: The version string from Project.toml, or "unknown" if not available
"""
function get_eegfun_version()
    try
        project_toml = TOML.parsefile("Project.toml")
        return project_toml["version"]
    catch
        return "unknown"
    end
end

"""
    print_config(config, [io=stdout])
    print_config(config, filename::String)

Print configuration in TOML format.

# Arguments
- `config`: Configuration dictionary (typically loaded from TOML)
- `io=stdout`: Optional IO object to print to (default: standard output)  
- `filename`: Optional filename to write TOML output to

# Examples
```julia
# Print to console in TOML format
print_config(config)

# Write to TOML file
print_config(config, "config_output.toml")
```
"""
function print_config(config, io::IO=stdout)
    # Use OrderedDict to ensure metadata appears first
    config_with_meta = OrderedDict{String, Any}()
    
    # Always add fresh metadata first
    config_with_meta["metadata"] = OrderedDict(
        "generated_at" => string(now()),
        "eegfun_version" => get_eegfun_version(),
        "git_commit" => get_git_commit()
    )
    
    # Copy config content, ensuring string keys and skipping any existing metadata
    for (key, value) in config
        if string(key) != "metadata"  # Skip any existing metadata
            config_with_meta[string(key)] = value
        end
    end
    
    TOML.print(io, config_with_meta)
end

# Convenience method for printing to a file
function print_config(config, filename::String)
    open(filename, "w") do file
        print_config(config, file)
    end
    println("Configuration written to: $filename")
end

