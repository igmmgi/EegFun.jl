"""
    print_vector(v::AbstractVector; max_length::Int = 10, n_ends::Int = 5)

Print a vector with a maximum length and number of ends.

# Arguments
- `v::AbstractVector`: The vector to print
- `max_length::Int`: The maximum length of the vector to print
- `n_ends::Int`: The number of ends to print
"""
function print_vector(v::AbstractVector; max_length::Int = 10, n_ends::Int = 5)
    isempty(v) && return "[]"
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join(string.(v), ", ")
end


"""
    get_package_version(; package_name::String = "EegFun")

Get the package version from Project.toml.

It locates the package using Base.find_package() and reads Project.toml
from the package root.

# Returns
- `String`: The version string from Project.toml, or "unknown" if not available/problem
"""
function get_package_version(; package_name::String = "EegFun")
    try
        pkg_path = Base.find_package(package_name)
        version = "unknown"
        if pkg_path !== nothing
            package_root = dirname(dirname(pkg_path))
            project_file = joinpath(package_root, "Project.toml")
            if isfile(project_file)
                project_toml = TOML.parsefile(project_file)
                if haskey(project_toml, "version")
                    version = project_toml["version"]
                end
            end
        end
        return version
    catch # TODO: can this actually fail if installed properly?
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
function print_config(config, io::IO = stdout)
    # Use OrderedDict to ensure package data appears first
    config_with_meta = OrderedDict{String,Any}()

    # Always add fresh metadata first
    config_with_meta["metadata"] = OrderedDict(
        "julia_version" => string(VERSION),
        "EegFun_version" => get_package_version(package_name = "EegFun"),
        "date" => string(now()),
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
    @info "Configuration written to: $filename"
end

"""
    EegFun_version_info() -> Dict{String,Any}

Get comprehensive version information for logging purposes.

This function returns a dictionary containing:
- `julia_version`: Julia version
- `EegFun_version`: EEGfun package version
- `git_commit`: Git commit hash (if available)
- `timestamp`: Current timestamp

# Returns
- `Dict{String,Any}`: Dictionary with version information

# Examples
```julia
# Get version info for logging
ver_info = EegFun_version_info()
@info "Starting analysis" ver_info...

# Log specific fields
ver_info = EegFun_version_info()
@info "Running EegFun \$(ver_info["EegFun_version"]) on Julia \$(ver_info["julia_version"])"

# Write to log file
ver_info = EegFun_version_info()
open("analysis.log", "a") do io
    println(io, "Version: \$(ver_info["EegFun_version"]), Commit: \$(ver_info["git_commit"])")
end
```
"""
function EegFun_version_info()
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "EegFun_version" => get_package_version(package_name = "EegFun"),
        "timestamp" => string(now()),
    )
end

