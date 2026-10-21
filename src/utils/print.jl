"""
    _print_vector(v::AbstractVector; max_length::Int = 10, n_ends::Int = 5)

Print a vector with a maximum length and number of ends.

# Arguments
- `v::AbstractVector`: The vector to print
- `max_length::Int`: The maximum length of the vector to print
- `n_ends::Int`: The number of ends to print
"""
function _print_vector(v::AbstractVector; max_length::Int = 10, n_ends::Int = 5)
    isempty(v) && return "[]"
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join(string.(v), ", ")
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
    config_with_meta["metadata"] =
        OrderedDict("julia_version" => string(VERSION), "EegFun_version" => string(pkgversion(@__MODULE__)), "date" => string(now()))

    # Copy config content, ensuring string keys and skipping any existing metadata
    for (key, value) in config
        if string(key) != "metadata"  # Skip any existing metadata
            config_with_meta[string(key)] = value
        end
    end
    TOML.print(x -> x isa Symbol ? string(x) : x, io, config_with_meta)
end

# Convenience method for printing to a file
function print_config(config, filename::String)
    open(filename, "w") do file
        print_config(config, file)
    end
    @info "Configuration written to: $filename"
end

"""
    version_info() -> Dict{String,Any}

Get version information for logging purposes.

Returns a dictionary with the keys:
- `"julia_version"`: Julia version string
- `"EegFun_version"`: EegFun package version string
- `"timestamp"`: Current timestamp string

# Examples
```julia
ver_info = version_info()
@info "Starting analysis" ver_info...

@info "Running EegFun \$(ver_info["EegFun_version"]) on Julia \$(ver_info["julia_version"])"
```
"""
function version_info()
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "EegFun_version" => string(pkgversion(@__MODULE__)),
        "timestamp" => string(now()),
    )
end

"""
    log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

Log a pretty table with specified log level. For general DataFrame logging.
Sets show_row_number=false and show_subheader=false by default for cleaner logs.

# Arguments
- `df::DataFrame`: The DataFrame to log
- `log_level::Symbol`: Log level (:debug, :info, :warn, :error) (default: :info)
- `kwargs...`: Additional arguments passed to pretty_table

# Examples
```julia
# Log with default info level
log_pretty_table(df; title = "My Table")

# Log with debug level
log_pretty_table(df; log_level = :debug, title = "Debug Table")

# Log with warn level
log_pretty_table(df; log_level = :warn, title = "Warning Table")
```
"""
function log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

    table_output = sprint() do output_io
        # Set a large display size to avoid terminal limitations (i.e., cropping)
        # when capturing the table string output.
        io_context = IOContext(output_io, :displaysize => (2000, 2000), :color => true)
        # Default column_labels to column names only — this suppresses the auto-generated
        # type subheader row that PrettyTables adds when given a DataFrame directly.
        # Callers can still override by passing their own column_labels kwarg.
        kw = Dict{Symbol,Any}(kwargs)
        get!(kw, :column_labels, names(df))
        pretty_table(io_context, df; kw...)
    end

    # Map symbol to standard LogLevel
    levels = Dict(:debug => Logging.Debug, :info => Logging.Info, :warn => Logging.Warn, :error => Logging.Error)

    lvl = get(levels, log_level, Logging.Info)
    if !haskey(levels, log_level)
        @minimal_warning "Unknown log level: $log_level, using :info instead"
    end

    # Log with minimal metadata (no file/line noise)
    @logmsg lvl "\n\n$table_output\n" _module=nothing _file=nothing _line=nothing

    return nothing
end
