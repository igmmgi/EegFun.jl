function _print_vector(v::AbstractVector; max_length::Int = 10, n_ends::Int = 5)
    isempty(v) && return "[]"
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join((string(x) for x in v), ", ")
end

"""
    print_config(config, [io::IO = stdout])
    print_config(config, filename::String)

Write `config` in TOML format (with fresh `metadata` prepended) to `io` or `filename`.
"""
function print_config(config, io::IO = stdout)
    config_with_meta = OrderedDict{String,Any}()
    config_with_meta["metadata"] =
        OrderedDict("julia_version" => string(VERSION), "EegFun_version" => string(pkgversion(@__MODULE__)), "date" => string(now()))
    for (key, value) in config
        sk = string(key)
        sk != "metadata" && (config_with_meta[sk] = value)
    end
    TOML.print(x -> x isa Symbol ? string(x) : x, io, config_with_meta)
end

function print_config(config, filename::String)
    open(filename, "w") do file
        print_config(config, file)
    end
    @info "Configuration written to: $filename"
end

"""
    version_info() -> Dict{String,Any}

Return a dictionary containing `"julia_version"`, `"EegFun_version"`, and `"timestamp"`.
"""
function version_info()
    return Dict{String,Any}(
        "julia_version" => string(VERSION),
        "EegFun_version" => string(pkgversion(@__MODULE__)),
        "timestamp" => string(now()),
    )
end

function log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)
    table_output = sprint() do output_io
        io_context = IOContext(output_io, :displaysize => (2000, 2000), :color => true)
        kw = Dict{Symbol,Any}(kwargs)
        get!(kw, :column_labels, names(df))
        pretty_table(io_context, df; kw...)
    end

    if !haskey(LOG_LEVEL_MAP, log_level)
        @minimal_warning "Unknown log level: $log_level, using :info instead"
    end
    lvl = get(LOG_LEVEL_MAP, log_level, Logging.Info)
    @logmsg lvl "\n\n$table_output\n" _module = nothing _file = nothing _line = nothing
    return nothing
end
