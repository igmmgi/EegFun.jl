mutable struct LoggingState
    log_handle::Union{Nothing,IO}
    log_start_time::Union{Nothing,DateTime}
    log_include_kwargs::Bool
    global_log_handle::Union{Nothing,IO}
    global_log_start_time::Union{Nothing,DateTime}
    global_log_include_kwargs::Bool
    saved_logger::Any
    log_level::Logging.LogLevel
end

const LOG_LOCK = ReentrantLock()

const LOG_STATE = LoggingState(nothing, nothing, true, nothing, nothing, false, nothing, Logging.Info)

_file_logging_active() = !isnothing(LOG_STATE.log_handle) || !isnothing(LOG_STATE.global_log_handle)

function _format_duration(duration::Millisecond)
    s = round(Int, duration.value / 1000)
    s == 0 && return "0 seconds"
    return string(Dates.canonicalize(Second(s)))
end

const LOG_LEVEL_MAP = Dict{Symbol,Logging.LogLevel}(
    :debug => Logging.Debug,
    :info => Logging.Info,
    :warn => Logging.Warn,
    :warning => Logging.Warn,
    :error => Logging.Error,
)

_to_loglevel(level::Logging.LogLevel) = level
_to_loglevel(level::Symbol) = get(LOG_LEVEL_MAP, level, Logging.Info)
_to_loglevel(level::String) = _to_loglevel(Symbol(lowercase(level)))

function _write_log_header(io::IO, start_time::DateTime)
    ver_info = version_info()
    println(io, "julia version: $(ver_info["julia_version"])")
    println(io, "EegFun version: $(ver_info["EegFun_version"])")
    println(io, "Start time: ", Dates.format(start_time, "dd-mm-yyyy HH:MM:SS"))
end

function _open_log_file(handle_field::Union{Nothing,IO}, filename::String)
    !isnothing(handle_field) && close(handle_field)
    handle = open(filename, "w")
    start_time = now()
    _write_log_header(handle, start_time)
    return handle, start_time
end

function _close_log_file(handle::Union{Nothing,IO}, start_time::Union{Nothing,DateTime}, label::String)
    if !isnothing(handle)
        println(handle, "\n$(label): ", _format_duration(now() - start_time))
        close(handle)
    end
end

function _create_file_logger(io::IO, level::Logging.LogLevel; include_kwargs::Bool = false)
    return FormatLogger(io) do io_inner, log
        lock(LOG_LOCK) do
            if isopen(io_inner) && log.level >= level
                println(io_inner, log.message)
                if include_kwargs && !isempty(log.kwargs)
                    for (key, val) in log.kwargs
                        println(io_inner, "$key = $val")
                    end
                end
                flush(io_inner)
            end
        end
    end
end

function _create_tee_logger(console_level::Logging.LogLevel, file_handle::IO, file_level::Logging.LogLevel; include_kwargs::Bool = false)
    return TeeLogger(ConsoleLogger(stdout, console_level), _create_file_logger(file_handle, file_level; include_kwargs = include_kwargs))
end

function _apply_active_logger!()
    use_local = !isnothing(LOG_STATE.log_handle)
    active_handle = use_local ? LOG_STATE.log_handle : LOG_STATE.global_log_handle
    active_include = use_local ? LOG_STATE.log_include_kwargs : LOG_STATE.global_log_include_kwargs
    level = LOG_STATE.log_level
    global_logger(_create_tee_logger(level, active_handle, level; include_kwargs = active_include))
end

"""
    setup_logging(log_file::String; log_level::Symbol = :info, is_global::Bool = false, include_kwargs::Bool = !is_global)
    setup_logging(f::Function, log_file::String; is_global::Bool = false, kwargs...)

Configure logging to both `stdout` and `log_file`. Local file logs take precedence over global logs.
"""
function setup_logging(log_file::String; log_level::Symbol = :info, is_global::Bool = false, include_kwargs::Bool = !is_global)
    return lock(LOG_LOCK) do
        LOG_STATE.log_level = _to_loglevel(log_level)
        if isnothing(LOG_STATE.saved_logger) && !(global_logger() isa TeeLogger)
            LOG_STATE.saved_logger = global_logger()
        end

        if is_global
            handle, start_time = _open_log_file(LOG_STATE.global_log_handle, log_file)
            LOG_STATE.global_log_handle = handle
            LOG_STATE.global_log_start_time = start_time
            LOG_STATE.global_log_include_kwargs = include_kwargs
        else
            handle, start_time = _open_log_file(LOG_STATE.log_handle, log_file)
            LOG_STATE.log_handle = handle
            LOG_STATE.log_start_time = start_time
            LOG_STATE.log_include_kwargs = include_kwargs
        end

        _apply_active_logger!()
        return handle
    end
end

function setup_logging(f::Function, log_file::String; is_global::Bool = false, kwargs...)
    handle = setup_logging(log_file; is_global = is_global, kwargs...)
    try
        return applicable(f, handle) ? f(handle) : f()
    finally
        close_logging(; is_global = is_global)
    end
end

"""
    setup_global_logging(log_file::String; kwargs...)
    setup_global_logging(f::Function, log_file::String; kwargs...)

Configure session-wide file logging (equivalent to `setup_logging(...; is_global = true)`).
"""
setup_global_logging(args...; kwargs...) = setup_logging(args...; is_global = true, kwargs...)

"""
    close_logging(; is_global::Bool = false)

Close the active local (or global) log file and restore the previous logger.
"""
function close_logging(; is_global::Bool = false)
    lock(LOG_LOCK) do
        if is_global
            _close_log_file(LOG_STATE.global_log_handle, LOG_STATE.global_log_start_time, "Global session duration")
            LOG_STATE.global_log_handle = nothing
            LOG_STATE.global_log_start_time = nothing
        else
            _close_log_file(LOG_STATE.log_handle, LOG_STATE.log_start_time, "File processing duration")
            LOG_STATE.log_handle = nothing
            LOG_STATE.log_start_time = nothing
        end

        if isnothing(LOG_STATE.log_handle) && isnothing(LOG_STATE.global_log_handle)
            if !isnothing(LOG_STATE.saved_logger)
                global_logger(LOG_STATE.saved_logger)
                LOG_STATE.saved_logger = nothing
            else
                global_logger(ConsoleLogger(stdout, LOG_STATE.log_level))
            end
        else
            _apply_active_logger!()
        end
    end
    return nothing
end

"""
    close_global_logging()

Close the active global session log file and restore the previous logger.
"""
close_global_logging() = close_logging(is_global = true)

function _format_log_call(name_str::String, local_vars::Dict{Symbol,Any})
    arg_strs = String[]
    if haskey(local_vars, :dat)
        push!(arg_strs, "dat")
    elseif haskey(local_vars, :data)
        push!(arg_strs, "data")
    end

    for (k, v) in local_vars
        if k === :dat || k === :data || k === :_eegfun_locals || k === Symbol("#self#") || startswith(string(k), "#") || v isa BatchContext
            continue
        end
        if v isa Union{Number,Symbol,AbstractString,Tuple} ||
           v isa AbstractVector{<:Number} ||
           v isa AbstractVector{Symbol} ||
           v isa AbstractVector{<:AbstractString}
            push!(arg_strs, "$k=$(repr(v))")
        elseif v isa Function
            push!(arg_strs, "$k=$(string(nameof(v)))")
        else
            push!(arg_strs, "$k=::$(typeof(v))")
        end
    end
    return string(name_str, "(", join(arg_strs, ", "), ")")
end

macro log_call(func_name)
    name_str = string(func_name)
    return quote
        local _eegfun_locals = Base.@locals
        @info "Function call: $(_format_log_call($name_str, _eegfun_locals))"
    end
end

_progress_stopped(::Nothing) = false
_progress_stopped(flag::Base.Threads.Atomic{Bool}) = flag[]
_progress_stopped(flag::Base.RefValue{Bool}) = flag[]

function track_progress_async(
    prog,
    counter::Base.Threads.Atomic{Int},
    total::Int;
    stop_flag::Union{Nothing,Base.Threads.Atomic{Bool},Base.RefValue{Bool}} = nothing,
)
    @async begin
        while counter[] < total && !_progress_stopped(stop_flag)
            ProgressMeter.update!(prog, counter[])
            sleep(0.5)
        end
        !_progress_stopped(stop_flag) && ProgressMeter.update!(prog, total)
    end
end

function track_progress_async(f::Function, prog, counter::Base.Threads.Atomic{Int}, total::Int)
    stop_flag = Base.Threads.Atomic{Bool}(false)
    task = track_progress_async(prog, counter, total; stop_flag = stop_flag)
    try
        return f()
    finally
        stop_flag[] = true
        try
            wait(task)
        catch
        end
    end
end
