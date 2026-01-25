# Keep track of logging state in a single mutable struct
mutable struct LoggingState
    log_handle::Union{Nothing,IO}
    log_start_time::Union{Nothing,DateTime}
    global_log_handle::Union{Nothing,IO}
    global_log_start_time::Union{Nothing,DateTime}
    saved_logger::Any
    log_level::Logging.LogLevel
end

# Initialize global logging state
const LOG_STATE = LoggingState(
    nothing,               # log_handle
    nothing,               # log_start_time
    nothing,               # global_log_handle
    nothing,               # global_log_start_time
    nothing,               # saved_logger (captured at runtime)
    Logging.Info,          # log_level
)

"""
    format_duration(duration::Millisecond)

Format a duration in a human-readable way, choosing appropriate units.
"""
function format_duration(duration::Millisecond)
    s = round(Int, duration.value / 1000) # nearest S
    s == 0 && return "0 seconds"
    return string(Dates.canonicalize(Second(s)))
end

# ============================================================================
# Helper Functions
# ============================================================================

# Internal: map symbols to log levels
const LOG_LEVEL_MAP = Dict{Symbol,Logging.LogLevel}(
    :debug => Logging.Debug,
    :info => Logging.Info,
    :warn => Logging.Warn,
    :warning => Logging.Warn,
    :error => Logging.Error,
)

# Internal: normalize log level
_to_loglevel(level::Logging.LogLevel) = level
_to_loglevel(level::Symbol) = get(LOG_LEVEL_MAP, level, Logging.Info)

# Back-compat overload for String input
_to_loglevel(level::String) = _to_loglevel(Symbol(lowercase(level)))

"""
    _write_log_header(io::IO, start_time::DateTime)

Write version info and start time to a log file.
"""
function _write_log_header(io::IO, start_time::DateTime)
    version_info = EegFun_version_info()
    println(io, "julia version: $(version_info["julia_version"])")
    println(io, "EegFun version: $(version_info["EegFun_version"])")
    println(io, "Start time: ", Dates.format(start_time, "dd-mm-yyyy HH:MM:SS"))
end

"""
    _open_log_file(handle_field::Union{Nothing,IO}, start_time_field::Union{Nothing,DateTime}, filename::String) -> Tuple{IO, DateTime}

Open a log file, closing any existing one, and record the start time.
Returns the opened file handle and start time.
"""
function _open_log_file(handle_field::Union{Nothing,IO}, start_time_field::Union{Nothing,DateTime}, filename::String)
    !isnothing(handle_field) && close(handle_field)

    handle = open(filename, "w")
    start_time = now()
    _write_log_header(handle, start_time)

    return handle, start_time
end

"""
    _close_log_file(handle::Union{Nothing,IO}, start_time::Union{Nothing,DateTime}, label::String)

Close a log file, write duration, and return nothing for cleanup.
"""
function _close_log_file(handle::Union{Nothing,IO}, start_time::Union{Nothing,DateTime}, label::String)
    if !isnothing(handle)
        duration = now() - start_time
        println(handle, "\n$(label): ", format_duration(duration))
        close(handle)
    end
end

"""
    _create_file_logger(io::IO, level::Logging.LogLevel; include_kwargs::Bool = false)

Create a FormatLogger that writes to the given IO handle.
"""
function _create_file_logger(io::IO, level::Logging.LogLevel; include_kwargs::Bool = false)
    return FormatLogger(io) do io_inner, log
        if log.level >= level
            println(io_inner, log.message)
            if include_kwargs && !isempty(log.kwargs)
                for (key, val) in log.kwargs
                    println(io_inner, "$key = $val")
                end
            end
        end
    end
end

"""
    _create_tee_logger(console_level::Logging.LogLevel, file_handle::IO, file_level::Logging.LogLevel; include_kwargs::Bool = false)

Create a TeeLogger that writes to both console and file.
"""
function _create_tee_logger(console_level::Logging.LogLevel, file_handle::IO, file_level::Logging.LogLevel; include_kwargs::Bool = false)
    console_logger = ConsoleLogger(stdout, console_level)
    file_logger = _create_file_logger(file_handle, file_level; include_kwargs = include_kwargs)
    return TeeLogger(console_logger, file_logger)
end

"""
    setup_logging(log_file::String; log_level::Symbol = :info, is_global::Bool = false, include_kwargs::Bool = !is_global)

Set up logging to both stdout and a file. 

# Arguments
- `log_file`: Path to the log file.
- `log_level`: Logging level (:debug, :info, :warn, :error).
- `is_global`: If true, sets up a session-wide log. If false (default), sets up for individual file processing.
- `include_kwargs`: Whether to include metadata in the log. Defaults to true for local logs.

# Returns
- The file handle (for direct writing if desired).
"""
function setup_logging(log_file::String; log_level::Symbol = :info, is_global::Bool = false, include_kwargs::Bool = !is_global)
    level = _to_loglevel(log_level)
    LOG_STATE.log_level = level

    h_field = is_global ? :global_log_handle : :log_handle
    t_field = is_global ? :global_log_start_time : :log_start_time

    # Capture original logger if we haven't yet (before replacing with our first file logger)
    if isnothing(LOG_STATE.log_handle) && isnothing(LOG_STATE.global_log_handle)
        LOG_STATE.saved_logger = global_logger()
    end

    # Open log file (closes existing if necessary)
    handle, start_time = _open_log_file(getproperty(LOG_STATE, h_field), getproperty(LOG_STATE, t_field), log_file)

    # Update state
    setproperty!(LOG_STATE, h_field, handle)
    setproperty!(LOG_STATE, t_field, start_time)

    # Re-apply logger.
    # If local is active, it always takes precedence. Otherwise global.
    active_handle = !isnothing(LOG_STATE.log_handle) ? LOG_STATE.log_handle : LOG_STATE.global_log_handle
    active_include = !isnothing(LOG_STATE.log_handle) ? true : false

    global_logger(_create_tee_logger(level, active_handle, level; include_kwargs = active_include))

    return handle
end

# Convenience aliases
setup_global_logging(args...; kwargs...) = setup_logging(args...; is_global = true, kwargs...)

"""
    close_logging(; is_global::Bool = false)

Close the currently active log file and restore previous logger.
"""
function close_logging(; is_global::Bool = false)
    h_field = is_global ? :global_log_handle : :log_handle
    t_field = is_global ? :global_log_start_time : :log_start_time
    label = is_global ? "Global session duration" : "File processing duration"

    _close_log_file(getproperty(LOG_STATE, h_field), getproperty(LOG_STATE, t_field), label)

    setproperty!(LOG_STATE, h_field, nothing)
    setproperty!(LOG_STATE, t_field, nothing)

    # If both file loggers are now closed, restore the original saved logger
    if isnothing(LOG_STATE.log_handle) && isnothing(LOG_STATE.global_log_handle)
        if !isnothing(LOG_STATE.saved_logger)
            global_logger(LOG_STATE.saved_logger)
            LOG_STATE.saved_logger = nothing
        end
    else
        # Still have one file logger active (must be the global one if were closing local, 
        # or we just closed global while local was active)
        active_handle = !isnothing(LOG_STATE.log_handle) ? LOG_STATE.log_handle : LOG_STATE.global_log_handle
        active_include = !isnothing(LOG_STATE.log_handle) ? true : false
        global_logger(_create_tee_logger(LOG_STATE.log_level, active_handle, LOG_STATE.log_level; include_kwargs = active_include))
    end
end

close_global_logging() = close_logging(is_global = true)


# ============================================================================
# FUNCTION CALL LOGGING
# ============================================================================
# TODO: Must be a better way of doing this. This hack is simpler than previous code mess but ...

"""Helper to get the last command from Julia REPL history file."""
function _get_last_history_line()::Union{String,Nothing}
    history_file = joinpath(homedir(), ".julia", "logs", "repl_history.jl")
    if isfile(history_file)
        try # Read last line (skip empty lines)
            lines = readlines(history_file)
            for i = length(lines):-1:1
                line = strip(lines[i])
                if !isempty(line) && !startswith(line, "#")
                    return line
                end
            end
        catch
            return nothing
        end
    end
    return nothing
end

"""
    @log_call func_name

Macro to log a function call by reading the last command from history file.
If history is unavailable, it logs the provided function name.

# Arguments
- `func_name`: String or Symbol representing the function name.
"""
macro log_call(func_name)
    # Ensure func_name is a string for logging
    name_str = string(func_name)
    return quote
        last_line = _get_last_history_line()
        if isnothing(last_line)
            @info "Function call: $($name_str)"
        else
            @info "Function call: $last_line"
        end
    end
end
