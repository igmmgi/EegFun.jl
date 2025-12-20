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
    ConsoleLogger(stdout), # saved_logger
    Logging.Info,          # log_level
)

"""
    format_duration(duration::Millisecond)

Format a duration in a human-readable way, choosing appropriate units.
"""
function format_duration(duration::Millisecond)
    total_seconds = Int(round(duration.value / 1000))

    if total_seconds < 60
        return "$total_seconds seconds"
    elseif total_seconds < 3600
        minutes = div(total_seconds, 60)
        seconds = rem(total_seconds, 60)
        return "$minutes minutes, $seconds seconds"
    else
        hours = div(total_seconds, 3600)
        minutes = div(rem(total_seconds, 3600), 60)
        seconds = rem(total_seconds, 60)
        return "$hours hours, $minutes minutes, $seconds seconds"
    end
end

# ============================================================================
# Helper Functions
# ============================================================================

# Internal: normalize log level
_to_loglevel(level::Symbol) =
    level === :debug ? Logging.Debug :
    level === :info  ? Logging.Info  :
    level === :warn  ? Logging.Warn  :
    level === :error ? Logging.Error : Logging.Info

# Back-compat overload for String input
_to_loglevel(level::String) = _to_loglevel(Symbol(lowercase(level)))

"""
    _write_log_header(io::IO, start_time::DateTime)

Write version info and start time to a log file.
"""
function _write_log_header(io::IO, start_time::DateTime)
    version_info = eegfun_version_info()
    println(io, "julia version: $(version_info["julia_version"])")
    println(io, "eegfun version: $(version_info["eegfun_version"])")
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
function _create_tee_logger(
    console_level::Logging.LogLevel,
    file_handle::IO,
    file_level::Logging.LogLevel;
    include_kwargs::Bool = false,
)
    console_logger = ConsoleLogger(stdout, console_level)
    file_logger = _create_file_logger(file_handle, file_level; include_kwargs = include_kwargs)
    return TeeLogger(console_logger, file_logger)
end

"""
    setup_global_logging(log_file::String)

Set up a global log file that will persist throughout the session. This is meant to be used
for overall program logging, while individual file processing can use setup_logging.

# Arguments
- `log_file::String`: Path to the global log file

# Returns
- The opened file handle for direct writing (using write/println)

# Example
```julia
global_log = setup_global_logging("global_analysis.log")
@info "This is logged to the global log and stdout"
write(global_log, "Direct write to global log\\n")
```
"""

function setup_global_logging(log_file::String; log_level::Symbol = :info)
    # Convert log level and store it
    level = _to_loglevel(log_level)
    LOG_STATE.log_level = level

    # Open log file (handles closing existing, recording time, writing header)
    handle, start_time = _open_log_file(LOG_STATE.global_log_handle, LOG_STATE.global_log_start_time, log_file)
    LOG_STATE.global_log_handle = handle
    LOG_STATE.global_log_start_time = start_time

    # Set up the global logger, but only if we're not already file logging
    if isnothing(LOG_STATE.log_handle)
        # Save current logger before replacing
        LOG_STATE.saved_logger = global_logger()

        # Create and set the combined logger
        logger = _create_tee_logger(level, handle, level)
        global_logger(logger)
    end

    return handle
end

"""
    close_global_logging()

Close the global log file and add a closing timestamp with duration.
"""
function close_global_logging()
    _close_log_file(LOG_STATE.global_log_handle, LOG_STATE.global_log_start_time, "Duration")
    LOG_STATE.global_log_handle = nothing
    LOG_STATE.global_log_start_time = nothing

    # Only restore logger if we're not in file logging mode
    isnothing(LOG_STATE.log_handle) && global_logger(LOG_STATE.saved_logger)

end

"""
    setup_logging(log_file::String)

Set up logging to both stdout and a file for processing an individual file.
This temporarily replaces any global logger and restores it when close_logging is called.

# Arguments
- `log_file::String`: Path to the log file

# Returns
- Nothing (use standard logging macros: @info, @warn, @error)

# Example
```julia
setup_logging("my_analysis.log")
@info "This will be logged to both stdout and the file"
```
"""
function setup_logging(log_file::String; log_level::Symbol = :info)
    # Convert log level and store it
    level = _to_loglevel(log_level)
    LOG_STATE.log_level = level

    # Check if we're entering logging for the first time (need to save logger)
    was_logging = !isnothing(LOG_STATE.log_handle)

    # Open log file (handles closing existing, recording time, writing header)
    handle, start_time = _open_log_file(LOG_STATE.log_handle, LOG_STATE.log_start_time, log_file)
    LOG_STATE.log_handle = handle
    LOG_STATE.log_start_time = start_time

    # Only save logger if we're entering logging for the first time
    !was_logging && (LOG_STATE.saved_logger = global_logger())

    # Create and set the combined logger (with kwargs support for file logging)
    logger = _create_tee_logger(level, handle, level; include_kwargs = true)
    global_logger(logger)

    return nothing
end

"""
    close_logging()

Close the individual log file and add a closing timestamp with duration.
Restores the previous logger (which may be the global logger).
"""
function close_logging()
    _close_log_file(LOG_STATE.log_handle, LOG_STATE.log_start_time, "File Processing Duration")
    LOG_STATE.log_handle = nothing
    LOG_STATE.log_start_time = nothing

    # Restore the saved logger (which was the global logger before file logging replaced it)
    global_logger(LOG_STATE.saved_logger)
end


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

# Notes
- Simply reads the last line from Julia's REPL history file
- Logs the actual command as entered
"""
macro log_call(func_name)
    func_name isa String || error("@log_call: argument must be a string (function name)")
    return quote
        last_line = _get_last_history_line()
    end
end
