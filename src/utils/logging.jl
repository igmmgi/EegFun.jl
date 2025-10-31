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

# Internal: normalize log level
_to_loglevel(level::Symbol) =
    level === :debug ? Logging.Debug :
    level === :info  ? Logging.Info  :
    level === :warn  ? Logging.Warn  :
    level === :error ? Logging.Error : Logging.Info

# Back-compat overload for String input
_to_loglevel(level::String) = _to_loglevel(Symbol(lowercase(level)))

function setup_global_logging(log_file::String; log_level::Symbol = :info)

    # Convert log level to LogLevel and store it
    level = _to_loglevel(log_level)
    LOG_STATE.log_level = level

    # Close any existing global log file
    if !isnothing(LOG_STATE.global_log_handle)
        close(LOG_STATE.global_log_handle)
    end

    # Open new global log file
    LOG_STATE.global_log_handle = open(log_file, "w")

    # Record start time
    LOG_STATE.global_log_start_time = now()

    # Some version info
    version_info = eegfun_version_info()
    println(LOG_STATE.global_log_handle, "julia version: $(version_info["julia_version"])")
    println(LOG_STATE.global_log_handle, "eegfun version: $(version_info["eegfun_version"])")
    println(LOG_STATE.global_log_handle, "Start time: ", Dates.format(LOG_STATE.global_log_start_time, "dd-mm-yyyy HH:MM:SS"))

    # Set up the global logger, but only if we're not already file logging
    if isnothing(LOG_STATE.log_handle)
        # Create a logger that writes to both stdout and file
        console_logger = ConsoleLogger(stdout, level)
        file_logger = FormatLogger(LOG_STATE.global_log_handle) do io, log
            if log.level >= level
                println(io, log.message)
            end
        end

        # Save current logger before replacing
        LOG_STATE.saved_logger = global_logger()

        # Combine both loggers and set as global
        logger = TeeLogger(console_logger, file_logger)
        global_logger(logger)
    end

    return LOG_STATE.global_log_handle
end

"""
    close_global_logging()

Close the global log file and add a closing timestamp with duration.
"""
function close_global_logging()
    if !isnothing(LOG_STATE.global_log_handle)
        duration = now() - LOG_STATE.global_log_start_time
        println(LOG_STATE.global_log_handle, "\nDuration: ", format_duration(duration))
        close(LOG_STATE.global_log_handle)
        LOG_STATE.global_log_handle = nothing
        LOG_STATE.global_log_start_time = nothing

        # Only restore logger if we're not in file logging mode
        if isnothing(LOG_STATE.log_handle)
            global_logger(LOG_STATE.saved_logger)
        end
    end
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
    # Convert log level to LogLevel
    level = _to_loglevel(log_level)

    # Store the current log level
    LOG_STATE.log_level = level

    # Check if we're entering logging for the first time (need to save logger)
    was_logging = !isnothing(LOG_STATE.log_handle)

    # Close any existing log file
    if was_logging
        close(LOG_STATE.log_handle)
    end

    # Open new log file (creates a new file each time)
    LOG_STATE.log_handle = open(log_file, "w")

    # Record start time
    LOG_STATE.log_start_time = now()

    # Some version info
    version_info = eegfun_version_info()
    println(LOG_STATE.log_handle, "julia version: $(version_info["julia_version"])")
    println(LOG_STATE.log_handle, "eegfun version: $(version_info["eegfun_version"])")
    println(LOG_STATE.log_handle, "Start time: ", Dates.format(LOG_STATE.log_start_time, "dd-mm-yyyy HH:MM:SS"))

    # Only save logger if we're entering logging for the first time
    if !was_logging
        LOG_STATE.saved_logger = global_logger()
    end

    # Create a logger that writes to both stdout and file
    console_logger = ConsoleLogger(stdout, level)
    file_logger = FormatLogger(LOG_STATE.log_handle) do io, log
        if log.level >= level
            println(io, log.message)
            if !isempty(log.kwargs)
                for (key, val) in log.kwargs
                    println(io, "$key = $val")
                end
            end
        end
    end

    # Combine both loggers
    logger = TeeLogger(console_logger, file_logger)
    global_logger(logger)

    return nothing
end

"""
    close_logging()

Close the individual log file and add a closing timestamp with duration.
Restores the previous logger (which may be the global logger).
"""
function close_logging()
    if !isnothing(LOG_STATE.log_handle)
        duration = now() - LOG_STATE.log_start_time
        println(LOG_STATE.log_handle, "\nFile Processing Duration: ", format_duration(duration))
        close(LOG_STATE.log_handle)
        LOG_STATE.log_handle = nothing
        LOG_STATE.log_start_time = nothing

        # Restore the saved logger (which was the global logger before file logging replaced it)
        global_logger(LOG_STATE.saved_logger)
    end
end
