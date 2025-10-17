# Keep track of the open file handles and start times
const log_file_handle = Ref{Union{Nothing,IO}}(nothing)
const log_start_time = Ref{Union{Nothing,DateTime}}(nothing)
const global_log_handle = Ref{Union{Nothing,IO}}(nothing)
const global_log_start_time = Ref{Union{Nothing,DateTime}}(nothing)
const saved_logger = Ref{Any}(ConsoleLogger(stdout))
const is_file_logging = Ref{Bool}(false)
const is_global_logging = Ref{Bool}(false)
const current_log_level = Ref{Logging.LogLevel}(Logging.Info)

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
function setup_global_logging(log_file::String; log_level::String = "info")
    # Convert string log level to LogLevel
    level = if log_level == "debug"
        Logging.Debug
    elseif log_level == "info"
        Logging.Info
    elseif log_level == "warn"
        Logging.Warn
    elseif log_level == "error"
        Logging.Error
    else
        Logging.Info  # Default fallback
    end
    
    # Store the current log level
    current_log_level[] = level
    
    # Close any existing global log file
    if !isnothing(global_log_handle[])
        close(global_log_handle[])
    end

    # Open new global log file
    global_log_handle[] = open(log_file, "w")

    # Record start time
    global_log_start_time[] = now()

    # Write timestamp as first line
    println(global_log_handle[], "Start time: ", Dates.format(global_log_start_time[], "dd-mm-yyyy HH:MM:SS\n"))

    # Set up the global logger, but only if we're not already file logging
    if !is_file_logging[]
        # Create a logger that writes to both stdout and file
        console_logger = ConsoleLogger(stdout, level)
        file_logger = FormatLogger(global_log_handle[]) do io, log
            if log.level >= level
                println(io, log.message)
            end
        end

        # Save current logger before replacing
        saved_logger[] = global_logger()

        # Combine both loggers and set as global
        logger = TeeLogger(console_logger, file_logger)
        global_logger(logger)
    end

    # Mark that global logging is active
    is_global_logging[] = true

    return global_log_handle[]
end

"""
    close_global_logging()

Close the global log file and add a closing timestamp with duration.
"""
function close_global_logging()
    if !isnothing(global_log_handle[])
        duration = now() - global_log_start_time[]
        println(global_log_handle[], "\nDuration: ", format_duration(duration))
        close(global_log_handle[])
        global_log_handle[] = nothing
        global_log_start_time[] = nothing

        # Only restore logger if we're not in file logging mode
        if !is_file_logging[] && is_global_logging[]
            global_logger(saved_logger[])
        end

        is_global_logging[] = false
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
function setup_logging(log_file::String; log_level::String = "info")
    # Convert string log level to LogLevel
    level = if log_level == "debug"
        Logging.Debug
    elseif log_level == "info"
        Logging.Info
    elseif log_level == "warn"
        Logging.Warn
    elseif log_level == "error"
        Logging.Error
    else
        Logging.Info  # Default fallback
    end
    
    # Store the current log level
    current_log_level[] = level
    
    # Close any existing log file
    if !isnothing(log_file_handle[])
        close(log_file_handle[])
    end

    # Open new log file (creates a new file each time)
    log_file_handle[] = open(log_file, "w")

    # Record start time
    log_start_time[] = now()

    # Write timestamp as first line
    println(log_file_handle[], "File Processing Started: ", Dates.format(log_start_time[], "dd-mm-yyyy HH:MM:SS\n"))

    # Only save logger if we're not already in file logging mode
    if !is_file_logging[]
        saved_logger[] = global_logger()
    end

    # Create a logger that writes to both stdout and file
    console_logger = ConsoleLogger(stdout, level)
    file_logger = FormatLogger(log_file_handle[]) do io, log
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

    # Mark that file logging is active
    is_file_logging[] = true

    return nothing
end

"""
    close_logging()

Close the individual log file and add a closing timestamp with duration.
Restores the previous logger (which may be the global logger).
"""
function close_logging()
    if !isnothing(log_file_handle[])
        duration = now() - log_start_time[]
        println(log_file_handle[], "\nFile Processing Duration: ", format_duration(duration))
        close(log_file_handle[])
        log_file_handle[] = nothing
        log_start_time[] = nothing

        # Mark that file logging is inactive
        is_file_logging[] = false

        # Restore the appropriate logger
        if is_global_logging[]
            # Recreate global logger
            console_logger = ConsoleLogger(stdout, current_log_level[])
            file_logger = FormatLogger(global_log_handle[]) do io, log
                if log.level >= current_log_level[]
                    println(io, log.message)
                end
            end
            global_logger(TeeLogger(console_logger, file_logger))
        else
            # No global logging, restore saved logger
            global_logger(saved_logger[])
        end
    end
end
