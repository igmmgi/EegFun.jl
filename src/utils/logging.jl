
# Keep track of the open file handle and start time
const log_file_handle = Ref{Union{Nothing,IO}}(nothing)
const log_start_time = Ref{Union{Nothing,DateTime}}(nothing)

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
    setup_logging(log_file::String)

Set up logging to both stdout and a file. Use standard Julia logging macros (@info, @warn, etc.)
and messages will be written to both places. Creates a new log file each time it's called.

# Arguments
- `log_file::String`: Path to the log file

# Example
```julia
setup_logging("my_analysis.log")
@info "This will be logged to both stdout and the file"
```
"""
function setup_logging(log_file::String)
    # Close any existing log file
    if !isnothing(log_file_handle[])
        close(log_file_handle[])
    end
    
    # Open new log file (creates a new file each time)
    log_file_handle[] = open(log_file, "w")
    
    # Record start time
    log_start_time[] = now()
    
    # Write timestamp as first line
    println(log_file_handle[], "Date: ", Dates.format(log_start_time[], "dd-mm-yyyy HH:MM:SS\n"))
    
    # Create a logger that writes to both stdout and file
    console_logger = ConsoleLogger(stdout)
    file_logger = FormatLogger(log_file_handle[]) do io, args
        println(io, args.message)
    end
    
    # Combine both loggers
    logger = TeeLogger(console_logger, file_logger)
    global_logger(logger)
end

"""
    close_logging()

Close the log file and add a closing timestamp with duration. This should be called when you're done logging.
"""
function close_logging()
    if !isnothing(log_file_handle[])
        duration = now() - log_start_time[]
        println(log_file_handle[], "\nDuration: ", format_duration(duration))
        close(log_file_handle[])
        log_file_handle[] = nothing
        log_start_time[] = nothing
    end
    global_logger(ConsoleLogger(stdout))
end 