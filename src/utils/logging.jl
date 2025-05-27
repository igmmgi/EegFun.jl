using Logging
using Dates
using LoggingExtras

# Keep track of the open file handle
const log_file_handle = Ref{Union{Nothing,IO}}(nothing)

"""
    setup_logging(log_file::String)

Set up logging to both stdout and a file. Use standard Julia logging macros (@info, @warn, etc.)
and messages will be written to both places.

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
    
    # Open new log file
    log_file_handle[] = open(log_file, "a")
    
    # Create a logger that writes to both stdout and file
    console_logger = ConsoleLogger(stdout)
    file_logger = FormatLogger(log_file_handle[]) do io, args
        println(io, "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] ", args.message)
    end
    
    # Combine both loggers
    logger = TeeLogger(console_logger, file_logger)
    global_logger(logger)
end

"""
    reset_logging()

Reset logging to default stdout only and close the log file.
"""
function reset_logging()
    if !isnothing(log_file_handle[])
        close(log_file_handle[])
        log_file_handle[] = nothing
    end
    global_logger(ConsoleLogger(stdout))
end

# Function to write to both file and console
function log_message(level::String, msg::String)
    if !isnothing(log_file_handle[])
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        println(log_file_handle[], "[$timestamp] $level: $msg")
        flush(log_file_handle[])
    end
    println(stderr, "$level: $msg")
end 