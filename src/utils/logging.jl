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
function _create_tee_logger(console_level::Logging.LogLevel, file_handle::IO, file_level::Logging.LogLevel; include_kwargs::Bool = false)
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

"""Helper to extract content between matching parentheses starting at given index."""
function _extract_paren_content(str::String, start_idx::Int)::String
    paren_count = 1
    for i in (start_idx + 1):length(str)
        if str[i] == '('
            paren_count += 1
        elseif str[i] == ')'
            paren_count -= 1
            paren_count == 0 && return str[(start_idx + 1):(i - 1)]
        end
    end
    return ""
end

"""Helper to parse closure representation and extract function name and args."""
function _parse_closure_repr(repr_str::String)::Union{Tuple{String,String},Nothing}
    # Pattern: eegfun.var"#samples##4#samples##5"{Tuple{Int64, Int64}}((0, 1))
    # Extract function name from var"#name##..."
    func_match = match(r"var\"#([^#]+)##", repr_str)
    isnothing(func_match) && return nothing
    
    # Extract arguments: find last }( (end of type params, start of args)
    brace_paren_range = findlast("}(", repr_str)
    if !isnothing(brace_paren_range)
        args_str = _extract_paren_content(repr_str, last(brace_paren_range))
    else
        # Fallback: try simple pattern without type parameters
        paren_match = match(r"\(([^)]*)\)$", repr_str)
        args_str = isnothing(paren_match) ? "" : paren_match.captures[1]
    end
    
    return (func_match.captures[1], args_str)
end

"""Helper to format a single kwarg value for logging."""
function _format_kwarg_value(k::Symbol, v)::String
    if v === nothing
        return "nothing"
    elseif isa(v, String)
        return "\"$v\""
    elseif isa(v, Function)
        repr_str = repr(v)
        parsed = _parse_closure_repr(repr_str)
        if !isnothing(parsed)
            func_name, args_str = parsed
            if isempty(args_str)
                return "eegfun.$func_name()"
            else
                return "eegfun.$func_name($args_str)"
            end
        end
        return repr_str
    else
        return string(v)
    end
end

"""
    _log_function_call(func_name::String, args::Vector, kwargs)

Log a function call in a generic way.

# Arguments
- `func_name::String`: Name of the function
- `args::Vector`: Positional arguments
- `kwargs`: Keyword arguments as pairs, named tuple, or dict
"""
function _log_function_call(
    func_name::String,
    args::Vector,
    kwargs::Union{Vector{Pair{Symbol,Any}},NamedTuple,Dict{Symbol,Any}},
)
    # Format positional arguments
    args_str = join(string.(args), ", ")

    # Convert to iterable pairs
    kw_pairs = kwargs isa NamedTuple ? pairs(kwargs) : kwargs

    # Format keyword arguments
    kwargs_str = join(["$k=$(_format_kwarg_value(k, v))" for (k, v) in kw_pairs], ", ")

    @info "Function call: $func_name($args_str; $kwargs_str)"
end

"""
    @log_call func_name (arg1, arg2, ...)

Macro to automatically log a function call by capturing local variables.
Any local variables not listed in the tuple are automatically captured as kwargs.

# Notes
- Positional arguments must be listed explicitly in the tuple
- Keyword arguments are automatically captured from local variables
- To capture kwargs, assign them to local variables before calling the macro
"""
macro log_call(func_name, args_tuple)
    # Validate inputs
    func_name isa String || error("@log_call: first argument must be a string (function name)")
    args_tuple isa Expr && args_tuple.head == :tuple || 
        error("@log_call: second argument must be a tuple of argument names")
    
    # Extract argument names and create set for filtering
    arg_names = Set(args_tuple.args)
    
    return quote
        local all_locals = Base.@locals()
        local kwargs_dict = Base.filter(p -> p.first ∉ $arg_names, all_locals)
        _log_function_call($(esc(func_name)), [$(esc.(args_tuple.args)...)], kwargs_dict)
    end
end

