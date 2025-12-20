"""
    @minimal_error(msg)

Displays an error message and stops execution without showing a full stacktrace.
"""
macro minimal_error(msg)
    quote
        @error "Error: " * $(esc(msg)) _module=nothing _file=nothing _line=nothing
        return nothing
    end
end

"""
    @minimal_error_throw(msg)

Displays an error message and throws an error without showing a full stacktrace.
"""
macro minimal_error_throw(msg)
    quote
        @error "Error: " * $(esc(msg)) _module=nothing _file=nothing _line=nothing
        error($(esc(msg)))
    end
end

"""
    @minimal_warning(msg)

Displays an error message and stops execution without showing a full stacktrace.
"""
macro minimal_warning(msg)
    quote
        @warn $(esc(msg)) _module=nothing _file=nothing _line=nothing
    end
end

"""
    @minimal_stacktrace(msg, e, max_lines=5)

Log an error with a limited stacktrace to avoid huge log files.

# Arguments
- `msg::String`: Error message
- `e`: The exception object
- `max_lines::Int`: Maximum number of stacktrace lines to include (default: 5)

# Example
```julia
catch e
    @minimal_stacktrace "Error processing file" e
    @minimal_stacktrace "Error processing file" e 10  # with custom max_lines
end
```
"""
macro minimal_stacktrace(msg, e, max_lines = 5)
    quote
        bt = catch_backtrace()
        error_msg = sprint(showerror, $(esc(e)), bt)
        st_lines = split(error_msg, '\n')
        limited_msg = join(st_lines[1:min($(esc(max_lines)), length(st_lines))], '\n')
        @error $(esc(msg)) limited_msg _module=nothing _file=nothing _line=nothing
    end
end
