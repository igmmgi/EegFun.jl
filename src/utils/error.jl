"""
    EegFunError <: Exception

Custom exception type for EegFun errors. Try and display a cleaner, user-friendly
error message without Julia's overly verbose stacktrace spew.

# Example
```julia
throw(EegFunError("No data found matching pattern 'xyz'"))
```
"""
struct EegFunError <: Exception
    msg::String
end

Base.showerror(io::IO, e::EegFunError) = print(io, e.msg)


"""
    @minimal_error(msg)

Throws an `EegFunError` with a clean user-facing error message.
Also logs the error via `@error` for diagnostics.
"""
macro minimal_error(msg)
    quote
        local _msg = string($(esc(msg)))
        @error "Error: " * _msg _module = nothing _file = nothing _line = nothing
        throw(EegFunError(_msg))
    end
end

"""
    @minimal_warning(msg)

Displays a warning message without showing module/file/line metadata.
"""
macro minimal_warning(msg)
    quote
        @warn $(esc(msg)) _module = nothing _file = nothing _line = nothing
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
        @error $(esc(msg)) limited_msg _module = nothing _file = nothing _line = nothing
    end
end
