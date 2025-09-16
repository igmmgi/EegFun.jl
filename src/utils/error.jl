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
