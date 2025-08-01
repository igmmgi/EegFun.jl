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
    @minimal_warning(msg)

Displays an error message and stops execution without showing a full stacktrace.
"""
macro minimal_warning(msg)
    quote
        @warn "Warning: " * $(esc(msg)) _module=nothing _file=nothing _line=nothing
    end
end
