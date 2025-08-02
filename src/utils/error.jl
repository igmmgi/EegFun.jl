"""
    @minimal_error(msg)

Displays an error message and stops execution without showing a full stacktrace.
"""
macro minimal_error(msg)
    # TODO: This is a hack to get a quick basic error message without a full stacktrace,
    # but only works if the error is not in a subfunction.
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
