"""
    EegFunError <: Exception

User-facing exception raised by `@minimal_error` for invalid inputs, files, selections, or configurations.
"""
struct EegFunError <: Exception
    msg::String
end

Base.showerror(io::IO, e::EegFunError) = print(io, e.msg)

_log_error(msg::AbstractString) = @error "Error: " * msg _module = nothing _file = nothing _line = nothing

macro minimal_error(msg)
    quote
        local _msg = string($(esc(msg)))
        _file_logging_active() && _log_error(_msg)
        throw(EegFunError(_msg))
    end
end

macro minimal_warning(msg)
    quote
        @warn $(esc(msg)) _module = nothing _file = nothing _line = nothing
    end
end

macro minimal_stacktrace(msg, e, max_lines = 5)
    quote
        bt = catch_backtrace()
        error_msg = sprint(showerror, $(esc(e)), bt)
        st_lines = split(error_msg, '\n')
        limited_msg = join(st_lines[1:min($(esc(max_lines)), length(st_lines))], '\n')
        @error $(esc(msg)) limited_msg _module = nothing _file = nothing _line = nothing
    end
end
