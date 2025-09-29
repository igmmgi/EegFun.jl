"""
    _log_function_call(func_name::String, args::Vector, kwargs::Vector{Pair{Symbol, Any}})

Log a function call in a generic way.

# Arguments
- `func_name::String`: Name of the function
- `args::Vector`: Positional arguments
- `kwargs::Vector{Pair{Symbol, Any}}`: Keyword arguments as pairs

# Example
```julia
_log_function_call("filter_data", ["epochs", 30.0], [:input_dir => "/path", :filter_type => "lp"])
# Logs: Function call: filter_data("epochs", 30.0; input_dir="/path", filter_type="lp")
```
"""
function _log_function_call(func_name::String, args::Vector, kwargs::Vector{Pair{Symbol, Any}})
    # Format positional arguments - use string interpolation for cleaner output
    args_str = join([string(arg) for arg in args], ", ")
    
    # Format keyword arguments properly
    kwargs_parts = String[]
    for (k, v) in kwargs
        if v === nothing
            push!(kwargs_parts, "$k=nothing")
        elseif isa(v, String)
            push!(kwargs_parts, "$k=\"$v\"")
        else
            push!(kwargs_parts, "$k=$v")
        end
    end
    kwargs_str = join(kwargs_parts, ", ")
    
    @info "Function call: $func_name($args_str; $kwargs_str)"
end
