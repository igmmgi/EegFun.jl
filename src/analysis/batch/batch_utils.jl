"""
    _log_function_call(func_name::String, args::Vector, kwargs::Vector{Pair{Symbol, Any}})
    _log_function_call(func_name::String, args::Vector, kwargs::NamedTuple)

Log a function call in a generic way.

# Arguments
- `func_name::String`: Name of the function
- `args::Vector`: Positional arguments
- `kwargs`: Keyword arguments as pairs or named tuple

# Example
```julia
_log_function_call("filter_data", ["epochs", 30.0], [:input_dir => "/path", :filter_type => "lp"])
_log_function_call("filter_data", ["epochs", 30.0], (; input_dir="/path", filter_type="lp"))
# Both log: Function call: filter_data("epochs", 30.0; input_dir="/path", filter_type="lp")
```
"""
# Helper to format a single kwarg value for logging
function _format_kwarg_value(k::Symbol, v)::String
    if v === nothing
        return "nothing"
    elseif isa(v, String)
        return "\"$v\""
    elseif isa(v, Function)
        # Special handling for common predicate functions
        if k in (:channel_selection, :component_selection, :epoch_selection, :sample_selection)
            return "<predicate>"
        else
            # Try to get a readable function name
            func_str = string(v)
            return occursin("#", func_str) ? "<function>" : func_str
        end
    else
        return string(v)
    end
end

function _log_function_call(func_name::String, args::Vector, kwargs::Union{Vector{Pair{Symbol, Any}}, NamedTuple, Dict{Symbol, Any}})
    # Format positional arguments
    args_str = join(string.(args), ", ")
    
    # Convert to iterable pairs
    kw_pairs = kwargs isa NamedTuple ? pairs(kwargs) : kwargs
    
    # Format keyword arguments
    kwargs_str = join(["$k=$(_format_kwarg_value(k, v))" for (k, v) in kw_pairs], ", ")
    
    @info "Function call: $func_name($args_str; $kwargs_str)"
end

"""
    @log_call
    @log_call func_name
    @log_call func_name n_args
    @log_call func_name args_tuple

Macro to automatically log a function call by capturing local variables.

# Arguments
- (none): Auto-detect function name from current scope, log all locals as args
- `func_name`: Name of the function as a string (optional, auto-detected if omitted)
- `n_args`: Number of positional arguments (integer), OR
- `args_tuple`: Tuple of positional argument variable names

# Examples
```julia
function my_function(x, y, z; opt1=1, opt2="test")
    # Simplest: auto-detect everything
    @log_call
    
    # Specify positional arg count
    @log_call "my_function" 3
    
    # Explicitly name args
    @log_call "my_function" (x, y, z)
end
```
"""
macro log_call(args...)
    func_name = nothing
    args_spec = nothing
    
    if length(args) == 0
        # No arguments: auto-detect from current scope
        func_name = String(__source__.file)  # Fallback
        return quote
            local all_locals = Base.@locals()
            _log_function_call("auto", collect(values(all_locals)), Dict{Symbol, Any}())
        end
    elseif length(args) == 1
        # Could be func_name only, or args_spec only
        if args[1] isa String
            # Just function name, log all locals as kwargs
            func_name = args[1]
            return quote
                local all_locals = Base.@locals()
                _log_function_call($(esc(func_name)), [], all_locals)
            end
        else
            error("Single argument must be function name string")
        end
    elseif length(args) == 2
        # func_name and args_spec
        func_name = args[1]
        args_spec = args[2]
        
        if args_spec isa Int
            # Integer: number of positional args
            n = args_spec
            return quote
                local all_locals = Base.@locals()
                local all_keys = collect(keys(all_locals))
                local args_keys = all_keys[1:$n]
                local kwargs_keys = all_keys[$n+1:end]
                
                local args_vals = [all_locals[k] for k in args_keys]
                local kwargs_dict = Dict(k => all_locals[k] for k in kwargs_keys)
                
                _log_function_call($(esc(func_name)), args_vals, kwargs_dict)
            end
        elseif args_spec isa Expr && args_spec.head == :tuple
            # Tuple: explicit argument names
            arg_names = Set(args_spec.args)
            return quote
                local all_locals = Base.@locals()
                local kwargs_dict = Base.filter(p -> p.first ∉ $arg_names, all_locals)
                _log_function_call(
                    $(esc(func_name)),
                    [$(esc.(args_spec.args)...)],
                    kwargs_dict
                )
            end
        else
            error("Second argument must be an integer or tuple of argument names")
        end
    else
        error("@log_call expects 0, 1, or 2 arguments")
    end
end
