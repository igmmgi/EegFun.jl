"""
Generic utility functions for batch processing operations.
These are reusable across different batch scripts.
"""

#=============================================================================
    DATA STRUCTURES
=============================================================================#

"""Immutable struct to hold batch processing results."""
struct BatchResult
    success::Bool
    filename::String
    message::String
end

"""Configuration for batch operations."""
struct BatchConfig
    file_pattern::String
    input_dir::String
    output_dir::String
    participants::Union{Int,Vector{Int},Nothing}
    conditions::Union{Int,Vector{Int},Nothing}
end

#=============================================================================
    PURE FUNCTIONS (no side effects)
=============================================================================#

"""
    find_batch_files(pattern::String, dir::String; participants=nothing)

Find JLD2 files matching pattern and optional participant filter.

Returns vector of filenames (not full paths).
"""
function _find_batch_files(pattern::String, dir::String; participants = nothing)
    all_files = readdir(dir)

    # Filter by pattern and extension
    files = Base.filter(all_files) do f
        endswith(f, ".jld2") && contains(f, pattern)
    end

    # Filter by participant if specified
    if participants !== nothing
        files = _filter_files(files; include = participants)
    end

    return files
end

"""
    load_eeg_data(filepath::String)

Load EEG data from JLD2 file, returning `(data_var, var_name)` or `nothing`.

Tries common variable names: "erps", "epochs".
"""
function _load_eeg_data(filepath::String)
    file_data = load(filepath)
    for var_name in ["erps", "epochs"]
        if haskey(file_data, var_name)
            return (file_data[var_name], var_name)
        end
    end
    return nothing
end

"""
    _condition_select(data, conditions)

Filter data by condition indices.

Returns filtered data or original if `conditions` is `nothing`.
"""
function _condition_select(data, conditions)
    isnothing(conditions) && return data
    isempty(data) && return data
    condition_nums = conditions isa Int ? [conditions] : conditions

    # Check if all requested conditions exist
    max_condition = length(data)
    invalid_conditions = condition_nums[condition_nums .> max_condition]
    if !isempty(invalid_conditions)
        throw(ArgumentError("Requested conditions $invalid_conditions exceed available conditions (1:$max_condition)"))
    end

    return data[condition_nums]
end

"""
    validate_input_dir(dir::String)

Validate input directory exists, returning error message or nothing.
"""
function _validate_input_dir(dir::String)
    !isdir(dir) && return "Input directory does not exist: $dir"
    return nothing
end

"""
    validate_channel_groups(groups::Vector{Vector{Symbol}})

Validate channel groups, returning error message or nothing.
Issues warnings for groups with < 2 channels.
"""
function _validate_channel_groups(groups::Vector{Vector{Symbol}})
    isempty(groups) && return "Channel groups cannot be empty"

    for (i, group) in enumerate(groups)
        isempty(group) && return "Channel group $i is empty"

        if length(group) < 2
            @minimal_warning "Channel group $i has only $(length(group)) channel(s): $group. Consider using more channels for meaningful averaging."
        end
    end

    return nothing
end

"""
    validate_condition_groups(groups::Vector{Vector{Int}})

Validate condition groups, returning error message or nothing.
Modifies groups in-place to remove duplicates and warns about overlaps.
"""
function _validate_condition_groups(groups::Vector{Vector{Int}})
    isempty(groups) && return "Condition groups cannot be empty"

    all_conditions = Int[]

    for (group_idx, group) in enumerate(groups)
        # Check for duplicates within the group
        if length(group) != length(unique(group))
            duplicates = group[findall(x -> count(==(x), group) > 1, group)]
            @minimal_warning "Group $group_idx contains duplicate conditions: $duplicates. Only unique conditions will be used."
        end

        # Remove duplicates and update the group
        unique_group = unique(group)
        if length(unique_group) != length(group)
            @info "  Group $group_idx: $(group) → $(unique_group) (removed duplicates)"
        end
        groups[group_idx] = unique_group

        # Check for overlaps between groups
        overlap = intersect(all_conditions, unique_group)
        if !isempty(overlap)
            @minimal_warning "Condition(s) $overlap appear in multiple groups. This may not be intended."
        end

        append!(all_conditions, unique_group)
    end

    return nothing
end

"""
    _validate_condition_pairs(pairs)

Validate condition pairs, returning error message or nothing.
Warns about pairs with identical conditions.
Accepts both `Vector{Tuple{Int, Int}}` and `Vector{Vector{Int}}`.
"""
function _validate_condition_pairs(pairs::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}})
    isempty(pairs) && return "Condition pairs cannot be empty"

    for (i, pair) in enumerate(pairs)
        if pair[1] == pair[2]
            @minimal_warning "Condition pair $i: ($(pair[1]), $(pair[2])) has identical conditions. Difference will be zero."
        end
    end

    return nothing
end

#=============================================================================
    ORCHESTRATION FUNCTIONS (with side effects: I/O, logging)
=============================================================================#

"""
    run_batch_operation(process_fn::Function, files::Vector{String}, 
                        input_dir::String, output_dir::String;
                        operation_name::String = "Processing")

Execute batch processing with logging and error handling.

Takes a pure processing function `process_fn(input_path, output_path)` and maps 
it over files, handling errors and logging progress.

Returns vector of `BatchResult`.
"""
function _run_batch_operation(
    process_fn::Function,
    files::Vector{String},
    input_dir::String,
    output_dir::String;
    operation_name::String = "Processing",
)
    n_files = length(files)
    results = Vector{BatchResult}(undef, n_files)

    for (i, file) in enumerate(files)
        @info "$operation_name: $file ($i/$n_files)"

        input_path = joinpath(input_dir, file)
        output_path = joinpath(output_dir, file)

        result = try
            process_fn(input_path, output_path)
        catch e
            @error "Error processing $file" exception=(e, catch_backtrace())
            BatchResult(false, file, "Exception: $(sprint(showerror, e))")
        end

        results[i] = result

        # Log individual result
        if result.success
            @info "  ✓ $(result.message)"
        else
            @minimal_warning "  ✗ $(result.message)"
        end
    end

    return results
end

"""
    log_batch_summary(results::Vector{BatchResult}, output_dir::String)

Log summary statistics from batch results.

Returns named tuple `(success=n, errors=n)`.
"""
function _log_batch_summary(results::Vector{BatchResult}, output_dir::String)
    n_success = count(r -> r.success, results)
    n_error = length(results) - n_success

    @info "Batch operation complete! Processed $n_success files successfully, $n_error errors"
    @info "Output saved to: $output_dir"

    return (success = n_success, errors = n_error)
end

"""
Cleanup global logging and move log file to output directory.
If output_dir is nothing, just closes logging without moving the file.
"""
function _cleanup_logging(log_file::String, output_dir::Union{String,Nothing})
    close_global_logging()

    if !isnothing(output_dir)
        log_dest = joinpath(output_dir, log_file)
        log_file != log_dest && mv(log_file, log_dest, force = true)
    end
end

#=============================================================================
    LOGGING FUNCTIONS
=============================================================================#

"""Helper to format a single kwarg value for logging."""
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

# Examples
```julia
function my_function(x, y, z; opt1=1, opt2="test")
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
            _log_function_call("auto", collect(values(all_locals)), Dict{Symbol,Any}())
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
                local args_keys = all_keys[1:($n)]
                local kwargs_keys = all_keys[($n+1):end]

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
                _log_function_call($(esc(func_name)), [$(esc.(args_spec.args)...)], kwargs_dict)
            end
        else
            error("Second argument must be an integer or tuple of argument names")
        end
    else
        error("@log_call expects 0, 1, or 2 arguments")
    end
end
