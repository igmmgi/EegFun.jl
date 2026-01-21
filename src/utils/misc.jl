"""
    make_output_filename(output_dir::String, input_file::String, suffix::String)

Create an output filename from input file path with given suffix.

# Arguments
- `output_dir::String`: Output directory path
- `input_file::String`: Input file path
- `suffix::String`: Suffix to add (e.g., "_ica", "_continuous")

# Returns
- `String`: Full output filename path

# Example
```julia
filename = make_output_filename("/output", "data/file.bdf", "_ica")
# Returns: "/output/file_ica.jld2"
```
"""
function make_output_filename(output_dir::String, input_file::String, suffix::String)
    base_name = basename_without_ext(input_file)
    return joinpath(output_dir, "$(base_name)$(suffix).jld2")
end


"""
    basename_without_ext(path::String)

Extract the base filename without extension from a file path.

# Arguments
- `path::String`: File path to process

# Returns
- `String`: Base filename without extension

# Example
```julia
filename = basename_without_ext("data/file.bdf")
```
"""
basename_without_ext(path::String) = splitext(basename(path))[1]



"""
    consecutive(f::Function, A::AbstractVector; step::Int=1) -> Vector

Apply function f to consecutive pairs of elements in vector A.

# Arguments
- `f::Function`: Function to apply to pairs
- `A::AbstractVector`: Input vector
- `step::Int=1`: Step size between pairs

# Returns
- `Vector`: Results of applying f to consecutive pairs
"""
function consecutive(f::Function, A::AbstractVector; step::Int = 1)
    step < 1 && @minimal_error "Step must be positive"
    (length(A) < step + 1) && @minimal_error "Vector too short for given step size"
    return [f(A[i+step], A[i]) for i = 1:(length(A)-step)]
end

"""
    splitgroups(v::AbstractVector) -> Tuple{Vector{Int64},Vector{Int64}}

Split vector into groups based on consecutive numbers.

# Returns
- `Tuple{Vector{Int64},Vector{Int64}}`: Start and end indices of groups
"""
function splitgroups(v::AbstractVector{<:Integer})

    isempty(v) && return Int64[], Int64[]

    start, start_idx, end_idx = 1, Int64[], Int64[]
    for stop in [findall(diff(v) .> 1); lastindex(v)]
        push!(start_idx, v[start])
        push!(end_idx, v[stop])
        start = stop + 1
    end
    return start_idx, end_idx
end




"""
    get_channel_indices(dat::DataFrame, channel_labels::AbstractVector{<:AbstractString}) -> Vector{Int}

Get column indices for specified channel labels.

# Arguments
- `dat::DataFrame`: Data containing channels
- `channel_labels::AbstractVector{<:AbstractString}`: Channel labels to find

# Returns
- `Vector{Int}`: Column indices for requested channels

# Throws
- `ArgumentError`: If no matching channels found
"""
function get_channel_indices(dat::DataFrame, channel_labels::AbstractVector{<:AbstractString})::Vector{Int}
    isempty(channel_labels) && @minimal_error "channel_labels cannot be empty"

    channel_indices = findall(col -> col in channel_labels, names(dat))
    isempty(channel_indices) && @minimal_error "No matching channel_labels found in the data frame"

    return channel_indices
end





"""
    find_idx_range(time::AbstractVector, start_time::Real, end_time::Real) -> UnitRange{Int}
    find_idx_range(time::AbstractVector, limits::AbstractVector) -> UnitRange{Int}

Find index range corresponding to time interval.
Assumes time vector is sorted in ascending order.

# Returns
- `UnitRange{Int}`: Range of indices
"""
find_idx_range(time::AbstractVector, start_time::Real, end_time::Real) =
    searchsortedfirst(time, start_time):searchsortedlast(time, end_time)
find_idx_range(time::AbstractVector, limits::AbstractVector) = find_idx_range(time, limits[1], limits[end])


"""
    find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) -> Tuple{Int,Int}
    find_idx_start_end(time::AbstractVector, limits::AbstractVector) -> Tuple{Int,Int}

Find start and end indices corresponding to time interval.
Assumes time vector is sorted in ascending order.

# Returns
- `Tuple{Int,Int}`: Start and end indices
"""
find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) =
    searchsortedfirst(time, start_time), searchsortedlast(time, end_time)
find_idx_start_end(time::AbstractVector, limits::AbstractVector) =
    searchsortedfirst(time, limits[1]), searchsortedlast(time, limits[end])


"""
    find_times(time::AbstractVector, requested_times::AbstractVector) -> (indices::Vector{Int}, times::Vector{Float64})

Find nearest time points in a sorted time vector for each requested time.

For each requested time, finds the nearest matching time point in the sorted time vector
and returns both the indices and the actual time values.

# Arguments
- `time::AbstractVector`: Sorted time vector (assumed to be in ascending order)
- `requested_times::AbstractVector`: Vector of requested time points

# Returns
- `indices::Vector{Int}`: Indices of nearest matching time points
- `times::Vector{Float64}`: Actual time values at those indices

# Example
```julia
time_vec = [0.0, 0.01, 0.02, 0.03, 0.04]
requested = [0.005, 0.015, 0.025]
indices, times = find_times(time_vec, requested)
# indices = [1, 2, 3]
# times = [0.0, 0.01, 0.02]
```
"""
function find_times(time::AbstractVector, requested_times::AbstractVector)::Tuple{Vector{Int}, Vector{Float64}}
    time_min = minimum(time)
    time_max = maximum(time)
    
    indices = Int[]
    times_out = Float64[]
    
    for t_requested in requested_times
        # Only include if within data range
        if t_requested >= time_min && t_requested <= time_max
            # Use searchsortedfirst (same as find_idx_range/find_idx_start_end)
            idx = searchsortedfirst(time, t_requested)
            # Find nearest (check previous index if closer)
            if idx > 1 && abs(time[idx-1] - t_requested) < abs(time[min(idx, length(time))] - t_requested)
                idx = idx - 1
            end
            idx = min(idx, length(time))  # Ensure within bounds
            
            # Avoid duplicates
            if isempty(indices) || indices[end] != idx
                push!(indices, idx)
                push!(times_out, time[idx])
            end
        end
    end
    
    return indices, times_out
end



"""
    detrend(x::AbstractVector, y::AbstractVector) -> Vector{Float64}

Remove linear trend from data using least squares regression.

# Arguments
- `x::AbstractVector`: Independent variable (e.g., time points)
- `y::AbstractVector`: Dependent variable to detrend

# Returns
- `Vector{Float64}`: Detrended data with linear trend removed

# Example
```julia
x = 1:10
y = 2 .* x .+ randn(10)  # Linear trend with noise
y_detrended = detrend(x, y)
```
"""
function detrend(x::AbstractVector, y::AbstractVector)::Vector{Float64}
    length(x) == length(y) || @minimal_error "x and y must have the same length"
    length(x) < 2 && @minimal_error "Need at least 2 points for detrending"

    X = hcat(ones(length(x)), x)  # Design matrix (with intercept)
    β = X \ y  # Solve for coefficients (m, b)
    return y - (X * β)
end


"""
    extract_int(s::String) -> Union{Int, Nothing}

Extract the first integer found in a string.

# Arguments
- `s::String`: Input string

# Returns
- `Union{Int, Nothing}`: First integer found, or `nothing` if no digits

# Example
```julia
extract_int("channel_123_data")  # Returns: 123
extract_int("no_numbers_here")   # Returns: nothing
```
"""
function extract_int(s::String)::Union{Int,Nothing}
    digits_only = Base.filter(isdigit, s)
    return isempty(digits_only) ? nothing : parse(Int, digits_only)
end


"""
    natural_sort_key(s::String) -> String

Generate a sort key for natural/numeric sorting of strings.
Pads numeric parts with zeros so that "file_3" sorts before "file_12".

# Arguments
- `s::String`: Input string to generate sort key for

# Returns
- `String`: Transformed string with numeric parts zero-padded

# Example
```julia
# Use with sort
files = ["file_10.jld2", "file_2.jld2", "file_1.jld2"]
sort(files, by=natural_sort_key)  # Returns: ["file_1.jld2", "file_2.jld2", "file_10.jld2"]

# Use with DataFrame sort
sort(df, :file, by=natural_sort_key)
```
"""
natural_sort_key(s::String)::String = replace(s, r"\d+" => m -> lpad(String(m), 10, '0'))



"""
    @add_nonmutating function_name!

Create non-mutating versions of all methods of a function that ends with !.
Preserves method signatures and documentation.

Example:
@add_nonmutating filter_data!
"""
macro add_nonmutating(func)

    func_name = string(func)
    if !endswith(func_name, "!")
        error("Function name must end with !")
    end

    # base names
    non_mut_name = Symbol(func_name[1:(end-1)])
    mut_name = Symbol(func_name)

    # expressions for each method
    exprs = Expr(:block)

    # Collect all method signatures first
    methods_list = String[]
    method_exprs = []

    # all methods
    for method in methods(eval(func))

        sig = method.sig
        types = sig.parameters[2:end]

        # original parameter names
        params = Base.method_argnames(method)[2:end]
        if isempty(params) || any(==(nothing), params)
            params = [Symbol("arg", i) for i = 1:length(types)]
        end

        # signature string
        sig_str = "$non_mut_name(" * join(["$p::$t" for (p, t) in zip(params, types)], ", ") * "; kwargs...)"
        push!(methods_list, sig_str)

        # method definition without docstring
        method_expr = quote
            function $non_mut_name($([:($p::$t) for (p, t) in zip(params, types)]...); kwargs...)
                data_copy = copy($(params[1]))
                $mut_name(data_copy, $(params[2:end]...); kwargs...)
                return data_copy
            end
        end

        push!(method_exprs, method_expr)

    end

    # Create the main docstring 
    doc = """
        $(join(methods_list, "\n    "))

    Non-mutating version of `$mut_name`. Creates a copy of the input data
    and applies the operation to the copy.
    """

    push!(exprs.args, :(Base.@doc $doc $non_mut_name))

    # add all method definitions
    append!(exprs.args, method_exprs)

    return esc(exprs)

end


# Custom copy functions for main types to avoid deepcopy
"""
    Base.copy(dat::ContinuousData) -> ContinuousData

Create a copy of ContinuousData with copied DataFrames and analysis info.
The data and layout DataFrames are copied with `copycols=true` to ensure 
independence, while immutable fields are shared.
"""
function Base.copy(dat::ContinuousData)::ContinuousData
    return ContinuousData(
        dat.file,
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )
end

"""
    Base.copy(dat::EpochData) -> EpochData

Create a copy of EpochData with copied epoch DataFrames and analysis info.
Each epoch DataFrame in the data vector is copied independently.
"""
function Base.copy(dat::EpochData)::EpochData
    return EpochData(
        dat.file,
        dat.condition,
        dat.condition_name,
        [copy(epoch, copycols = true) for epoch in dat.data],
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
    )
end

"""
    Base.copy(dat::ErpData) -> ErpData

Create a copy of ErpData with copied data DataFrame and analysis info.
"""
function Base.copy(dat::ErpData)::ErpData
    return ErpData(
        dat.file,
        dat.condition,
        dat.condition_name,
        copy(dat.data, copycols = true),
        copy(dat.layout),
        dat.sample_rate,
        copy(dat.analysis_info),
        dat.n_epochs,
    )
end

"""
    Base.copy(info::AnalysisInfo) -> AnalysisInfo

Create a copy of AnalysisInfo. Since all fields are immutable (Symbol and Float64),
this creates a new instance with the same field values.
"""
function Base.copy(info::AnalysisInfo)::AnalysisInfo
    return AnalysisInfo(info.reference, info.hp_filter, info.lp_filter)
end

# Function to parse component input text into a list of component indices
"""
    parse_string_to_ints(text::String, total_components::Int)

Parses a comma-separated string potentially containing ranges (e.g., "1,3-5,8")
into a sorted, unique vector of valid component indices.

# Arguments
- `text::String`: The input string.
- `total_components::Int`: The maximum valid component index.

# Returns
- `Vector{Int}`: Sorted, unique vector of valid component indices found in the text.
                 Returns an empty vector if input is empty or invalid.
"""
function parse_string_to_ints(text::String)

    components = Int[]
    if isempty(text)
        return components
    end

    # Check for decimal points and error
    if occursin('.', text)
        throw(ArgumentError("Decimal points not allowed in int selection: '$text'"))
    end

    # Split by comma or semicolon and filter empty parts
    parts = Base.filter(!isempty, strip.(split(text, r"[,;]")))

    # Filter non-numeric parts (except :) and warn
    numeric_parts = []
    for part in parts
        if all(c -> isdigit(c) || c == ':', part)
            push!(numeric_parts, part)
        else
            @minimal_warning "Skipping non-numeric component: '$part'"
        end
    end

    for part in numeric_parts
        if occursin(':', part) # Handle ranges like "1:5"
            range_parts = strip.(split(part, ':'))
            if length(range_parts) == 2
                start_num = parse(Int, range_parts[1])
                end_num = parse(Int, range_parts[2])
                if start_num <= end_num
                    append!(components, start_num:end_num)
                end
            end
        else # Handle single numbers
            num = parse(Int, part)
            push!(components, num)
        end
    end

    # Remove duplicates and sort
    unique!(sort!(components))
    return components
end

function parse_string_to_ints(text::String, max_count::Int)
    all_components = parse_string_to_ints(text)
    return all_components[1:min(length(all_components), max_count)]
end


# Helper function to extract just the defaults
function _get_defaults(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::Dict{Symbol,Any}
    return Dict(key => value[1] for (key, value) in kwargs_dict)
end

"""
    _merge_plot_kwargs(defaults_dict, user_kwargs; validate=true)

Helper function to merge user keyword arguments with defaults and optionally validate parameter names.

# Arguments
- `defaults_dict`: Dictionary with (default_value, description) tuples
- `user_kwargs`: User-provided keyword arguments (NamedTuple or Dict)
- `validate`: Whether to validate unknown parameters (default: true)

# Returns
- `Dict{Symbol,Any}`: Merged and validated keyword arguments

# Example
```julia
# With validation (default):
plot_kwargs = _merge_plot_kwargs(PLOT_KWARGS, kwargs)

# Without validation (for multi-component functions):
plot_kwargs = _merge_plot_kwargs(PLOT_KWARGS, kwargs; validate=false)
```
"""
function _merge_plot_kwargs(
    defaults_dict::Dict{Symbol,Tuple{Any,String}},
    user_kwargs::NamedTuple;
    validate::Bool = true,
)::Dict{Symbol,Any}

    # Get default and user kwargs
    defaults = _get_defaults(defaults_dict)
    user_dict = Dict{Symbol,Any}(pairs(user_kwargs))

    # Check for unknown parameters (only if validation is enabled)
    if validate
        valid_keys = keys(defaults)
        unknown_keys = setdiff(keys(user_dict), valid_keys)
        if !isempty(unknown_keys)
            @minimal_error_throw "Unknown keyword arguments: $(join(unknown_keys, ", ")). Valid arguments: $(join(valid_keys, ", "))"
        end
    end

    # Merge defaults with user kwargs
    merged_kwargs = merge(defaults, user_dict)

    return merged_kwargs
end

# Convenience function for the common pattern
function _merge_plot_kwargs(
    defaults_dict::Dict{Symbol,Tuple{Any,String}},
    user_kwargs::Dict;
    validate::Bool = true,
)::Dict{Symbol,Any}
    return _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)
end

# Handle empty keyword arguments (Base.Pairs)
function _merge_plot_kwargs(
    defaults_dict::Dict{Symbol,Tuple{Any,String}},
    user_kwargs::Base.Pairs;
    validate::Bool = true,
)::Dict{Symbol,Any}
    return _merge_plot_kwargs(defaults_dict, NamedTuple(user_kwargs); validate = validate)
end



"""
    _orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64}) -> Int

Calculate the orientation of three points (p1, p2, p3).
Returns 0 if collinear, 1 if clockwise, 2 if counterclockwise.

This is a computational geometry utility function used in convex hull algorithms.
"""
function _orientation(p1::Vector{Float64}, p2::Vector{Float64}, p3::Vector{Float64})::Int
    val = (p2[2] - p1[2]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[2] - p2[2])
    if abs(val) < 1e-10
        return 0  # collinear
    elseif val > 0
        return 1  # clockwise
    else
        return 2  # counterclockwise
    end
end

# Helper function to generate documentation
function generate_kwargs_doc(kwargs_dict::Dict{Symbol,Tuple{Any,String}})::String
    doc_lines = ["# Keyword Arguments"]
    push!(doc_lines, "All keyword arguments below have sensible defaults defined in `DEFAULT_CHANNEL_SUMMARY_KWARGS`.")
    push!(doc_lines, "You can override any of these defaults by passing the corresponding keyword argument.")
    push!(doc_lines, "")
    for (param_name, (default_val, desc)) in kwargs_dict
        type_info = typeof(default_val)
        push!(doc_lines, "- `$(param_name)::$(type_info)=$(default_val)`: $(desc)")
    end
    return join(doc_lines, "\n")
end


"""
    combine_boolean_columns!(dat::ContinuousData, columns::Vector{Symbol}, operation::Symbol; output_column::Symbol = :combined_flags)

Combine multiple boolean columns using logical operations.

# Arguments
- `dat::ContinuousData`: The continuous EEG data object
- `columns::Vector{Symbol}`: Vector of column names to combine
- `operation::Symbol`: Logical operation to apply (:and, :or, :nand, :nor)
- `output_column::Symbol`: Name of the output column (default: :combined_flags)

# Modifies
- `dat`: Adds the combined boolean column to the data

# Examples
```julia
# Combine multiple artifact detection columns with AND operation
combine_boolean_columns!(dat, [:is_extreme_value_100, :is_eog_onset], :and)

# Combine with OR operation
combine_boolean_columns!(dat, [:is_extreme_value_100, :is_eog_onset], :or, output_column = :any_artifact)
```
"""
function combine_boolean_columns!(
    dat::ContinuousData,
    columns::Vector{Symbol},
    operation::Symbol;
    output_column::Symbol = :combined_flags,
)
    # Input validation
    @assert !isempty(columns) "Must specify at least one column to combine"
    @assert all(col -> hasproperty(dat.data, col), columns) "All specified columns must exist in the data"
    @assert operation in [:and, :or, :nand, :nor] "Invalid operation. Must be one of: :and, :or, :nand, :nor"

    # Get the boolean columns
    bool_columns = [dat.data[!, col] for col in columns]

    # Apply the logical operation
    result = if operation == :and
        all.(zip(bool_columns...))
    elseif operation == :or
        any.(zip(bool_columns...))
    elseif operation == :nand
        .!(all.(zip(bool_columns...)))
    elseif operation == :nor
        .!(any.(zip(bool_columns...)))
    end

    # Store the result
    dat.data[!, output_column] = result

    @info "combine_boolean_columns!: Combined $(length(columns)) columns using :$operation operation into column :$output_column"
end
