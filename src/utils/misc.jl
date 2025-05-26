
"""
    channel_number_to_channel_label(channel_labels::Vector{Symbol}, channel_numbers::Union{Int,Vector{Int},UnitRange}) -> Vector{Symbol}

Convert channel numbers to their corresponding labels.

# Arguments
- `channel_labels::Vector{Symbol}`: List of all channel labels
- `channel_numbers`: Channel number(s) to convert (can be Int, Vector{Int}, or UnitRange)

# Returns
- `Vector{Symbol}`: Channel labels corresponding to the input numbers
"""
function channel_number_to_channel_label(channel_labels, channel_numbers::Union{Int,Vector{Int},UnitRange})
    return channel_labels[channel_numbers]
end




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
    if step < 1
        throw(ArgumentError("Step must be positive"))
    end
    if length(A) < step + 1
        throw(ArgumentError("Vector too short for given step size"))
    end
    return [f(A[i+step], A[i]) for i = 1:length(A)-step]
end

"""
    splitgroups(v::AbstractVector) -> Tuple{Vector{Int64},Vector{Int64}}

Split vector into groups based on consecutive numbers.

# Returns
- `Tuple{Vector{Int64},Vector{Int64}}`: Start and end indices of groups
"""
function splitgroups(v::AbstractVector{<:Integer})
    if isempty(v)
        return Int64[], Int64[]
    end

    start = 1
    start_idx = Int64[]
    end_idx = Int64[]

    for stop in [findall(diff(v) .> 1); lastindex(v)]
        push!(start_idx, v[start])
        push!(end_idx, v[stop])
        start = stop + 1
    end
    return start_idx, end_idx
end

# data limits




"""
    validate_baseline_interval(time::AbstractVector, baseline_interval::Union{IntervalIdx,IntervalTime}) -> IntervalIdx

Validate and convert baseline interval to index format.

# Arguments
- `time::AbstractVector`: Time points vector
- `baseline_interval::Union{IntervalIdx,IntervalTime}`: Interval to validate

# Returns
- `IntervalIdx`: Validated interval in index format

# Throws
- `ArgumentError`: If interval is invalid
"""
function validate_baseline_interval(
    time::AbstractVector,
    baseline_interval::Union{IntervalIdx,IntervalTime},
)::IntervalIdx
    if baseline_interval isa IntervalTime
        baseline_interval =
            IntervalIdx(find_idx_start_end(time, baseline_interval.interval_start, baseline_interval.interval_end)...)
    end

    if !(1 <= baseline_interval.interval_start <= length(time)) ||
       !(1 <= baseline_interval.interval_end <= length(time)) ||
       !(baseline_interval.interval_start <= baseline_interval.interval_end)
        throw(ArgumentError("Invalid baseline_interval: $baseline_interval"))
    end

    return baseline_interval
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
    if isempty(channel_labels)
        throw(ArgumentError("channel_labels cannot be empty"))
    end

    channel_indices = findall(col -> col in channel_labels, names(dat))
    if isempty(channel_indices)
        throw(ArgumentError("No matching channel_labels found in the data frame"))
    end

    return channel_indices
end


"""
    search_sequence(array::AbstractVector, sequence::Int) -> Vector{Int}

Find starting indices of a sequence in an array.

# Returns
- `Vector{Int}`: Indices where sequence starts
"""
search_sequence(array::AbstractVector, sequence::Int) =
    intersect(findall(array .== sequence), findall(diff(vcat(0, array)) .>= 1))


"""
    find_idx_range(time::AbstractVector, start_time::Real, end_time::Real) -> UnitRange{Int}
    find_idx_range(time::AbstractVector, limits::AbstractVector) -> UnitRange{Int}

Find index range corresponding to time interval.

# Returns
- `UnitRange{Int}`: Range of indices
"""
find_idx_range(time::AbstractVector, start_time::Real, end_time::Real) =
    findmin(abs.(time .- start_time))[2]:findmin(abs.(time .- end_time))[2]
find_idx_range(time::AbstractVector, limits::AbstractVector) = find_idx_range(time, limits[1], limits[end])


"""
    find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) -> Tuple{Int,Int}
    find_idx_start_end(time::AbstractVector, limits::AbstractVector) -> Tuple{Int,Int}

Find start and end indices corresponding to time interval.

# Returns
- `Tuple{Int,Int}`: Start and end indices
"""
find_idx_start_end(time::AbstractVector, start_time::Real, end_time::Real) =
    findmin(abs.(time .- start_time))[2], findmin(abs.(time .- end_time))[2]
find_idx_start_end(time::AbstractVector, limits::AbstractVector) =
    findmin(abs.(time .- limits[1]))[2], findmin(abs.(time .- limits[end]))[2]



function detrend(x, y)
    X = hcat(ones(length(x)), x)  # Design matrix (with intercept)
    β = X \ y  # Solve for coefficients (m, b)
    return y - (X * β)
end


function extract_int(s::String)
    digits_only = filter(isdigit, s)
    return isempty(digits_only) ? nothing : parse(Int, digits_only)
end



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
    non_mut_name = Symbol(func_name[1:end-1])
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
            params = [Symbol("arg", i) for i in 1:length(types)]
        end
        
        # signature string
        sig_str = "$non_mut_name(" * join(["$p::$t" for (p,t) in zip(params, types)], ", ") * "; kwargs...)"
        push!(methods_list, sig_str)
        
        # method definition without docstring
        method_expr = quote
            function $non_mut_name($([:($p::$t) for (p,t) in zip(params, types)]...); kwargs...)
                data_copy = deepcopy($(params[1]))
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

    See also: [`$func_name`](@ref)
    """
    
    push!(exprs.args, :(Base.@doc $doc $non_mut_name))
    
    # add all method definitions
    append!(exprs.args, method_exprs)
    
    return esc(exprs)

end


function best_rect(n)
    dim1 = ceil(Int, sqrt(n))
    dim2 = ceil(Int, n ./ dim1)
    return [dim1, dim2]
end


