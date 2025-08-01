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
# Returns: "file"
```
"""
function basename_without_ext(path::String)
    return splitext(basename(path))[1]
end





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
    return [f(A[i+step], A[i]) for i = 1:(length(A)-step)]
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


function best_rect(n)
    dim1 = ceil(Int, sqrt(n))
    dim2 = ceil(Int, n ./ dim1)
    return [dim1, dim2]
end

"""
    orientation(p::Vector{Float64}, q::Vector{Float64}, r::Vector{Float64})

Helper function to find orientation of triplet (p, q, r).
Returns:
 0 --> p, q and r are collinear
 1 --> Clockwise
 2 --> Counterclockwise
"""
function orientation(p::Vector{Float64}, q::Vector{Float64}, r::Vector{Float64})
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val ≈ 0
        return 0
    end
    return val > 0 ? 1 : 2
end

"""
    create_convex_hull(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)

Create a convex hull around a set of 2D points with a specified border size.
Uses Graham's Scan algorithm for convex hull computation.

# Arguments
- `xpos`: Array of x-coordinates
- `ypos`: Array of y-coordinates
- `border_size`: Size of the border around points

# Returns
- A Vector of 2D points forming the convex hull
"""
function create_convex_hull(xpos::Vector{<:Real}, ypos::Vector{<:Real}, border_size::Real)
    # Generate points around each electrode with the border
    circle_points = 0:(2π/361):2π
    xs = (border_size .* sin.(circle_points) .+ transpose(xpos))[:]
    ys = (border_size .* cos.(circle_points) .+ transpose(ypos))[:]

    # Convert to array of points
    points = [[xs[i], ys[i]] for i in eachindex(xs)]
    n = length(points)

    # Find the bottommost point (and leftmost if tied)
    ymin = minimum(p -> p[2], points)
    p0 = points[findfirst(p -> p[2] == ymin, points)]

    # Sort points by polar angle with respect to p0
    sort!(points, by = p -> begin
        if p == p0
            return -Inf
        end
        return atan(p[2] - p0[2], p[1] - p0[1])
    end)

    # Initialize stack for Graham's scan
    stack = Vector{Vector{Float64}}()
    push!(stack, points[1])
    push!(stack, points[2])

    # Process remaining points
    for i = 3:n
        while length(stack) > 1 && orientation(stack[end-1], stack[end], points[i]) != 2
            pop!(stack)
        end
        push!(stack, points[i])
    end

    # Close the hull by connecting back to the first point
    push!(stack, stack[1])

    return stack
end


# Custom copy functions for main types to avoid deepcopy
"""
    Base.copy(dat::ContinuousData) -> ContinuousData

Create a copy of ContinuousData with copied DataFrames and analysis info.
The data and layout DataFrames are copied with `copycols=true` to ensure 
independence, while immutable fields are shared.
"""
function Base.copy(dat::ContinuousData)::ContinuousData
    return ContinuousData(copy(dat.data, copycols = true), copy(dat.layout), dat.sample_rate, copy(dat.analysis_info))
end

"""
    Base.copy(dat::EpochData) -> EpochData

Create a copy of EpochData with copied epoch DataFrames and analysis info.
Each epoch DataFrame in the data vector is copied independently.
"""
function Base.copy(dat::EpochData)::EpochData
    return EpochData(
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
