"""
    check_files_exist(conditions::Union{Vector{Int}, Int} filetype::String) -> Bool

Check if files exist for all given conditions with specified filetype.

# Arguments
- `conditions::Vector{String}`: List of condition names
- `filetype::String`: Type of file to check

# Returns
- `Bool`: true if all files exist, false otherwise
"""
function check_files_exist(conditions::Union{Vector{Int},Int}, filetype::String)
    for condition in conditions
        fname = "$(condition)_$(filetype).jld2"
        if !isfile(fname)
            @warn "File not found: $(fname)"
            all_files_exist = false
        end
    end
    return all_files_exist
end

"""
    check_files_exist(subjects::Union{Vector{Int}, Int}, conditions::Union{Vector{Int}, Int},, filetype::String) -> Bool

Check if files exist for all combinations of subjects and conditions.

# Arguments
- `subjects::Vector{String}`: List of subject identifiers
- `conditions::Vector{String}`: List of condition names
- `filetype::String`: Type of file to check

# Returns
- `Bool`: true if all files exist, false otherwise
"""
function check_files_exist(subjects::Union{Vector{Int},Int}, conditions::Union{Vector{Int},Int}, filetype::String)
    all_files_exist = true
    for subject in subjects
        for condition in conditions
            fname = "$(subject)_$(condition)_$(filetype).jld2"
            if !isfile(fname)
                @warn "File not found: $(fname)"
                all_files_exist = false
            end
        end
    end
    return all_files_exist
end

"""
    channel_number_to_channel_label(channel_labels::Vector{String}, channel_numbers::Union{Int64,Vector{Int64}}) -> Vector{String}

Convert channel numbers to their corresponding labels.

# Arguments
- `channel_labels::Vector{String}`: List of all channel labels
- `channel_numbers::Union{Int64,Vector{Int64}}`: Channel number(s) to convert

# Returns
- `Vector{String}`: Channel labels corresponding to the input numbers
"""
channel_number_to_channel_label(channel_labels, channel_numbers::Int64) = [channel_labels[channel_numbers]]
channel_number_to_channel_label(channel_labels, channel_numbers::Vector{Int64}) = channel_labels[channel_numbers]
channel_number_to_channel_label(channel_labels, channel_numbers::UnitRange) = channel_labels[channel_numbers]


function print_vector_(v::Vector; max_length::Int = 10, n_ends::Int = 5)
    if length(v) > max_length
        v = vcat(first(v, n_ends), "...", last(v, n_ends))
    end
    return join(v, ", ")
end
   
function print_vector(v::UnitRange; max_length::Int = 10, n_ends::Int = 5)
    print_vector_(collect(v), max_length=max_length, n_ends=n_ends)
end


"""
    datarange(x::AbstractVector) -> Float64

Calculate the range of data (maximum - minimum).

# Returns
- `Float64`: Difference between maximum and minimum values

# Example
```julia
datarange([1.0, 2.0, 3.0]) # returns 2.0
```
"""
datarange(x::AbstractVector) = -(-(extrema(x)...))


"""
    colmeans(df::DataFrame, cols) -> Vector{Float64}
    colmeans(df::Matrix) -> Vector{Float64}
    colmeans(df::Matrix, cols) -> Vector{Float64}

Calculate the mean of specified columns in a DataFrame.

# Arguments
- `df::DataFrame`: The DataFrame containing the data.
- `cols`: The columns for which to calculate the mean. This can be a vector of column names or indices.

# Returns
- `Vector{Float64}`: A vector containing the mean of each specified column.
"""
colmeans(df::DataFrame, cols) = reduce(+, eachcol(df[!, cols])) ./ length(cols)
colmeans(df::Matrix) = reduce(+, eachrow(df)) ./ size(df)[1]
colmeans(df::Matrix, cols) = reduce(+, eachrow(df[:, cols])) ./ size(df)[1]


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
    data_limits_x(dat::DataFrame) -> Tuple{Float64,Float64}

Get the time range of the data.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum time values
"""
data_limits_x(dat::DataFrame) = extrema(dat.time)

"""
    data_limits_y(dat::DataFrame, col) -> Vector{Float64}

Get the value range for specified columns.

# Returns
- `Vector{Float64}`: [minimum, maximum] across specified columns
"""
data_limits_y(dat::DataFrame, col) = [minimum(Matrix(dat[!, col])), maximum(Matrix(dat[!, col]))]



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


function to_data_frame(dat::EpochData)
    return vcat(dat.data...)
end

function to_data_frame(dat::Vector{EpochData})
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end
