
function viewer(dat)
    ENV["TERM_PROGRAM"] == "vscode" ? vscodedisplay(dat) : display(dat)
end

function viewer(dat::EegData)
    viewer(data(dat))
end

function head(dat::EegData; n=nothing)
    isnothing(n) && (n=5)
    viewer(data(dat)[1:n, :])
end

function tail(dat::EegData; n=nothing)
    isnothing(n) && (n=5)
    viewer(data(dat)[end-n+1:end, :])
end

function to_data_frame(dat::EpochData)
    return vcat(dat.data...)
end

function to_data_frame(dat::Vector{EpochData})
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
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
    data_limits_x(dat::DataFrame) -> Tuple{Float64,Float64}

Get the time range of the data.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum time values
"""
data_limits_x(dat::DataFrame; col = :time) = extrema(dat[!, col])

"""
    data_limits_y(dat::DataFrame, col) -> Vector{Float64}

Get the value range for specified columns.

# Returns
- `Vector{Float64}`: [minimum, maximum] across specified columns
"""
data_limits_y(dat::DataFrame, col) = [minimum(Matrix(dat[!, col])), maximum(Matrix(dat[!, col]))]
