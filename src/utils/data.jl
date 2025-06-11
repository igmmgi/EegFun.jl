function viewer(dat)
    if ENV["TERM_PROGRAM"] == "vscode"
        try
            vscodedisplay(dat)
        catch
            display(dat)
        end
    else
        display(dat)
    end
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
    isempty(dat.data) && return DataFrame()
    return vcat(dat.data...)
end

function to_data_frame(dat::Vector{EpochData})
    isempty(dat) && return DataFrame()
    isempty(dat[1].data) && return DataFrame()
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
colmeans(df::Matrix) = reduce(+, eachcol(df)) ./ size(df)[2]
colmeans(df::Matrix, cols) = reduce(+, eachcol(df[:, cols])) ./ length(cols)


"""
    data_limits_x(dat::DataFrame) -> Union{Tuple{Float64,Float64}, Nothing}

Get the time range of the data.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum time values
- `Nothing`: If the DataFrame is empty
"""
function data_limits_x(dat::DataFrame; col = :time)
    isempty(dat) && return nothing
    return extrema(dat[!, col])
end

"""
    data_limits_y(dat::DataFrame, col) -> Union{Vector{Float64}, Nothing}

Get the value range for specified columns.

# Returns
- `Vector{Float64}`: [minimum, maximum] across specified columns
- `Nothing`: If the DataFrame is empty
"""
function data_limits_y(dat::DataFrame, col)
    isempty(dat) && return nothing
    return [minimum(Matrix(dat[!, col])), maximum(Matrix(dat[!, col]))]
end
