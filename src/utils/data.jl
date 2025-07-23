function viewer(dat)
    if ENV["TERM_PROGRAM"] == "vscode"
        try
            Main.vscodedisplay(dat)
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

"""
    subset(dat::ContinuousData; 
           channel_selection::Function = channels(), 
           sample_selection::Function = samples())

Create a subset of ContinuousData by applying channel and sample predicates.

# Arguments
- `dat`: ContinuousData object to subset
- `channel_selection`: Function that returns channel labels to include (default: all channels)
- `sample_selection`: Function that returns sample mask to include (default: all samples)

# Returns
- New ContinuousData object with filtered channels and samples

# Examples
```julia
# Subset by channels only
filtered_dat = subset(dat, channel_selection = channels([:Fp1, :Fp2]))

# Subset by samples only
filtered_dat = subset(dat, sample_selection = x -> x.sample .< 1000)

# Subset by both
filtered_dat = subset(dat, 
    channel_selection = channels([:Fp1, :Fp2]), 
    sample_selection = x -> x.sample .< 1000)
```
"""
function subset(dat::ContinuousData; 
               channel_selection::Function = channels(), 
               sample_selection::Function = samples())
    
    # Get selected channels and samples directly
    selected_channels = get_selected_channels(dat, channel_selection)
    selected_samples = get_selected_samples(dat, sample_selection)
    
    # Filter data by samples and channels
    dat_subset = dat.data[selected_samples, :]
    dat_subset = select(dat_subset, vcat([:time, :sample, :triggers], selected_channels))
    
    # Filter layout to match selected channels
    layout_subset = filter(:label => in(selected_channels), dat.layout)
    
    # Create new ContinuousData object
    return ContinuousData(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info)
end
