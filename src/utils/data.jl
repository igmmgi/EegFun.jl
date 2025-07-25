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
    subset_dataframe(dat::SingleDataFrameEeg; 
           channel_selection::Function = channels(), 
           sample_selection::Function = samples())

Create a subset of SingleDataFrameEeg (ContinuousData or ErpData) by applying channel and sample predicates.

# Arguments
- `dat`: SingleDataFrameEeg object to subset (ContinuousData or ErpData)
- `channel_selection`: Function that returns channel labels to include (default: all channels)
- `sample_selection`: Function that returns sample mask to include (default: all samples)

# Returns
- New SingleDataFrameEeg object with filtered channels and samples

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
function subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})
    dat_subset = df[selected_samples, :]
    dat_subset = select(dat_subset, selected_channels)
    return dat_subset
end

function subset(dat::SingleDataFrameEeg; 
               channel_selection::Function = channels(), 
               sample_selection::Function = samples(), 
               include_extra_channels::Bool = false)
   
    @info "subset: Subsetting $(typeof(dat)) ..."
    # Get selected channels and samples 
    selected_channels = get_selected_channels(dat, channel_selection, include_extra_channels = include_extra_channels)
    selected_samples = get_selected_samples(dat, sample_selection)
    
    # Filter data by samples and channels and match layout
    dat_subset = subset_dataframe(dat.data, selected_channels, selected_samples)
    layout_subset = subset_layout(dat.layout, channel_selection = channel_selection)
    
    # Create new SingleDataFrameEeg object
    return typeof(dat)(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info)

end

"""
    subset(dat::EpochData; 
           channel_selection::Function = channels(), 
           sample_selection::Function = samples(),
           epoch_selection::Function = epochs())

Create a subset of EpochData by applying channel, sample, and epoch predicates.

# Arguments
- `dat`: EpochData object to subset
- `channel_selection`: Function that returns channel labels to include (default: all channels)
- `sample_selection`: Function that returns sample mask to include (default: all samples)
- `epoch_selection`: Function that returns epoch mask to include (default: all epochs)

# Returns
- New EpochData object with filtered channels, samples, and epochs

# Examples
```julia
# Subset by channels only
filtered_dat = subset(dat, channel_selection = channels([:Fp1, :Fp2]))

# Subset by samples only
filtered_dat = subset(dat, sample_selection = x -> x.sample .< 1000)

# Subset by epochs only
filtered_dat = subset(dat, epoch_selection = epochs(1:10))

# Subset by all three
filtered_dat = subset(dat, 
    channel_selection = channels([:Fp1, :Fp2]), 
    sample_selection = x -> x.sample .< 1000,
    epoch_selection = epochs([1, 3, 5]))
```
"""
function subset(dat::EpochData; 
               channel_selection::Function = channels(), 
               sample_selection::Function = samples(),
               epoch_selection::Function = epochs(), 
               include_extra_channels::Bool = false)
    
    @info "subset: Subsetting $(typeof(dat)) ..."

    # Get selected channels, samples, and epochs
    selected_channels = get_selected_channels(dat, channel_selection, include_extra_channels = include_extra_channels)
    selected_samples = get_selected_samples(dat, sample_selection)
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    
    # Filter epochs first, then apply channel/sample filtering and match layout
    epochs_subset = subset_dataframe.(dat.data[selected_epochs], Ref(selected_channels), Ref(selected_samples))
    layout_subset = subset_layout(dat.layout, channel_selection = channel_selection)
    
    # Create new EpochData object
    return EpochData(epochs_subset, layout_subset, dat.sample_rate, dat.analysis_info)
end
