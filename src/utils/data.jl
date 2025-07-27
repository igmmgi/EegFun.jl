# === DATAFRAME METADATA UTILITIES ===
"""
    _get_cols_by_group(df::DataFrame, group::Symbol) -> Vector{Symbol}

Get column names that belong to a specific metadata group.

# Arguments
- `df::DataFrame`: The DataFrame to search
- `group::Symbol`: The metadata group to search for

# Returns
- `Vector{Symbol}`: Column names that belong to the specified group

# Examples
```julia
cols = _get_cols_by_group(df, :channels)
```
"""
function _get_cols_by_group(df::DataFrame, group::Symbol)
    cols = propertynames(df)
    return [col for col in cols if haskey(metadata(df), string(col)) && metadata(df, string(col)) == ("group" => group)]
end

# === EEG DATA ACCESS FUNCTIONS ===
"""
    all_columns_data(dat::EegData) -> DataFrame

Get the complete DataFrame with all columns.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `DataFrame`: Complete DataFrame with all columns

# Examples
```julia
complete_df = all_columns_data(dat)
```
"""
all_columns_data(dat::SingleDataFrameEeg) = dat.data # single data frame
all_columns_data(dat::MultiDataFrameEeg) = to_data_frame(dat) # single data frame with all epochs

"""
    all_columns_labels(dat::EegData) -> Vector{Symbol}

Get all column names from the complete DataFrame.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: All column names

# Examples
```julia
all_cols = all_columns_labels(dat)
```
"""
all_columns_labels(dat::SingleDataFrameEeg) = propertynames(dat.data)
all_columns_labels(dat::MultiDataFrameEeg) = propertynames(dat.data[1])

# === EEG METADATA GROUP ACCESSORS ===
"""
    meta_data(eeg_data::EegData) -> DataFrame

Get metadata columns (time, sample, triggers) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing metadata columns

# Examples
```julia
meta_data = meta_data(dat)
```
"""
meta_data(dat::SingleDataFrameEeg) = dat.data[:, _get_cols_by_group(dat.data, :metadata)]
meta_data(dat::MultiDataFrameEeg, epoch::Int) = dat.data[epoch][:, _get_cols_by_group(dat.data, :metadata)]
meta_data(dat::MultiDataFrameEeg) = to_data_frame(dat)[:, _get_cols_by_group(dat.data, :metadata)] 

"""
    channel_data(eeg_data::EegData) -> DataFrame

Get EEG channel data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing EEG channel columns

# Examples
```julia
channel_data = channel_data(dat)
```
"""
channel_data(dat::SingleDataFrameEeg) = dat.data[:, _get_cols_by_group(dat.data, :channels)]
channel_data(dat::MultiDataFrameEeg, epoch::Int) = dat.data[epoch][:, _get_cols_by_group(dat.data, :channels)]
channel_data(dat::MultiDataFrameEeg) = to_data_frame(dat)[:, _get_cols_by_group(dat.data, :channels)] 

"""
    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing extra/derived columns

# Examples
```julia
extra_data = extra_data(dat)
```
"""
extra_data(dat::SingleDataFrameEeg) = dat.data[:, _get_cols_by_group(dat.data, :derived)]
extra_data(dat::MultiDataFrameEeg, epoch::Int) = dat.data[epoch][:, _get_cols_by_group(dat.data, :derived)]
extra_data(dat::MultiDataFrameEeg) = to_data_frame(dat)[:, _get_cols_by_group(dat.data, :derived)]

# === EEG METADATA LABEL ACCESSORS ===
"""
    meta_data_labels(eeg_data::EegData) -> Vector{Symbol}

Get metadata column names (time, sample, triggers) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of metadata column names

# Examples
```julia
meta_labels = meta_data_labels(dat)
```
"""
meta_data_labels(dat::SingleDataFrameEeg) = _get_cols_by_group(dat.data, :metadata)
meta_data_labels(dat::MultiDataFrameEeg, epoch::Int) = _get_cols_by_group(dat.data[epoch], :metadata)
meta_data_labels(dat::MultiDataFrameEeg) = _get_cols_by_group(to_data_frame(dat), :metadata)

"""
    channel_column_labels(eeg_data::EegData) -> Vector{Symbol}

Get EEG channel column names from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of channel column names

# Examples
```julia
channel_labels = channel_column_labels(dat)
```
"""
channel_column_labels(dat::SingleDataFrameEeg) = _get_cols_by_group(dat.data, :channels)
channel_column_labels(dat::MultiDataFrameEeg, epoch::Int) = _get_cols_by_group(dat.data[epoch], :channels)
channel_column_labels(dat::MultiDataFrameEeg) = _get_cols_by_group(to_data_frame(dat), :channels)

"""
    extra_column_labels(eeg_data::EegData) -> Vector{Symbol}

Get extra/derived column names (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of extra column names

# Examples
```julia
extra_labels = extra_column_labels(dat)
```
"""
extra_column_labels(dat::SingleDataFrameEeg) = _get_cols_by_group(dat.data, :derived)
extra_column_labels(dat::MultiDataFrameEeg, epoch::Int) = _get_cols_by_group(dat.data[epoch], :derived)
extra_column_labels(dat::MultiDataFrameEeg) = _get_cols_by_group(to_data_frame(dat), :derived)

# === EEG CONVENIENCE FUNCTIONS ===
# Basic information functions
"""
    sample_rate(dat::EegData) -> Int

Get the sample rate of the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Sample rate in Hz

# Examples
```julia
fs = sample_rate(dat)
```
"""
sample_rate(dat::EegData) = dat.sample_rate
sample_rate(dat::DataFrame) = Int(1 / mean(diff(dat.time)))

"""
    reference(dat::EegData) -> String

Get the reference information from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `String`: Reference information

# Examples
```julia
ref = reference(dat)
```
"""
reference(dat::EegData) = dat.analysis_info.reference
reference(dat::AnalysisInfo) = dat.reference

"""
    filter_info(dat::AnalysisInfo) -> Vector

Get filter information from the analysis info.

# Arguments
- `dat::AnalysisInfo`: The analysis info object

# Returns
- `Vector`: Filter information [hp_filter, lp_filter]

# Examples
```julia
filters = filter_info(dat.analysis_info)
```
"""
filter_info(dat::AnalysisInfo) = [dat.hp_filter, dat.lp_filter]

# Unique convenience functions (not available through metadata)
"""
    n_epochs(dat::EegData) -> Int

Get the number of epochs in the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Number of epochs

# Examples
```julia
n = n_epochs(dat)
```
"""
n_epochs(dat::SingleDataFrameEeg) = 1
n_epochs(dat::MultiDataFrameEeg) = length(dat.data)

"""
    duration(dat::EegData) -> Float64

Get the duration of the EEG data in seconds.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Float64`: Duration in seconds

# Examples
```julia
dur = duration(dat)
```
"""
duration(dat::EegData) = last(meta_data(dat).time) - first(meta_data(dat).time)

"""
    n_average(dat::ErpData) -> Int

Get the number of averaged epochs in ERP data.

# Arguments
- `dat::ErpData`: The ERP data object

# Returns
- `Int`: Number of averaged epochs

# Examples
```julia
n = n_average(dat)
```
"""
n_average(dat::ErpData) = dat.n_epochs

# EEG channel information
"""
    has_channels(dat::EegData, chans::Vector{Symbol}) -> Bool

Check if the EEG data contains all specified channels.

# Arguments
- `dat::EegData`: The EEG data object
- `chans::Vector{Symbol}`: Vector of channel symbols to check

# Returns
- `Bool`: True if all channels are present

# Examples
```julia
has_all = has_channels(dat, [:Fp1, :Fp2])
```
"""
has_channels(dat::EegData, chans::Vector{Symbol}) = all(in(channel_column_labels(dat)), chans)

"""
    common_channels(dat1::EegData, dat2::EegData) -> Vector{Symbol}

Find common channels between two EEG data objects.

# Arguments
- `dat1::EegData`: First EEG data object
- `dat2::EegData`: Second EEG data object

# Returns
- `Vector{Symbol}`: Vector of common channel symbols

# Examples
```julia
common = common_channels(dat1, dat2)
```
"""
common_channels(dat1::EegData, dat2::EegData) = intersect(channel_column_labels(dat1), channel_column_labels(dat2))

# === DATA VIEWING FUNCTIONS ===
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
