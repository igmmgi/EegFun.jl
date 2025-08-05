# === COLUMN IDENTIFICATION SYSTEM ===
"""
    get_cols_by_group(dat::EegData, group::Symbol) -> Vector{Symbol}

Get columns by group type for EegData objects using layout-based identification.

# Arguments
- `dat::EegData`: The EEG data object
- `group::Symbol`: The group type (:metadata, :channels, :extra)

# Returns
- `Vector{Symbol}`: Column names of the specified group

# Group Types
- `:channels`: EEG channel columns (intersection of layout labels and DataFrame columns)
- `:metadata`: System columns (all columns before first layout label)
- `:extra`: Derived columns (all columns after last layout label)

"""
function get_cols_by_group(dat::EegData, group::Symbol)::Vector{Symbol}
    labels = all_labels(dat)
    channel_labels = dat.layout.data.label

    if group == :channels
        return intersect(channel_labels, labels)
    elseif group == :metadata
        if isempty(channel_labels)
            return Symbol[]
        end
        first_channel_idx = findfirst(col -> col == channel_labels[1], labels)
        isnothing(first_channel_idx) && @minimal_error "First channel label not found in data"
        return labels[1:(first_channel_idx-1)]
    elseif group == :extra
        if isempty(channel_labels)
            return Symbol[]
        end
        last_channel_idx = findlast(col -> col == channel_labels[end], labels)
        isnothing(last_channel_idx) && @minimal_error "Last channel label not found in data"
        return labels[(last_channel_idx+1):end]
    else
        @minimal_error "Unknown group type: $group"
    end
end


# === EEG DATA ACCESS FUNCTIONS ===
"""
    all_data(dat::EegData) -> DataFrame

Get the complete DataFrame with all columns.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `DataFrame`: Complete DataFrame with all columns
"""
all_data(dat::SingleDataFrameEeg)::DataFrame = dat.data # single data frame
all_data(dat::MultiDataFrameEeg)::DataFrame = to_data_frame(dat) # single data frame with all epochs

"""
    all_labels(dat::EegData) -> Vector{Symbol}

Get all column names from the complete DataFrame.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: All column names
"""
all_labels(dat::SingleDataFrameEeg)::Vector{Symbol} = propertynames(dat.data)
all_labels(dat::MultiDataFrameEeg)::Vector{Symbol} = propertynames(dat.data[1])
all_labels(dat::MultiDataFrameEeg, epoch::Int)::Vector{Symbol} = propertynames(dat.data[epoch])
all_labels(dat::DataFrame)::Vector{Symbol} = propertynames(dat)


"""
    meta_labels(dat::EegData) -> Vector{Symbol}

Get metadata column names from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of metadata column names
"""
meta_labels(dat::EegData)::Vector{Symbol} = get_cols_by_group(dat, :metadata)


"""
    meta_data(eeg_data::EegData) -> DataFrame

Get meta data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing EEG channel columns
"""
meta_data(dat::SingleDataFrameEeg)::DataFrame = dat.data[:, get_cols_by_group(dat, :metadata)]
meta_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame = dat.data[epoch][:, get_cols_by_group(dat, :metadata)]
meta_data(dat::MultiDataFrameEeg)::DataFrame = to_data_frame(dat)[:, get_cols_by_group(dat, :metadata)]


"""
    channel_data(eeg_data::EegData) -> DataFrame

Get EEG channel data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing EEG channel columns
"""
channel_labels(dat::EegData)::Vector{Symbol} = get_cols_by_group(dat, :channels)
channel_labels(dat::EegData, channel_numbers::Vector{<:UnitRange})::Vector{Symbol} =
    channel_labels(dat)[channel_numbers...]
channel_labels(dat::EegData, channel_numbers)::Vector{Symbol} = channel_labels(dat)[channel_numbers]
channel_labels(dat::EegData, channel_numbers::Int)::Vector{Symbol} = channel_labels(dat)[[channel_numbers]]


"""
    channel_data(eeg_data::EegData) -> DataFrame

Get EEG channel data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing EEG channel columns
"""
channel_data(dat::SingleDataFrameEeg)::DataFrame = dat.data[:, get_cols_by_group(dat, :channels)]
channel_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame = dat.data[epoch][:, get_cols_by_group(dat, :channels)]
channel_data(dat::MultiDataFrameEeg)::DataFrame = to_data_frame(dat)[:, get_cols_by_group(dat, :channels)]

"""

    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing extra/derived columns
"""
extra_labels(dat::EegData)::Vector{Symbol} = get_cols_by_group(dat, :extra)



"""
    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing extra/derived columns
"""
extra_data(dat::SingleDataFrameEeg)::DataFrame = dat.data[:, get_cols_by_group(dat, :extra)]
extra_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame = dat.data[epoch][:, get_cols_by_group(dat, :extra)]
extra_data(dat::MultiDataFrameEeg)::DataFrame = to_data_frame(dat)[:, get_cols_by_group(dat, :extra)]



# === EEG CONVENIENCE FUNCTIONS ===
# Basic information functions
"""
    sample_rate(dat::EegData) -> Int

Get the sample rate of the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Sample rate in Hz
"""
sample_rate(dat::EegData)::Int = dat.sample_rate
sample_rate(dat::DataFrame)::Int = Int(1 / mean(diff(dat.time)))

"""
    reference(dat::EegData) -> String

Get the reference information from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `String`: Reference information
"""
reference(dat::EegData)::Symbol = dat.analysis_info.reference
reference(dat::AnalysisInfo)::Symbol = dat.reference

"""
    filter_info(dat::AnalysisInfo) -> Vector

Get filter information from the analysis info.

# Arguments
- `dat::AnalysisInfo`: The analysis info object

# Returns
- `Vector`: Filter information [hp_filter, lp_filter]
"""
filter_info(dat::AnalysisInfo)::Vector{Float64} = [dat.hp_filter, dat.lp_filter]

# Unique convenience functions (not available through metadata)
"""
    n_samples(dat::EegData) -> Int

Get the number of samples in the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Number of samples
"""
n_samples(dat::SingleDataFrameEeg)::Int = nrow(dat.data)
n_samples(dat::MultiDataFrameEeg)::Int = nrow(dat.data[1])
n_samples(dat::MultiDataFrameEeg, epoch::Int)::Int = nrow(dat.data[epoch])
n_samples(dat::DataFrame)::Int = nrow(dat)

"""
    n_channels(dat::EegData) -> Int

Get the number of channels in the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Number of channels
"""
n_channels(dat::EegData)::Int = length(channel_labels(dat))
n_channels(dat::DataFrame)::Int = length(channel_labels(dat))
n_layout(layout::Layout)::Int = nrow(layout.data)

"""
    n_epochs(dat::EegData) -> Int

Get the number of epochs in the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Int`: Number of epochs
"""
n_epochs(dat::SingleDataFrameEeg)::Int = 1
n_epochs(dat::MultiDataFrameEeg)::Int = length(dat.data)
n_epochs(dat::ErpData)::Int = dat.n_epochs

"""
    duration(dat::EegData) -> Float64

Get the duration of the EEG data in seconds.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Float64`: Duration in seconds
"""
duration(dat::SingleDataFrameEeg) = isempty(dat.data.time) ? 0.0 : last(dat.data.time) - first(dat.data.time)
duration(dat::MultiDataFrameEeg) = isempty(dat.data[1].time) ? 0.0 : last(dat.data[1].time) - first(dat.data[1].time)
duration(dat::MultiDataFrameEeg, epoch::Int) =
    isempty(dat.data[epoch].time) ? 0.0 : last(dat.data[epoch].time) - first(dat.data[epoch].time)


"""
    has_channels(dat::EegData, chans::Vector{Symbol}) -> Bool

Check if the EEG data contains all specified channels.

# Arguments
- `dat::EegData`: The EEG data object
- `chans::Vector{Symbol}`: Vector of channel symbols to check

# Returns
- `Bool`: True if all channels are present
"""
has_channels(dat::EegData, chans::Vector{Symbol})::Bool = all(in(channel_labels(dat)), chans)

"""
    common_channels(dat1::EegData, dat2::EegData) -> Vector{Symbol}

Find common channels between two EEG data objects.

# Arguments
- `dat1::EegData`: First EEG data object
- `dat2::EegData`: Second EEG data object

# Returns
- `Vector{Symbol}`: Vector of common channel symbols
"""
common_channels(dat1::EegData, dat2::EegData)::Vector{Symbol} = intersect(channel_labels(dat1), channel_labels(dat2))

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

function head(dat::EegData; n = nothing)
    isnothing(n) && (n=5)
    result = all_data(dat)[1:n, :]
    viewer(result)
    return result
end

function tail(dat::EegData; n = nothing)
    isnothing(n) && (n=5)
    result = all_data(dat)[(end-n+1):end, :]
    viewer(result)
    return result
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
"""
datarange(x::AbstractVector)::Float64 = -(-(extrema(x)...))


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
colmeans(df::DataFrame, cols)::Vector{Float64} = reduce(+, eachcol(df[!, cols])) ./ length(cols)
colmeans(df::Matrix)::Vector{Float64} = reduce(+, eachcol(df)) ./ size(df)[2]
colmeans(df::Matrix, cols)::Vector{Float64} = reduce(+, eachcol(df[:, cols])) ./ length(cols)


"""
    data_limits_x(dat::DataFrame) -> Union{Tuple{Float64,Float64}, Nothing}

Get the time range of the data.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum time values
- `Nothing`: If the DataFrame is empty
"""
function data_limits_x(dat::DataFrame; col::Symbol = :time)
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
function data_limits_y(dat::DataFrame, col::Symbol)
    isempty(dat) && return nothing
    return [minimum(dat[!, col]), maximum(dat[!, col])]
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
"""
function subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame
    return df[selected_samples, selected_channels]
end


function subset(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)

    @info "subset: Subsetting $(typeof(dat)) ..."
    # Get selected channels/samples 
    selected_channels = get_selected_channels(dat, channel_selection, include_extra = include_extra)
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
"""
function subset(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::EpochData

    @info "subset: Subsetting $(typeof(dat)) ..."

    # Get selected channels, samples, and epochs
    selected_channels = get_selected_channels(dat, channel_selection, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    selected_epochs = get_selected_epochs(dat, epoch_selection)

    # Filter epochs first, then apply channel/sample filtering and match layout
    epochs_subset = subset_dataframe.(dat.data[selected_epochs], Ref(selected_channels), Ref(selected_samples))
    layout_subset = subset_layout(dat.layout, channel_selection = channel_selection)

    # Create new EpochData object
    return EpochData(epochs_subset, layout_subset, dat.sample_rate, dat.analysis_info)

end
