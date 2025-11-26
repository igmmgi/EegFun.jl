#  === EEG DATA UTILITIES ===
#
# This file provides utilities for accessing and manipulating EEG data structures.
# It includes functions for column identification, data access, convenience functions,
# data viewing, mathematical utilities, and subsetting operations.
#
# Organization:
# - Column Identification: Functions to identify metadata, channel, and extra columns
# - Data Access: Functions to extract different types of data (all_data, meta_data, etc.)
# - Convenience Functions: Basic information and size functions
# - Data Viewing: head, tail, viewer functions
# - Mathematical Utilities: datarange, colmeans, data_limits
# - Subsetting: Functions to subset EEG data objects
#

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
function get_cols_by_group(dat::EegData, group::Symbol)

    if !(group in [:channels, :metadata, :extra])
        @minimal_error "Unknown group type: $group"
    end

    labels = all_labels(dat)
    layout_channels = dat.layout.data.label

    if group == :channels
        return intersect(layout_channels, labels)
    elseif group == :metadata
        isempty(layout_channels) && return Symbol[]
        first_channel_idx = findfirst(col -> col == layout_channels[1], labels)
        isnothing(first_channel_idx) && @minimal_error "First channel label not found in data"
        return labels[1:(first_channel_idx-1)]
    elseif group == :extra
        isempty(layout_channels) && return Symbol[]
        last_channel_idx = findlast(col -> col == layout_channels[end], labels)
        if isnothing(last_channel_idx)
            @warn "Last channel label $(layout_channels[end]) not found in data columns. Available columns: $(labels)"
            return Symbol[]  # Return empty array instead of throwing error
        end
        return labels[(last_channel_idx+1):end]
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
all_data(dat::Vector{<:MultiDataFrameEeg})::DataFrame = to_data_frame(dat) # single data frame with all epochs from all objects

function all_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    return to_data_frame(subset(dat, epoch_selection = epoch_selection))
end

function meta_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    meta_cols = get_cols_by_group(dat, :metadata)
    return isempty(meta_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[:, meta_cols]
end

function channel_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    channel_cols = get_cols_by_group(dat, :channels)
    return isempty(channel_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[:, channel_cols]
end

function extra_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    extra_cols = get_cols_by_group(dat, :extra)
    return isempty(extra_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[:, extra_cols]
end

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

"""
    all_labels(dat::DataFrame) -> Vector{Symbol}

Get all column names from a DataFrame.

# Arguments
- `dat::DataFrame`: The DataFrame

# Returns
- `Vector{Symbol}`: All column names
"""
all_labels(dat::DataFrame)::Vector{Symbol} = propertynames(dat)


"""
    meta_labels(dat::EegData) -> Vector{Symbol}

Get metadata column names from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of metadata column names
"""
meta_labels(dat::EegData) = get_cols_by_group(dat, :metadata)


"""
    meta_data(eeg_data::EegData) -> DataFrame

Get meta data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing metadata columns
"""
function meta_data(dat::SingleDataFrameEeg)
    meta_cols = get_cols_by_group(dat, :metadata)
    return isempty(meta_cols) ? DataFrame() : dat.data[:, meta_cols]
end

function meta_data(dat::MultiDataFrameEeg, epoch::Int)
    meta_cols = get_cols_by_group(dat, :metadata)
    return isempty(meta_cols) ? DataFrame() : dat.data[epoch][:, meta_cols]
end

function meta_data(dat::MultiDataFrameEeg)
    meta_cols = get_cols_by_group(dat, :metadata)
    return isempty(meta_cols) ? DataFrame() : to_data_frame(dat)[:, meta_cols]
end


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

# Handle collections of EEG data
# TODO: is it possible that all do not have the same channel labels?
channel_labels(dat::Vector{<:EegData})::Vector{Symbol} = channel_labels(first(dat))


"""
    channel_data(eeg_data::EegData) -> DataFrame

Get EEG channel data columns from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing EEG channel columns
"""
function channel_data(dat::SingleDataFrameEeg)::DataFrame
    channel_cols = get_cols_by_group(dat, :channels)
    return isempty(channel_cols) ? DataFrame() : dat.data[:, channel_cols]
end

function channel_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame
    channel_cols = get_cols_by_group(dat, :channels)
    return isempty(channel_cols) ? DataFrame() : dat.data[epoch][:, channel_cols]
end

function channel_data(dat::MultiDataFrameEeg)::DataFrame
    channel_cols = get_cols_by_group(dat, :channels)
    return isempty(channel_cols) ? DataFrame() : to_data_frame(dat)[:, channel_cols]
end

"""

    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing extra/derived columns
"""
extra_labels(dat::EegData) = get_cols_by_group(dat, :extra)

# Handle collections of EEG data
# TODO: is it possible that all do not have the same channel labels?
extra_labels(dat::Vector{<:EegData})::Vector{Symbol} = extra_labels(first(dat))



"""
    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.

# Arguments
- `eeg_data::EegData`: The EEG data object

# Returns
- `DataFrame`: DataFrame containing extra/derived columns
"""
function extra_data(dat::SingleDataFrameEeg)::DataFrame
    extra_cols = get_cols_by_group(dat, :extra)
    return isempty(extra_cols) ? DataFrame() : dat.data[:, extra_cols]
end

function extra_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame
    extra_cols = get_cols_by_group(dat, :extra)
    return isempty(extra_cols) ? DataFrame() : dat.data[epoch][:, extra_cols]
end

function extra_data(dat::MultiDataFrameEeg)::DataFrame
    extra_cols = get_cols_by_group(dat, :extra)
    return isempty(extra_cols) ? DataFrame() : to_data_frame(dat)[:, extra_cols]
end



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
    reference(dat::EegData) -> Symbol

Get the reference information from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Symbol`: Reference information
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

"""
    n_layout(layout::Layout) -> Int

Get the number of channels in the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Int`: Number of channels in the layout
"""
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


condition_number(dat::ContinuousData)::String = "Raw Data"
condition_number(dat::ErpData)::Int = dat.condition
condition_number(dat::EpochData)::Int = dat.condition

condition_name(dat::ContinuousData)::String = "Raw Data"
condition_name(dat::ErpData)::String = dat.condition_name
condition_name(dat::EpochData)::String = dat.condition_name

file_name(dat::EpochData)::String = dat.file


"""
    duration(dat::EegData) -> Float64

Get the duration of the EEG data in seconds.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Float64`: Duration in seconds
"""
duration(dat::SingleDataFrameEeg)::Float64 =
    hasproperty(dat.data, :time) && !isempty(dat.data.time) ? last(dat.data.time) - first(dat.data.time) : 0.0
duration(dat::MultiDataFrameEeg)::Float64 =
    hasproperty(dat.data[1], :time) && !isempty(dat.data[1].time) ? last(dat.data[1].time) - first(dat.data[1].time) :
    0.0
duration(dat::MultiDataFrameEeg, epoch::Int)::Float64 =
    hasproperty(dat.data[epoch], :time) && !isempty(dat.data[epoch].time) ?
    last(dat.data[epoch].time) - first(dat.data[epoch].time) : 0.0


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

"""
    viewer(dat)

Display data using VS Code viewer if available, otherwise use standard display.

# Arguments
- `dat`: Any data object to display
"""
function viewer(dat)
    if get(ENV, "TERM_PROGRAM", "") == "vscode"
        try
            Main.vscodedisplay(dat)
        catch
            display(dat)
        end
    else
        display(dat)
    end
end

"""
    viewer(dat::EegData)

Display EEG data by extracting all data and using the appropriate viewer.

# Arguments
- `dat::EegData`: The EEG data object to display
"""
function viewer(dat::EegData)
    viewer(all_data(dat))
end

function head(dat::EegData; n = nothing)
    isnothing(n) && (n = 5)
    data = all_data(dat)
    nrows = nrow(data)
    n = min(n, nrows)  # Don't exceed available rows
    result = n > 0 ? data[1:n, :] : DataFrame()
    viewer(result)
    return result
end

function tail(dat::EegData; n = nothing)
    isnothing(n) && (n = 5)
    data = all_data(dat)
    nrows = nrow(data)
    n = min(n, nrows)  # Don't exceed available rows
    result = n > 0 ? data[max(1, nrows-n+1):nrows, :] : DataFrame()
    viewer(result)
    return result
end

"""
    to_data_frame(dat::EpochData) -> DataFrame

Convert EpochData to a single DataFrame by concatenating all epochs.

# Arguments
- `dat::EpochData`: The epoch data to convert

# Returns
- `DataFrame`: Single DataFrame with all epochs concatenated vertically
"""
function to_data_frame(dat::EpochData)
    isempty(dat.data) && return DataFrame()
    return vcat(dat.data...)
end

"""
    to_data_frame(dat::Vector{EpochData}) -> DataFrame

Convert a vector of EpochData objects to a single DataFrame.

# Arguments
- `dat::Vector{EpochData}`: Vector of epoch data objects to convert

# Returns
- `DataFrame`: Single DataFrame with all epochs from all objects concatenated
"""
function to_data_frame(dat::Vector{EpochData})
    isempty(dat) && return DataFrame()
    isempty(dat[1].data) && return DataFrame()
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end

# === MATHEMATICAL UTILITIES ===

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

# Arguments
- `dat::DataFrame`: The DataFrame containing the data
- `col::Symbol`: The column to get limits for

# Returns
- `Vector{Float64}`: [minimum, maximum] across specified column
- `Nothing`: If the DataFrame is empty
"""
function data_limits_y(dat::DataFrame, col::Symbol)
    isempty(dat) && return nothing
    mn, mx = extrema(dat[!, col])
    return [mn, mx]
end

"""
    data_limits_y(dat::DataFrame, cols::Vector{Symbol}) -> Union{Vector{Float64}, Nothing}

Get the value range across multiple specified columns.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data
- `cols::Vector{Symbol}`: The columns to get limits for

# Returns
- `Vector{Float64}`: [minimum, maximum] across all specified columns
- `Nothing`: If the DataFrame is empty
"""
function data_limits_y(dat::DataFrame, cols::Vector{Symbol})
    isempty(dat) && return nothing
    global_min, global_max = Inf, -Inf
    @inbounds for col in cols
        mn, mx = extrema(dat[!, col])
        global_min = min(global_min, mn)
        global_max = max(global_max, mx)
    end
    return [global_min, global_max]
end



# === DATAFRAME SUBSETTING UTILITIES ===

"""
    subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> DataFrame

Create a subset of a DataFrame by selecting specific channels and samples.

# Arguments
- `df::DataFrame`: The DataFrame to subset
- `selected_channels::Vector{Symbol}`: Column names to include
- `selected_samples::Vector{Int}`: Row indices to include

# Returns
- `DataFrame`: Subset DataFrame with selected channels and samples
"""
function subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame
    return df[selected_samples, selected_channels]
end


# === DEFAULT Y-RANGE HELPERS ===
"""
    yrange(dat::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1)

Compute a padded y-range for ERP data using channel selection predicate.
Returns (min, max) padded by `buffer` proportion.
"""
function ylimits(
    dat::ErpData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}
    dat_sub = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = include_extra,
    )
    chs = channel_labels(dat_sub)
    lims = data_limits_y(dat_sub.data, chs)
    return (lims[1], lims[2])
end

"""
    yrange(dat::EpochData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1, average_only::Bool = false)

Compute a padded y-range for epoch data using channel selection predicate.
If `average_only=true`, uses the averaged waveform across epochs.
"""
function ylimits(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}

    # Apply predicates via subset first
    dat_sub = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = include_extra,
    )
    # Determine which value columns to use
    chs = channel_labels(dat_sub)
    # Compute limits per epoch and combine
    limits = map(df -> data_limits_y(df, chs), dat_sub.data)
    limits = Base.filter(!isnothing, limits)
    isempty(limits) && return (0.0, 1.0)
    min_val = minimum(lim[1] for lim in limits)
    max_val = maximum(lim[2] for lim in limits)
    return (min_val, max_val)
end




"""
    subset_dataframes(dataframes::Vector{DataFrame}, selected_epochs::Vector{Int}, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> Vector{DataFrame}

Create subsets of multiple DataFrames by selecting specific epochs, channels, and samples.

# Arguments
- `dataframes::Vector{DataFrame}`: Vector of DataFrames to subset
- `selected_epochs::Vector{Int}`: Indices of epochs (DataFrames) to include
- `selected_channels::Vector{Symbol}`: Column names to include in each DataFrame
- `selected_samples::Vector{Int}`: Row indices to include in each DataFrame

# Returns
- `Vector{DataFrame}`: Vector of subset DataFrames
"""
function subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame}
    return subset_dataframe.(dataframes[selected_epochs], Ref(selected_channels), Ref(selected_samples))
end




# === INTERNAL SUBSET HELPERS ===

"""
    _subset_common(dat, channel_selection, sample_selection, include_extra)

Internal helper for common subsetting operations on EEG data.

# Arguments
- `dat`: EEG data object
- `channel_selection::Function`: Channel selection predicate
- `sample_selection::Function`: Sample selection predicate
- `include_extra::Bool`: Whether to include extra columns

# Returns
- `Tuple`: (selected_channels, selected_samples, layout_subset)
"""
function _subset_common(dat, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting channels and samples"
    # Get subset selected channels, samples, and layout
    selected_channels = get_selected_channels(dat, channel_selection, include_extra = include_extra)
    selected_samples = get_selected_samples(dat, sample_selection)
    layout_subset = subset_layout(dat.layout, channel_selection = channel_selection)

    return selected_channels, selected_samples, layout_subset
end

"""
    _subset_common(dat::EpochData, epoch_selection, channel_selection, sample_selection, include_extra)

Internal helper for subsetting EpochData with epoch selection.

# Arguments
- `dat::EpochData`: Epoch data object
- `epoch_selection::Function`: Epoch selection predicate
- `channel_selection::Function`: Channel selection predicate
- `sample_selection::Function`: Sample selection predicate
- `include_extra::Bool`: Whether to include extra columns

# Returns
- `Tuple`: (selected_epochs, selected_channels, selected_samples, layout_subset)
"""
function _subset_common(dat::EpochData, epoch_selection, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting epochs, channels and samples"
    # Get subset selected epochs, channels, samples, and layout
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    selected_channels, selected_samples, layout_subset =
        _subset_common(dat, channel_selection, sample_selection, include_extra)
    return selected_epochs, selected_channels, selected_samples, layout_subset
end

"""
    _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info) -> ContinuousData

Internal helper to create ContinuousData from subset DataFrame.
"""
function _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info, file::String)
    return ContinuousData(file, data_subset, layout, sample_rate, analysis_info)
end

"""
    _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info, n_epochs::Int) -> ErpData

Internal helper to create ErpData from subset DataFrame.
"""
function _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info, n_epochs::Int, condition::Int, condition_name::String, file::String)
    return ErpData(file, condition, condition_name, data_subset, layout, sample_rate, analysis_info, n_epochs)
end

"""
    _create_subset(data_subset::Vector{DataFrame}, layout, sample_rate::Int, analysis_info) -> EpochData

Internal helper to create EpochData from subset DataFrames.
"""
function _create_subset(data_subset::Vector{DataFrame}, layout, sample_rate::Int, analysis_info, condition::Int, condition_name::String, file::String)
    return EpochData(file, condition, condition_name, data_subset, layout, sample_rate, analysis_info)
end

# === SUBSET IMPLEMENTATIONS ===

function subset(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::SingleDataFrameEeg
    selected_channels, selected_samples, layout_subset =
        _subset_common(dat, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info, dat.file)
end

function subset(
    dat::ErpData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::ErpData
    selected_channels, selected_samples, layout_subset =
        _subset_common(dat, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info, dat.n_epochs, dat.condition, dat.condition_name, dat.file)
end

function subset(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::EpochData
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframes(dat.data, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info, dat.condition, dat.condition_name, dat.file)
end

function subset(
    datasets::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::Vector{ErpData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]
    
    # Then apply channel and sample selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::Vector{EpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::Vector{EpochData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]
    
    # Then apply channel, sample, and epoch selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
end



"""
    log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

Log a pretty table with specified log level. For general DataFrame logging.
Sets show_row_number=false and show_subheader=false by default for cleaner logs.

# Arguments
- `df::DataFrame`: The DataFrame to log
- `log_level::Symbol`: Log level (:debug, :info, :warn, :error) (default: :info)
- `kwargs...`: Additional arguments passed to pretty_table

# Examples
```julia
# Log with default info level
log_pretty_table(df; title = "My Table")

# Log with debug level
log_pretty_table(df; log_level = :debug, title = "Debug Table")

# Log with warn level
log_pretty_table(df; log_level = :warn, title = "Warning Table")
```
"""
function log_pretty_table(df::DataFrame; log_level::Symbol = :info, kwargs...)

    table_output = sprint() do output_io
        # TODO: better way of doing this?
        # Set a large display size to avoid terminal limitations (i.e., cropping!)
        io_context = IOContext(output_io, :displaysize => (2000, 2000))
        pretty_table(io_context, df; kwargs...)
    end
    
    # Log with specified level
    if log_level == :debug
        @debug "\n\n$table_output\n"
    elseif log_level == :info
        @info "\n\n$table_output\n"
    elseif log_level == :warn
        @warn "\n\n$table_output\n"
    elseif log_level == :error
        @error "\n\n$table_output\n"
    else
        @warn "Unknown log level: $log_level, using :info instead"
        @info "\n\n$table_output\n"
    end
    
    return nothing
end


# Helper function predicates for easier channel filtering
channels() = x -> fill(true, length(x))  # Default: select all channels given
channels(channel_names::Vector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_name::Symbol) = x -> x .== channel_name
channels(channel_number::Int) = channels([channel_number])
channels(channel_numbers::Union{Vector{Int},UnitRange}) = x -> [i in channel_numbers for i = 1:length(x)]
channels(channel_ranges::Vector{UnitRange{Int}}) = x -> [i in union(channel_ranges...) for i = 1:length(x)]
channels_not(channel_names::Vector{Symbol}) = x -> .!(x .∈ Ref(channel_names))
channels_not(channel_name::Symbol) = x -> .!(x .== channel_name)
channels_not(channel_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in channel_numbers for i = 1:length(x)])
channels_not(mixed::Vector) = x -> begin
    # Handle mixed Int and UnitRange{Int} (e.g., [-2, 1:10])
    combined = Set{Int}()
    for item in mixed
        if item isa Int
            push!(combined, item)
        elseif item isa UnitRange{Int}
            union!(combined, item)
        else
            throw(ArgumentError("channels_not() only accepts Int or UnitRange{Int}, got $(typeof(item))"))
        end
    end
    .!([i in combined for i = 1:length(x)])
end

# Helper function predicates for easier component filtering
components() = x -> fill(true, length(x))  # Default: select all components given
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> [i in component_numbers for i = 1:length(x)]
components(component_number::Int) = x -> x .== component_number
components_not(component_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in component_numbers for i = 1:length(x)])
components_not(component_number::Int) = x -> .!(x .== component_number)

# Helper function predicates for easier sample filtering
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
samples_or(columns::Vector{Symbol}) = x -> any(x[!, col] for col in columns)
samples_and(columns::Vector{Symbol}) = x -> all(x[!, col] for col in columns)
samples_not(column::Symbol) = x -> .!(x[!, column])
samples_or_not(columns::Vector{Symbol}) = x -> .!(any(x[!, col] for col in columns))
samples_and_not(columns::Vector{Symbol}) = x -> .!(all(x[!, col] for col in columns))

# Helper function predicates for easier epoch filtering
epochs() = x -> fill(true, length(x))  # Default: select all epochs given
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> [i in epoch_numbers for i in x]
epochs(epoch_number::Int) = x -> x .== epoch_number
epochs_not(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in epoch_numbers for i in x])
epochs_not(epoch_number::Int) = x -> .!(x .== epoch_number)

# Helper to extract condition name from ErpData or EpochData
_get_condition_name(dat::ErpData)::String = dat.condition_name
_get_condition_name(dat::EpochData)::String = dat.condition_name

# Helper function predicates for easier condition filtering (for Vector{ErpData} and Vector{EpochData})
conditions() = x -> fill(true, length(x))  # Default: select all conditions given
conditions(condition_indices::Union{Vector{Int},UnitRange}) = x -> [i in condition_indices for i = 1:length(x)]
conditions(condition_index::Int) = x -> [i == condition_index for i = 1:length(x)]
conditions(condition_names::Vector{String}) = x -> [_get_condition_name(dat) in condition_names for dat in x]
conditions(condition_name::String) = x -> [_get_condition_name(dat) == condition_name for dat in x]
conditions_not(condition_indices::Union{Vector{Int},UnitRange}) = x -> .!([i in condition_indices for i = 1:length(x)])
conditions_not(condition_index::Int) = x -> .!([i == condition_index for i = 1:length(x)])
conditions_not(condition_names::Vector{String}) = x -> .!([_get_condition_name(dat) in condition_names for dat in x])
conditions_not(condition_name::String) = x -> .!([_get_condition_name(dat) == condition_name for dat in x])

# Helper to select channels/columns based on a predicate (+ which to include)
function get_selected_channels(dat, channel_selection::Function; include_meta::Bool = true, include_extra::Bool = true)

    # Columns/channels in dataframe to include
    metadata_cols = include_meta ? meta_labels(dat) : Symbol[]
    selectable_cols = include_extra ? vcat(channel_labels(dat), extra_labels(dat)) : channel_labels(dat)

    # Apply channel selection to non-metadata columns
    selection_mask = channel_selection(selectable_cols)
    selected = selectable_cols[selection_mask]
    
    # Preserve user-specified order if available
    selection_type = typeof(channel_selection)
    if hasfield(selection_type, :channel_names)
        # channels([Symbol...]) or channels_not([Symbol...]) - validate and optionally preserve order
        user_order = getfield(channel_selection, :channel_names)
        existing_in_order = [ch for ch in user_order if ch in selectable_cols]
        existing_in_selected = [ch for ch in existing_in_order if ch in selected]
        
        # Validate: check for missing channels and duplicates
        seen = Set{Symbol}()
        for ch in user_order
            if ch ∉ selectable_cols
                @minimal_warning "Channel $(ch) not found in data!"
            elseif ch ∈ seen
                @minimal_warning "Channel $(ch) already specified!"
            else
                push!(seen, ch)
            end
        end
        
        # If all existing channels from user_order are in selected, it's channels([...]) - preserve order
        if !isempty(existing_in_order) && length(existing_in_selected) == length(existing_in_order)
            selected = Symbol[]
            for ch in user_order
                if ch ∉ selectable_cols || ch ∈ selected
                    continue
                end
                push!(selected, ch)
            end
        end
        # Otherwise it's channels_not([...]) - use default order from selection_mask
        
    elseif hasfield(selection_type, :channel_numbers)
        # channels([Int...]) or channels_not([Int...]) - validate and optionally preserve order
        user_order_numbers = getfield(channel_selection, :channel_numbers)
        
        # Validate: check for invalid indices and duplicates
        seen = Set{Int}()
        for i in user_order_numbers
            if i < 1 || i > length(selectable_cols)
                @minimal_warning "Channel index $(i) out of range (valid: 1:$(length(selectable_cols)))!"
            elseif i ∈ seen
                @minimal_warning "Channel index $(i) already specified!"
            else
                push!(seen, i)
            end
        end
        
        # Check if this is channels([...]) or channels_not([...])
        # For channels([...]), the indices should be in selected
        # For channels_not([...]), the indices should NOT be in selected
        valid_indices = [i for i in user_order_numbers if 1 <= i <= length(selectable_cols)]
        selected_indices = [i for i in valid_indices if selectable_cols[i] in selected]
        
        # If all valid indices are in selected, it's channels([...]) - preserve order
        if !isempty(valid_indices) && length(selected_indices) == length(valid_indices)
            selected = Symbol[]
            for i in user_order_numbers
                if i < 1 || i > length(selectable_cols) || selectable_cols[i] ∈ selected
                    continue
                end
                push!(selected, selectable_cols[i])
            end
        end
        # Otherwise it's channels_not([...]) - use default order from selection_mask
    end

    # Return metadata + selected channels
    return vcat(metadata_cols, selected)
end


# Helper to select components based on a predicate
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end

# Helper to select samples based on a predicate
function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data))
end

# Helper to select samples based on a predicate
function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data[1])) # assume all data have the same samples
end

# Helper to select samples from a DataFrame
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    return findall(sample_selection(dat))
end

# Helper to select epochs based on a predicate
function get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection::Function)
    all_epochs = 1:length(dat.data)
    return findall(epoch_selection(all_epochs))
end

# Helper to select conditions from Vector{ErpData} based on a predicate
function get_selected_conditions(datasets::Vector{ErpData}, condition_selection::Function)
    all_indices = 1:length(datasets)
    return findall(condition_selection(datasets))
end

# Helper to select conditions from Vector{EpochData} based on a predicate
function get_selected_conditions(datasets::Vector{EpochData}, condition_selection::Function)
    all_indices = 1:length(datasets)
    return findall(condition_selection(datasets))
end




"""
    convert(dat::MultiDataFrameEeg, epoch_idx::Int) -> SingleDataFrameEeg

Convert a single epoch from MultiDataFrameEeg to SingleDataFrameEeg.

# Arguments
- `dat::MultiDataFrameEeg`: The multi-DataFrame EEG data
- `epoch_idx::Int`: Index of the epoch to convert (1-based)

# Returns
- `SingleDataFrameEeg`: Single DataFrame containing only the specified epoch

# Examples
```julia
# Convert epoch 3 to single DataFrame
single_dat = convert(dat, 3)

# Now you can use single DataFrame functions
summary = channel_summary(single_dat)
```
"""
function convert(dat::T, epoch_idx::Int)::ContinuousData where {T<:MultiDataFrameEeg}
    # Validate epoch index
    if epoch_idx < 1 || epoch_idx > length(dat.data)
        @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    end
    return ContinuousData(dat.file, dat.data[epoch_idx], dat.layout, dat.sample_rate, dat.analysis_info) # TODO: should we use SingleDataFrameEeg?
end


# === CHANNEL RENAMING UTILITIES ===
"""
    rename_channel!(dat::EegData, rename_dict::Dict{Symbol, Symbol})

Rename channels in EEG data using a dictionary mapping old names to new names.
Modifies the data in place by updating both the data columns and the layout.

# Arguments
- `dat::EegData`: The EEG data object to modify
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `nothing` (modifies the data in place)

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
rename_channel!(dat, rename_dict)

# Rename a single channel
rename_channel!(dat, Dict(:Cz => :Cz_new))
```

# Notes
- Only channels that exist in the data will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- Updates both the data columns and the layout labels
- Clears any cached neighbour information in the layout since channel names have changed
- Properly handles swaps (e.g., Dict(:A => :B, :B => :A) correctly exchanges the channels)
"""
function rename_channel!(dat::EegData, rename_dict::Dict{Symbol,Symbol})
    # Capture original channel names before renaming layout
    original_channels = Set(dat.layout.data.label)

    # First rename channels in the layout
    rename_channel!(dat.layout, rename_dict)

    # Get the list of channels that were actually renamed in the layout
    layout_channels = dat.layout.data.label
    channels_to_rename = keys(rename_dict)
    channels_found = intersect(original_channels, channels_to_rename)

    if isempty(channels_found)
        @info "rename_channel!: No channels found to rename in data"
        return nothing
    end

    # Now rename the corresponding data columns using multiple dispatch
    _rename_data_columns!(dat, rename_dict, original_channels)

    @info "rename_channel!: Renamed $(length(channels_found)) channels in data and layout"
    return nothing
end

# Multiple dispatch for different EEG data types
function _rename_data_columns!(df::DataFrame, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    # Check for potential duplicate names before applying any renames
    final_names = Symbol[]
    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels && old_name ∈ propertynames(df)
            push!(final_names, new_name)
        end
    end

    # Check for duplicates in final names
    if length(final_names) != length(unique(final_names))
        duplicate_names = Base.filter(x -> count(==(x), final_names) > 1, unique(final_names))
        @minimal_error_throw "Cannot rename channels to duplicate names: $(join(duplicate_names, ", "))"
    end

    # Apply the renaming with proper swap handling
    # First, collect all the final rename mappings to avoid interference
    final_renames = Dict{Symbol,Symbol}()  # old_name => final_name

    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels && old_name ∈ propertynames(df)
            final_renames[old_name] = new_name
        end
    end

    # Nothing to do
    isempty(final_renames) && return

    # Two-phase rename to avoid collisions (e.g., swaps)
    # Phase 1: rename old names to unique temporary names
    temp_renames = Pair{Symbol,Symbol}[]
    used_names = Set(propertynames(df))  # currently present in df
    union!(used_names, Set(values(final_renames)))  # also avoid targeting final names

    for (old_name, new_name) in final_renames
        # Propose a unique temporary name
        base_tmp = Symbol(string(new_name), "__tmp__")
        tmp = base_tmp
        counter = 1
        while tmp ∈ used_names
            tmp = Symbol(string(base_tmp), "_", counter)
            counter += 1
        end
        push!(temp_renames, old_name => tmp)
        push!(used_names, tmp)
    end

    # Execute phase 1 renames
    for p in temp_renames
        rename!(df, p)
    end

    # Phase 2: rename temporary names to final names
    for (old_name, new_name) in final_renames
        # Find the temporary assigned to this old_name
        tmp = only(last.(Base.filter(p -> first(p) == old_name, temp_renames)))
        rename!(df, tmp => new_name)
    end
end

function _rename_data_columns!(
    dat::SingleDataFrameEeg,
    rename_dict::Dict{Symbol,Symbol},
    existing_channels::Set{Symbol},
)
    _rename_data_columns!(dat.data, rename_dict, existing_channels)
end

function _rename_data_columns!(dat::MultiDataFrameEeg, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
    # Use broadcasting to apply the DataFrame method to all DataFrames
    _rename_data_columns!.(dat.data, Ref(rename_dict), Ref(existing_channels))
end




"""
    rename_channel(dat::EegData, rename_dict::Dict{Symbol, Symbol})

Create a renamed copy of EEG data using a dictionary mapping old names to new names.

# Arguments
- `dat::EegData`: The EEG data object to rename
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `EegData`: A new EEG data object with renamed channels

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
new_dat = rename_channel(dat, rename_dict)
```

# Notes
- Only channels that exist in the data will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- Updates both the data columns and the layout labels
- The original data is not modified
- Properly handles swaps (e.g., Dict(:A => :B, :B => :A) correctly exchanges the channels)
"""
function rename_channel(dat::EegData, rename_dict::Dict{Symbol,Symbol})
    # Create a copy of the data using the existing copy function
    new_dat = copy(dat)
    rename_channel!(new_dat, rename_dict)
    return new_dat
end


# =============================================================================
# DATA CREATION FUNCTIONS
# =============================================================================

"""
    create_eeg_dataframe(dat::BiosemiDataFormat.BiosemiData)::DataFrame

Creates a DataFrame from a BiosemiDataFormat data structure.

# Arguments
- `dat::BiosemiDataFormat.BiosemiData`: The BiosemiDataFormat data structure containing EEG data.

# Returns
- `DataFrame`: DataFrame containing the EEG data with time, sample, triggers, and channel columns.

# Examples
```julia
# Create DataFrame from BiosemiDataFormat data
df = create_eeg_dataframe(biosemi_data)
```
"""
function create_eeg_dataframe(dat::BiosemiDataFormat.BiosemiData)::DataFrame
    @info "create_eeg_dataframe: Creating EEG DataFrame"
    df = hcat(
        DataFrame(
            time = dat.time,
            sample = 1:length(dat.time),
            triggers = _clean_triggers(dat.triggers.raw),
        ),
        DataFrame(Float64.(dat.data), Symbol.(dat.header.channel_labels[1:(end-1)])),  # assumes last channel is trigger
    )
    return df
end

"""
    create_eeg_dataframe(dat::BiosemiDataFormat.BiosemiData, layout::Layout)::ContinuousData

Creates a ContinuousData object from a BiosemiDataFormat data structure and a layout.

# Arguments
- `dat::BiosemiDataFormat.BiosemiData`: The BiosemiDataFormat data structure containing EEG data.
- `layout::Layout`: The layout object containing electrode information.

# Returns
- `ContinuousData`: ContinuousData object containing the EEG data and layout information.

# Examples
```julia
# Create ContinuousData from BiosemiDataFormat data and layout
eeg_data = create_eeg_dataframe(biosemi_data, layout)
```
"""
function create_eeg_dataframe(dat::BiosemiDataFormat.BiosemiData, layout::Layout)::ContinuousData
    file_name = filename(dat)
    df = create_eeg_dataframe(dat)
    return ContinuousData(file_name, df, layout, dat.header.sample_rate[1], AnalysisInfo())
end

"""
    create_eeg_dataframe(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame

Creates a DataFrame from a BrainVisionDataFormat data structure.

# Arguments
- `dat::BrainVisionDataFormat.BrainVisionData`: The BrainVisionDataFormat data structure containing EEG data.

# Returns
- `DataFrame`: DataFrame containing the EEG data with time, sample, triggers, and channel columns.

# Examples
```julia
# Create DataFrame from BrainVisionDataFormat data
df = create_eeg_dataframe(brainvision_data)
```
"""
function create_eeg_dataframe(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame
    @info "create_eeg_dataframe: Creating EEG DataFrame from BrainVision data"

    # Check if data is available
    if isnothing(dat.data)
        @minimal_error_throw "BrainVision data is empty (data field is nothing)"
    end

    if isnothing(dat.header)
        @minimal_error_throw "BrainVision header is empty (header field is nothing)"
    end

    # Extract basic information
    # BrainVision data is now samples × channels (after package modification)
    n_samples = size(dat.data, 1)
    n_channels = size(dat.data, 2)
    sample_rate = dat.header.Fs

    # Create time vector
    time = collect(0:(n_samples-1)) ./ sample_rate

    # Create sample vector
    sample = 1:n_samples

    # Extract channel labels from header
    channel_labels = dat.header.label

    # Verify channel count matches
    if length(channel_labels) != n_channels
        @minimal_error_throw "Number of channel labels ($(length(channel_labels))) does not match number of channels ($n_channels)"
    end

    # Create triggers and marker strings columns from markers
    if isnothing(dat.markers) || isempty(dat.markers)
        @info "create_eeg_dataframe: No markers found, creating empty trigger columns"
        triggers = zeros(Int, n_samples)
        triggers_info = fill("", n_samples)
    else
        @info "create_eeg_dataframe: Found $(length(dat.markers)) markers"
        triggers, triggers_info = _extract_triggers_from_markers(dat.markers, n_samples)
    end

    # Create the DataFrame 
    df = hcat(
        DataFrame(
            time = time,
            sample = sample,
            triggers = triggers,
            triggers_info = triggers_info,
        ),
        DataFrame(dat.data, channel_labels),
    );

    return df
end

"""
    create_eeg_dataframe(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData

Creates a ContinuousData object from a BrainVisionDataFormat data structure and a layout.

# Arguments
- `dat::BrainVisionDataFormat.BrainVisionData`: The BrainVisionDataFormat data structure containing EEG data.
- `layout::Layout`: The layout object containing electrode information.

# Returns
- `ContinuousData`: ContinuousData object containing the EEG data and layout information.

# Examples
```julia
# Create ContinuousData from BrainVisionDataFormat data and layout
eeg_data = create_eeg_dataframe(brainvision_data, layout)
```
"""
function create_eeg_dataframe(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData
    file_name = basename_without_ext(dat.filename)
    df = create_eeg_dataframe(dat)
    return ContinuousData(file_name, df, layout, dat.header.Fs, AnalysisInfo())
end

"""
    _extract_triggers_from_markers(markers::Vector{BrainVisionDataFormat.BrainVisionMarker}, n_samples::Int)::Tuple{Vector{Int}, Vector{String}}

Extract trigger values from BrainVision markers and create trigger and marker string vectors.

# Arguments
- `markers::Vector{BrainVisionDataFormat.BrainVisionMarker}`: Vector of BrainVision markers
- `n_samples::Int`: Number of samples in the data

# Returns
- `Tuple{Vector{Int}, Vector{String}}`: (trigger vector, marker string vector) with values at appropriate sample positions
"""
function _extract_triggers_from_markers(
    markers::Vector{BrainVisionDataFormat.BrainVisionMarker},
    n_samples::Int,
)::Tuple{Vector{Int},Vector{String}}
    triggers = zeros(Int, n_samples)
    triggers_info = fill("", n_samples)

    @info "Processing $(length(markers)) markers for $n_samples samples"

    # Debug: Show first few marker sample indices
    if length(markers) > 0
        first_few_samples = [marker.sample for marker in markers[1:min(5, length(markers))]]
        @info "First few marker sample indices: $first_few_samples"
    end

    # First pass: extract all unique trigger values (including empty strings for system markers)
    # Check if markers are 0-based or 1-based by looking at the first marker
    is_zero_based = length(markers) > 0 && markers[1].sample == 0
    @info "Detected $(is_zero_based ? "0-based" : "1-based") indexing for marker samples"

    unique_values = Set{String}()
    valid_markers = 0
    for marker in markers
        # Convert to 1-based if needed
        sample_idx = is_zero_based ? marker.sample + 1 : marker.sample

        if 1 <= sample_idx <= n_samples
            # Include all markers, even those with empty values (like "New Segment")
            push!(unique_values, marker.value)
            valid_markers += 1
        else
            @warn "Marker sample $sample_idx out of bounds (1:$n_samples), skipping"
        end
    end

    @info "Found $valid_markers valid markers with $(length(unique_values)) unique values: $(collect(unique_values))"

    # Create mapping from original values to sequential integers (1, 2, 3, ...)
    value_to_trigger = Dict{String,Int}()
    for (i, value) in enumerate(sort(collect(unique_values)))
        value_to_trigger[value] = i
    end

    # Second pass: assign sequential trigger values and original marker strings
    for marker in markers
        # Convert to 1-based if needed
        sample_idx = is_zero_based ? marker.sample + 1 : marker.sample

        if 1 <= sample_idx <= n_samples
            # Include all markers, even those with empty values
            triggers[sample_idx] = value_to_trigger[marker.value]
            triggers_info[sample_idx] = marker.value
        end
    end

    non_zero_triggers = count(x -> x != 0, triggers)
    non_empty_strings = count(x -> x != "", triggers_info)
    @info "Created triggers with $non_zero_triggers non-zero values and $non_empty_strings non-empty marker strings"

    return triggers, triggers_info
end
