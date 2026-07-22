#  === EEG DATA UTILITIES ===
#
# This file provides utilities for accessing EegFun data structures.
# It includes functions for column identification, data access, convenience functions,
# data viewing, mathematical utilities, and subsetting operations.
#
# Organization:
# - Column Identification: Functions to identify metadata, channel, and extra columns
# - Data Access: Functions to extract different types of data (all_data, meta_data, etc.)
# - Convenience Functions: Basic information and size functions
# - Data Viewing: head, tail, viewer functions
# - Mathematical Utilities: _datarange, _colmeans, data_limits
# - Subsetting: Functions to subset EEG data objects

# === COLUMN IDENTIFICATION SYSTEM ===
"""
    _get_cols_by_group(dat::EegData, group::Symbol) -> Vector{Symbol}

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
function _get_cols_by_group(dat::EegData, group::Symbol)

    if !(group in [:channels, :metadata, :extra])
        @minimal_error "Unknown group type: $group"
    end


    labels = all_labels(dat)
    layout_channels = dat.layout.data.label

    if group == :channels
        return intersect(layout_channels, labels)
    elseif group == :metadata
        isempty(layout_channels) && return Symbol[]
        valid_channels = intersect(layout_channels, labels)
        isempty(valid_channels) && return Symbol[]
        first_channel_idx = findfirst(col -> col == valid_channels[1], labels)
        return labels[1:(first_channel_idx-1)]
    elseif group == :extra
        isempty(layout_channels) && return Symbol[]
        valid_channels = intersect(layout_channels, labels)
        isempty(valid_channels) && return Symbol[]
        last_channel_idx = findlast(col -> col == valid_channels[end], labels)
        return labels[(last_channel_idx+1):end]
    end
end


"""
    _get_cols_by_group(dat::Vector{<:EegData}, group::Symbol) -> Vector{Symbol}

Get columns by group type for a vector of EegData objects.
Delegates to the first element (assumes all have the same structure).

# Arguments
- `dat::Vector{<:EegData}`: Vector of EEG data objects
- `group::Symbol`: The group type (:metadata, :channels, :extra)

# Returns
- `Vector{Symbol}`: Column names of the specified group
"""
function _get_cols_by_group(dat::Vector{<:EegData}, group::Symbol)
    isempty(dat) && return Symbol[]
    return _get_cols_by_group(first(dat), group)
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
all_data(dat::DataFrame)::DataFrame = dat
all_data(dat::Vector{DataFrame})::DataFrame = isempty(dat) ? DataFrame() : vcat(dat...)

function all_data(dat::Vector{<:SingleDataFrameEeg})::DataFrame
    isempty(dat) && return DataFrame()
    dfs = map(dat) do d
        df = copy(d.data)
        insertcols!(df, 1, :condition => condition_number(d), :condition_name => condition_name(d); makeunique = true)
        df
    end
    return vcat(dfs...)
end

all_data(dat::Layout)::DataFrame = dat.data
all_data(dat::BiosemiDataFormat.BiosemiData)::DataFrame = _create_eegfun_dataframe(dat)
all_data(dat::EuropeanDataFormat.EdfData)::DataFrame = _create_eegfun_dataframe(dat)
all_data(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame = _create_eegfun_dataframe(dat)
all_data(dat::FunctionalImageFormat.FifData)::DataFrame = _create_eegfun_dataframe(dat)
all_data(dat::ExtensibleDataFormat.XdfData)::DataFrame = _create_eegfun_dataframe(dat)
all_data(dat::ErpMeasurementsResult)::DataFrame = dat.data
all_data(dat::TriggerInfo)::DataFrame = dat.data

function all_data(ica::InfoIca)::DataFrame
    removed = Set(keys(ica.removed_activations))
    return DataFrame(
        component = ica.ica_label,
        variance_pct = round.(ica.variance .* 100; digits = 2),
        is_sub_gaussian = ica.is_sub_gaussian,
        removed = [i in removed for i in eachindex(ica.ica_label)],
    )
end

# Handle collections and epoch selections
function all_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    return to_data_frame(subset(dat, epoch_selection = epoch_selection))
end

function meta_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    meta_cols = _get_cols_by_group(dat, :metadata)
    return isempty(meta_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[!, meta_cols]
end

function channel_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    channel_cols = _get_cols_by_group(dat, :channels)
    return isempty(channel_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[!, channel_cols]
end

function extra_data(dat::Union{MultiDataFrameEeg,Vector{<:MultiDataFrameEeg}}; epoch_selection::Function = epochs())
    extra_cols = _get_cols_by_group(dat, :extra)
    return isempty(extra_cols) ? DataFrame() : all_data(dat, epoch_selection = epoch_selection)[!, extra_cols]
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
all_labels(dat::DataFrame)::Vector{Symbol} = propertynames(dat)

"""
    time_vector(dat::EegData) -> Vector{Float64}
    time_vector(dat::EegData, epoch::Int) -> Vector{Float64}
    time_vector(df::DataFrame) -> Vector{Float64}

Get the time column from the EEG data or DataFrame as a vector.

# Arguments
- `dat::EegData`: The EEG data object
- `df::DataFrame`: A DataFrame with a `:time` column

# Returns
- `Vector{Float64}`: Time vector in seconds

# Notes
- For MultiDataFrameEeg types (EpochData, etc.), returns time from first epoch
- All epochs are assumed to have identical time vectors
"""
time_vector(dat::SingleDataFrameEeg)::Vector{Float64} = hasproperty(dat.data, :time) ? dat.data[!, :time] : Float64[]
time_vector(dat::MultiDataFrameEeg)::Vector{Float64} = hasproperty(dat.data[1], :time) ? dat.data[1][!, :time] : Float64[]
time_vector(dat::MultiDataFrameEeg, epoch::Int)::Vector{Float64} =
    hasproperty(dat.data[epoch], :time) ? dat.data[epoch][!, :time] : Float64[]
time_vector(dat::Vector{T}) where {T<:EegData} = isempty(dat) ? Float64[] : time_vector(dat[1])
time_vector(df::DataFrame)::Vector{Float64} = hasproperty(df, :time) ? df[!, :time] : Float64[]


"""
    meta_labels(dat::EegData) -> Vector{Symbol}

Get metadata column names from the EEG data.

# Arguments
- `dat::EegData`: The EEG data object

# Returns
- `Vector{Symbol}`: Vector of metadata column names
"""
meta_labels(dat::EegData) = _get_cols_by_group(dat, :metadata)

"""Return a DataFrame of just the columns in `group` (`:metadata`, `:channels`, `:extra`)."""
function _get_cols_data(dat::SingleDataFrameEeg, group::Symbol)
    cols = _get_cols_by_group(dat, group)
    return isempty(cols) ? DataFrame() : dat.data[!, cols]
end

function _get_cols_data(dat::MultiDataFrameEeg, group::Symbol, epoch::Int)
    cols = _get_cols_by_group(dat, group)
    return isempty(cols) ? DataFrame() : dat.data[epoch][!, cols]
end

function _get_cols_data(dat::MultiDataFrameEeg, group::Symbol)
    cols = _get_cols_by_group(dat, group)
    return isempty(cols) ? DataFrame() : vcat([df[!, cols] for df in dat.data]...)
end

"""
    meta_data(eeg_data::EegData) -> DataFrame

Get meta data columns from the EEG data.
"""
meta_data(dat::SingleDataFrameEeg) = _get_cols_data(dat, :metadata)
meta_data(dat::MultiDataFrameEeg, epoch::Int) = _get_cols_data(dat, :metadata, epoch)
meta_data(dat::MultiDataFrameEeg) = _get_cols_data(dat, :metadata)


"""
    channel_labels(dat::EegData) -> Vector{Symbol}

Get EEG channel column names from the EEG data.
"""
channel_labels(dat::EegData)::Vector{Symbol} = _get_cols_by_group(dat, :channels)
channel_labels(dat::EegData, channel_numbers::Vector{<:UnitRange})::Vector{Symbol} = channel_labels(dat)[channel_numbers...]
channel_labels(dat::EegData, channel_numbers)::Vector{Symbol} = channel_labels(dat)[channel_numbers]
channel_labels(dat::EegData, channel_numbers::Int)::Vector{Symbol} = channel_labels(dat)[[channel_numbers]]

# Handle collections of EEG data
# TODO: is it possible that all do not have the same channel labels?
channel_labels(dat::Vector{<:EegData})::Vector{Symbol} = channel_labels(first(dat))


"""
    channel_data(eeg_data::EegData) -> DataFrame

Get EEG channel data columns from the EEG data.
"""
channel_data(dat::SingleDataFrameEeg)::DataFrame = _get_cols_data(dat, :channels)
channel_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame = _get_cols_data(dat, :channels, epoch)
channel_data(dat::MultiDataFrameEeg)::DataFrame = _get_cols_data(dat, :channels)

# Extract channel signals as vectors (returns vector of signals for processing)
channel_data(dat::MultiDataFrameEeg, channel::Symbol) = reduce(vcat, (trial[!, channel] for trial in dat.data))
channel_data(dat::SingleDataFrameEeg, channel::Symbol)::Vector{Float64} = dat.data[!, channel]

"""
    extra_labels(dat::EegData) -> Vector{Symbol}

Get extra/derived column names (EOG, flags, etc.) from the EEG data.
"""
extra_labels(dat::EegData) = _get_cols_by_group(dat, :extra)
extra_labels(dat::Vector{<:EegData})::Vector{Symbol} = extra_labels(first(dat))

"""
    extra_data(eeg_data::EegData) -> DataFrame

Get extra/derived columns (EOG, flags, etc.) from the EEG data.
"""
extra_data(dat::SingleDataFrameEeg)::DataFrame = _get_cols_data(dat, :extra)
extra_data(dat::MultiDataFrameEeg, epoch::Int)::DataFrame = _get_cols_data(dat, :extra, epoch)
extra_data(dat::MultiDataFrameEeg)::DataFrame = _get_cols_data(dat, :extra)



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
"""
n_epochs(dat::SingleDataFrameEeg)::Int = 1
n_epochs(dat::MultiDataFrameEeg)::Int = length(dat.data)
n_epochs(dat::ErpData)::Int = dat.n_epochs

"""
    condition_number(dat::EegFunData) -> Union{Int, String}

Get the condition number or identifier.
"""
condition_number(dat::EegFunData) = hasproperty(dat, :condition) ? dat.condition : "N/A"
condition_number(dat::ContinuousData)::String = "Raw Data"

"""
    condition_name(dat::EegFunData) -> String

Get the condition name.
"""
condition_name(dat::EegFunData) = hasproperty(dat, :condition_name) ? dat.condition_name : "N/A"
condition_name(dat::ContinuousData)::String = "Raw Data"

"""
    file_name(dat::EegFunData) -> String

Get the source filename.
"""
function file_name(dat::EegFunData)
    hasproperty(dat, :file) && return dat.file
    hasproperty(dat, :filename) && return dat.filename
    return "Unknown"
end


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
    hasproperty(dat.data[1], :time) && !isempty(dat.data[1].time) ? last(dat.data[1].time) - first(dat.data[1].time) : 0.0
duration(dat::MultiDataFrameEeg, epoch::Int)::Float64 =
    hasproperty(dat.data[epoch], :time) && !isempty(dat.data[epoch].time) ? last(dat.data[epoch].time) - first(dat.data[epoch].time) : 0.0
duration(dat::Vector{T}) where {T<:EegData} = isempty(dat) ? 0.0 : duration(dat[1])


"""
    _have_same_structure(dat1::EegData, dat2::EegData) -> Bool
    _have_same_structure(dats::Vector{<:EegData}) -> Bool

Check if EegData objects have the same structure (sample rate, number of samples, channel labels, and time vectors).

# Arguments
- `dat1::EegData`, `dat2::EegData`: Two EegData objects to compare
- `dats::Vector{<:EegData}`: Vector of EegData objects to compare

# Returns
- `Bool`: `true` if all objects have the same sample rate, number of samples, channel labels, and time vectors

# Examples
```julia
# Check if two ERPs have the same structure
_have_same_structure(erp1, erp2)

# Check if all ERPs in a vector have the same structure
_have_same_structure(erps)
```
"""
function _have_same_structure(dat1::EegData, dat2::EegData)::Bool
    if sample_rate(dat1) != sample_rate(dat2)
        @minimal_error("Sample rates do not match: $(dat1.file)/$(dat1.condition) vs. $(dat2.file)/$(dat2.condition)")
        return false
    end
    if n_samples(dat1) != n_samples(dat2)
        @minimal_error("Number of samples do not match: $(dat1.file)/$(dat1.condition) vs. $(dat2.file)/$(dat2.condition)")
        return false
    end
    if channel_labels(dat1) != channel_labels(dat2)
        @minimal_error("Channel labels do not match: $(dat1.file)/$(dat1.condition) vs. $(dat2.file)/$(dat2.condition)")
        return false
    end
    # Check time vectors match (if they exist)
    time1 = time_vector(dat1)
    time2 = time_vector(dat2)
    if !isempty(time1) && !isempty(time2) && !all(time1 .≈ time2)
        @minimal_error("Time vectors do not match: $(dat1.file)/$(dat1.condition) vs. $(dat2.file)/$(dat2.condition)")
        return false
    end
    return true
end

function _have_same_structure(dats::Vector{<:EegData})::Bool
    isempty(dats) && return true
    length(dats) == 1 && return true

    # Check all against the first one
    template = first(dats)
    for i = 2:length(dats)
        !_have_same_structure(template, dats[i]) && return false
    end
    return true
end


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




"""
    to_data_frame(dat::MultiDataFrameEeg) -> DataFrame

Convert a multi-dataframe EEG object to a single DataFrame by concatenating all epochs.
"""
function to_data_frame(dat::MultiDataFrameEeg)
    isempty(dat.data) && return DataFrame()
    return vcat(dat.data...)
end
function to_data_frame(dat::Vector{EpochData})
    isempty(dat) && return DataFrame()
    isempty(dat[1].data) && return DataFrame()
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end

# === MATHEMATICAL UTILITIES ===
"""
    find_closest_time_index(times::AbstractVector{<:Real}, target::Real) -> Int

Find the index in a strictly sorted time vector that is closest to `target`.
"""
function find_closest_time_index(times::AbstractVector{<:Real}, target::Real)
    idx = searchsortedfirst(times, target)
    if idx == 1
        return 1
    elseif idx > length(times)
        return length(times)
    else
        return abs(times[idx] - target) < abs(times[idx-1] - target) ? idx : idx - 1
    end
end

"""
    _datarange(x::AbstractVector) -> Float64

Calculate the range of data (maximum - minimum).

# Returns
- `Float64`: Difference between maximum and minimum values
"""
_datarange(x::AbstractVector)::Float64 = -(-(extrema(x)...))


"""
    _colmeans(df::DataFrame, cols) -> Vector{Float64}
    _colmeans(df::Matrix) -> Vector{Float64}
    _colmeans(df::Matrix, cols) -> Vector{Float64}

Calculate the mean of specified columns in a DataFrame.

# Arguments
- `df::DataFrame`: The DataFrame containing the data.
- `cols`: The columns for which to calculate the mean. This can be a vector of column names or indices.

# Returns
- `Vector{Float64}`: A vector containing the mean of each specified column.
"""
function _colmeans(df::DataFrame, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = nrow(df)
    out = zeros(Float64, N)
    for c in cols
        out .+= df[!, c]
    end
    return out ./ length(cols)
end

function _colmeans(mat::Matrix)::Vector{Float64}
    N, C = size(mat)
    C == 0 && return Float64[]
    out = zeros(Float64, N)
    for c = 1:C
        @views out .+= mat[:, c]
    end
    return out ./ C
end

function _colmeans(mat::Matrix, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = size(mat, 1)
    out = zeros(Float64, N)
    for c in cols
        @views out .+= mat[:, c]
    end
    return out ./ length(cols)
end


"""
    _data_limits_x(dat::DataFrame) -> Union{Tuple{Float64,Float64}, Nothing}

Get the time range of the data.

# Returns
- `Tuple{Float64,Float64}`: Minimum and maximum time values
- `Nothing`: If the DataFrame is empty
"""
function _data_limits_x(dat::DataFrame; col::Symbol = :time)
    isempty(dat) && return nothing
    return extrema(dat[!, col])
end

"""
    _data_limits_y(dat::DataFrame, col) -> Union{Vector{Float64}, Nothing}

Get the value range for specified columns.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data
- `col::Symbol`: The column to get limits for

# Returns
- `Vector{Float64}`: [minimum, maximum] across specified column
- `Nothing`: If the DataFrame is empty
"""
function _data_limits_y(dat::DataFrame, col::Symbol)
    isempty(dat) && return nothing
    mn, mx = extrema(dat[!, col])
    return [mn, mx]
end

"""
    _data_limits_y(dat::DataFrame, cols::Vector{Symbol}) -> Union{Vector{Float64}, Nothing}

Get the value range across multiple specified columns.

# Arguments
- `dat::DataFrame`: The DataFrame containing the data
- `cols::Vector{Symbol}`: The columns to get limits for

# Returns
- `Vector{Float64}`: [minimum, maximum] across all specified columns
- `Nothing`: If the DataFrame is empty
"""
function _data_limits_y(dat::DataFrame, cols::Vector{Symbol})
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
    _subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> DataFrame

Create a subset of a DataFrame by selecting specific channels and samples.

# Arguments
- `df::DataFrame`: The DataFrame to subset
- `selected_channels::Vector{Symbol}`: Column names to include
- `selected_samples::Vector{Int}`: Row indices to include

# Returns
- `DataFrame`: Subset DataFrame with selected channels and samples
"""
function _subset_dataframe(df::DataFrame, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::DataFrame
    return df[selected_samples, selected_channels]
end


# === DEFAULT Y-RANGE HELPERS ===
"""
    yrange(dat::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1)

Compute a padded y-range for ERP data using channel selection predicate.
Returns (min, max) padded by `buffer` proportion.
"""
function _ylimits(
    dat::ErpData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}
    dat_sub = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
    chs = channel_labels(dat_sub)
    lims = _data_limits_y(dat_sub.data, chs)
    return (lims[1], lims[2])
end

"""
    yrange(dat::EpochData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1, average_only::Bool = false)

Compute a padded y-range for epoch data using channel selection predicate.
If `average_only=true`, uses the averaged waveform across epochs.
"""
function _ylimits(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}

    # Apply predicates via subset first
    dat_sub = subset(
        dat;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
    # Determine which value columns to use
    chs = channel_labels(dat_sub)
    # Compute limits per epoch and combine
    limits = map(df -> _data_limits_y(df, chs), dat_sub.data)
    limits = filter(!isnothing, limits)
    isempty(limits) && return (0.0, 1.0)
    min_val = minimum(lim[1] for lim in limits)
    max_val = maximum(lim[2] for lim in limits)
    return (min_val, max_val)
end




"""
    _subset_dataframes(dataframes::Vector{DataFrame}, selected_epochs::Vector{Int}, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> Vector{DataFrame}

Create subsets of multiple DataFrames by selecting specific epochs, channels, and samples.

# Arguments
- `dataframes::Vector{DataFrame}`: Vector of DataFrames to subset
- `selected_epochs::Vector{Int}`: Indices of epochs (DataFrames) to include
- `selected_channels::Vector{Symbol}`: Column names to include in each DataFrame
- `selected_samples::Vector{Int}`: Row indices to include in each DataFrame

# Returns
- `Vector{DataFrame}`: Vector of subset DataFrames
"""
function _subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame}
    return _subset_dataframe.(dataframes[selected_epochs], Ref(selected_channels), Ref(selected_samples))
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
function _subset_common(dat::MultiDataFrameEeg, epoch_selection, channel_selection, sample_selection, include_extra)
    @debug "Subsetting $(typeof(dat)): selecting epochs, channels and samples"
    # Get subset selected epochs, channels, samples, and layout
    selected_epochs = get_selected_epochs(dat, epoch_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    return selected_epochs, selected_channels, selected_samples, layout_subset
end

"""
    _interval_to_samples(interval::Interval)

Internal helper to convert a Interval (nothing, tuple, or range)
into a samples() predicate function for use in subset().

# Arguments
- `interval::Interval`: Time interval specification

# Returns
- `Function`: samples() predicate that selects the specified time range
"""
function _interval_to_samples(sel::TimeSelection)
    start_time, stop_time = sel.start, sel.stop
    return x -> begin
        t = x[!, :time]
        mask = (t .>= start_time) .& (t .<= stop_time)
        selected_times = t[mask]
        if !isempty(selected_times)
            actual_start = round(first(selected_times); digits = 4)
            actual_stop = round(last(selected_times); digits = 4)
            n = length(selected_times)
            @debug "times(): Requested $(start_time)s–$(stop_time)s, actual sample range $(actual_start)s–$(actual_stop)s ($(n) sample$(n > 1 ? "s" : ""))"
        elseif start_time == stop_time
            # Single time point with no exact match: snap to nearest sample
            _, idx = findmin(abs.(t .- start_time))
            @debug "times(): Requested $(start_time)s, snapped to nearest sample at $(round(t[idx]; digits=4))s"
            mask = fill(false, nrow(x))
            mask[idx] = true
        else
            @warn "times(): No samples found in range $(start_time)s–$(stop_time)s"
        end
        mask
    end
end

_interval_to_samples(::AllSelection) = samples()
_interval_to_samples(::Nothing) = samples()
_interval_to_samples(t::Tuple{Real, Real}) = _interval_to_samples(TimeSelection(t[1], t[2]))
_interval_to_samples(f::Function) = f

"""
    _combine_interval_sample(interval::Interval, sample_selection::Function)

Combine an interval selection with a sample selection predicate.
"""
function _combine_interval_sample(interval::Interval, sample_selection::Function)
    if interval isa Nothing || interval isa AllSelection
        return sample_selection
    else
        interval_sel = _interval_to_samples(interval)
        return x -> interval_sel(x) .& sample_selection(x)
    end
end

"""Construct a new data object of the same type with subsetted data and layout."""
_create_subset(dat::ContinuousData, ds, ls) = ContinuousData(dat.file, ds, ls, dat.sample_rate, copy(dat.analysis_info))
_create_subset(dat::ErpData, ds, ls) =
    ErpData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info), dat.n_epochs)
_create_subset(dat::EpochData, ds, ls) =
    EpochData(dat.file, dat.condition, dat.condition_name, ds, ls, dat.sample_rate, copy(dat.analysis_info))
_create_subset(dat::TimeFreqData, ds_power::DataFrame, ds_phase::DataFrame, ls) = TimeFreqData(
    dat.file,
    dat.condition,
    dat.condition_name,
    ds_power,
    ds_phase,
    ls,
    dat.sample_rate,
    dat.method,
    dat.baseline,
    copy(dat.analysis_info),
)
_create_subset(dat::TimeFreqEpochData, ds_power::Vector{DataFrame}, ds_phase::Vector{DataFrame}, ls) = TimeFreqEpochData(
    dat.file,
    dat.condition,
    dat.condition_name,
    ds_power,
    ds_phase,
    ls,
    dat.sample_rate,
    dat.method,
    dat.baseline,
    copy(dat.analysis_info),
)

# === SUBSET IMPLEMENTATIONS ===

"""
    subset(dat::SingleDataFrameEeg; channel_selection, sample_selection, interval_selection, include_extra)

Create a subset of single-DataFrame EEG data (ContinuousData, ErpData).

# Arguments
- `channel_selection::Function`: Channel predicate (default: `channels()` for all)
- `sample_selection::Function`: Sample predicate (default: `samples()` for all)
- `interval_selection::Interval`: Time interval (default: `times()` for all)
- `include_extra::Bool`: Include extra columns (default: `false`)

# Examples
```julia
# Select specific channels and time window
sub = subset(erp, channel_selection = channels(:Fz, :Cz), interval_selection = times(-0.2, 0.5))
```
"""
function subset(
    dat::T;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::T where {T<:SingleDataFrameEeg}
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, combined_sel, include_extra)
    dat_subset = _subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat, dat_subset, layout_subset)
end

function subset(
    dat::TimeFreqData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::TimeFreqData
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, combined_sel, include_extra)
    # Subset BOTH power and phase DataFrames
    power_subset = _subset_dataframe(dat.data_power, selected_channels, selected_samples)
    phase_subset = _subset_dataframe(dat.data_phase, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

"""
    subset(dat::MultiDataFrameEeg; channel_selection, sample_selection, interval_selection, epoch_selection, include_extra)

Create a subset of epoched EEG data (e.g., `EpochData`). 

# Arguments
- `epoch_selection::Function`: Epoch predicate. You can supply numerical indices using `epochs(1:5)` or a custom DataFrame predicate via `epochs(df -> nrow(df) > 100)` to efficiently filter trials by their contents (default: `epochs()` for all).
- `channel_selection::Function`: Channel predicate (default: `channels()` for all)
- `sample_selection::Function`: Sample predicate (default: `samples()` for all)
- `interval_selection::Interval`: Time interval (default: `times()` for all)
- `include_extra::Bool`: Include extra columns (default: `false`)

# Examples
```julia
# Drop trial outliers by examining the interval channel inside the epoch
cleaned_epochs = subset(epochs_data, epoch_selection = epochs(df -> 0.2 < df.interval[1] < 1.0))
```
"""
function subset(
    dat::T;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::T where {T<:MultiDataFrameEeg}
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    dat_subset = _subset_dataframes(dat.data, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, dat_subset, layout_subset)
end

"""
    subset(dat::TimeFreqEpochData; channel_selection, sample_selection, interval_selection, epoch_selection, include_extra)

Create a subset of time-frequency epoched data. Accepts the same arguments as epoched data, including functional `epoch_selection`.
"""
function subset(
    dat::TimeFreqEpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::TimeFreqEpochData
    # Combine interval and sample selection
    combined_sel = _combine_interval_sample(interval_selection, sample_selection)
    selected_epochs, selected_channels, selected_samples, layout_subset =
        _subset_common(dat, epoch_selection, channel_selection, combined_sel, include_extra)
    # Subset BOTH power and phase Vector{DataFrame}
    power_subset = _subset_dataframes(dat.data_power, selected_epochs, selected_channels, selected_samples)
    phase_subset = _subset_dataframes(dat.data_phase, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat, power_subset, phase_subset, layout_subset)
end

function subset(
    datasets::Vector{ErpData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Vector{ErpData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, and interval selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::Vector{EpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::Vector{EpochData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, interval, and epoch selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::Vector{TimeFreqData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Vector{TimeFreqData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, and interval selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        include_extra = include_extra,
    )
end

function subset(
    datasets::Vector{TimeFreqEpochData};
    condition_selection::Function = conditions(),
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    interval_selection::Interval = times(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::Vector{TimeFreqEpochData}
    # First filter by condition_selection
    selected_conditions = get_selected_conditions(datasets, condition_selection)
    datasets_filtered = datasets[selected_conditions]

    # Then apply channel, sample, interval, and epoch selection to each dataset
    return subset.(
        datasets_filtered;
        channel_selection = channel_selection,
        sample_selection = sample_selection,
        interval_selection = interval_selection,
        epoch_selection = epoch_selection,
        include_extra = include_extra,
    )
end



# === BATCH PROCESSING ===

function _default_subset_output_dir(input_dir::String, pattern::String)
    joinpath(input_dir, "subset_$(pattern)")
end

function _process_subset_file(filepath::String, output_path::String; kwargs...)
    filename = basename(filepath)
    dat = read_data(filepath)

    if isnothing(dat)
        return BatchResult(false, filename, "No data found in file")
    end

    try
        subset_data = subset(dat; kwargs...)
        jldsave(output_path; data = subset_data)

        # Determine stats for logging
        if subset_data isa Vector && !isempty(subset_data) && hasproperty(subset_data[1], :data)
            n_items = sum(length(cond.data) for cond in subset_data)
            return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
        elseif hasproperty(subset_data, :data) && subset_data.data isa Vector
            n_items = length(subset_data.data)
            return BatchResult(true, filename, "Processed subset. Elements remaining: $(n_items)")
        else
            return BatchResult(true, filename, "Processed subset.")
        end
    catch e
        return BatchResult(false, filename, "Error: $(sprint(showerror, e))")
    end
end

"""
    subset(file_pattern::String; kwargs...)

Batch process JLD2 files that match `file_pattern`, applying the `subset` operation.
Supports all `subset` kwargs (`epoch_selection`, `channel_selection`, `sample_selection`, `interval_selection`).
"""
function subset(
    file_pattern::String;
    input_dir::String = pwd(),
    participant_selection::Function = participants(),
    output_dir::Union{String,Nothing} = nothing,
    kwargs...,
)
    log_file = "subset.log"
    setup_global_logging(log_file)

    try
        @info "Batch subset started at $(now())"

        error_msg = _validate_input_dir(input_dir)
        if !isnothing(error_msg)
            @minimal_error(error_msg)
        end

        output_dir = something(output_dir, _default_subset_output_dir(input_dir, file_pattern))
        mkpath(output_dir)

        files = _find_batch_files(file_pattern, input_dir, participant_selection)

        if isempty(files)
            @minimal_warning "No JLD2 files found matching pattern '$(file_pattern)' in $(input_dir)"
            return nothing
        end

        @info "Found $(length(files)) JLD2 files matching pattern '$(file_pattern)'"

        process_fn = (input_path, output_path) -> _process_subset_file(input_path, output_path; kwargs...)

        results = _run_batch_operation(process_fn, files, input_dir, output_dir; operation_name = "Subsetting parameters")

        _log_batch_summary(results, output_dir)

    finally
        _cleanup_logging(log_file, output_dir)
    end
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
        # Default column_labels to column names only — this suppresses the auto-generated
        # type subheader row that PrettyTables adds when given a DataFrame directly.
        # Callers can still override by passing their own column_labels kwarg.
        kw = Dict{Symbol,Any}(kwargs)
        get!(kw, :column_labels, names(df))
        pretty_table(io_context, df; kw...)
    end

    # Log with specified level
    if log_level == :debug
        @debug "\n\n$table_output\n"
    elseif log_level == :info
        @info "\n\n$table_output\n"
    elseif log_level == :warn
        @minimal_warning "\n\n$table_output\n"
    elseif log_level == :error
        @error "\n\n$table_output\n"
    else
        @minimal_warning "Unknown log level: $log_level, using :info instead"
        @info "\n\n$table_output\n"
    end

    return nothing
end


"""
    channels()
    channels(name::Symbol)
    channels(names::Symbol...)
    channels(names::Vector{Symbol})
    channels(index::Int)
    channels(indices::Union{Vector{Int}, UnitRange})
    channels(ranges::Vector{UnitRange{Int}})
    channels(predicate::Function)

Predicate generators for channel selection in subsetting and plotting operations.

Returns a function that, given a vector of channel labels, produces a boolean mask
selecting the matching channels. Used with the `channel_selection` keyword argument.

# Examples
```julia
channels()                          # Select all channels (default)
channels(:Cz)                       # Select a single channel by name
channels(:Fz, :Cz, :Pz)            # Select multiple channels by name
channels([:Fz, :Cz, :Pz])          # Same, from a Vector
channels(1:10)                      # Select first 10 channels by index
channels(x -> startswith.(string.(x), "F"))  # Custom predicate
```
"""
channels() = x -> fill(true, length(x))  # Default: select all channels given
channels(channel_names::Vector{Symbol}) = x -> x .∈ Ref(channel_names)
channels(channel_name::Symbol) = x -> x .== channel_name
channels(channel_names::Symbol...) = channels(collect(channel_names))
channels(channel_number::Int) = channels([channel_number])
channels(channel_numbers::Union{Vector{Int},UnitRange}) = x -> [i in channel_numbers for i = 1:length(x)]
channels(channel_ranges::Vector{UnitRange{Int}}) = x -> [i in union(channel_ranges...) for i = 1:length(x)]
channels(predicate::Function) = predicate  # Allow custom function predicates
"""
    channels_not(name::Symbol)
    channels_not(names::Symbol...)
    channels_not(names::Vector{Symbol})
    channels_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified channels.
Inverse of `channels(...)`.

# Examples
```julia
channels_not(:Fp1, :Fp2)           # Exclude two frontal channels
channels_not(1:4)                   # Exclude first 4 channels by index
```
"""
channels_not(channel_names::Vector{Symbol}) = x -> .!(x .∈ Ref(channel_names))
channels_not(channel_name::Symbol) = x -> .!(x .== channel_name)
channels_not(channel_names::Symbol...) = channels_not(collect(channel_names))
channels_not(channel_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in channel_numbers for i = 1:length(x)])

"""
    components()
    components(index::Int)
    components(indices::Int...)
    components(indices::Union{Vector{Int}, UnitRange})
    components(predicate::Function)

Predicate generators for ICA component selection.

# Examples
```julia
components()                        # Select all components (default)
components(1, 2, 3)                 # Select components 1, 2, and 3
components(1:5)                     # Select first 5 components
```
"""
components() = x -> fill(true, length(x))  # Default: select all components given
components(component_numbers::Union{Vector{Int},UnitRange}) = x -> [i in component_numbers for i = 1:length(x)]
components(component_number::Int) = x -> x .== component_number
components(component_numbers::Int...) = components(collect(component_numbers))
components(predicate::Function) = predicate  # Allow custom function predicates
"""
    components_not(index::Int)
    components_not(indices::Int...)
    components_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified ICA components.

# Examples
```julia
components_not(1)                   # Exclude component 1
components_not(1:3)                 # Exclude first 3 components
```
"""
components_not(component_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in component_numbers for i = 1:length(x)])
components_not(component_number::Int) = x -> .!(x .== component_number)
components_not(component_numbers::Int...) = components_not(collect(component_numbers))

"""
    samples()
    samples(column::Symbol)
    samples(predicate::Function)

Predicate generators for sample (row) selection. Returns a function that,
given a DataFrame, produces a boolean mask over rows.

# Examples
```julia
samples()                           # Select all samples (default)
samples(:is_clean)                  # Select rows where :is_clean is true
samples(df -> df.time .> 0)         # Custom predicate on the DataFrame
```
"""
samples() = x -> fill(true, nrow(x))
samples(column::Symbol) = x -> x[!, column]
# Allow custom function predicates
samples(predicate::Function) = predicate
"""
    samples_or(columns::Vector{Symbol})
    samples_and(columns::Vector{Symbol})
    samples_not(column::Symbol)
    samples_or_not(columns::Vector{Symbol})
    samples_and_not(columns::Vector{Symbol})

Combine or negate sample predicates across boolean columns.

# Examples
```julia
samples_or([:is_extreme, :is_drift])    # Rows where ANY flag is true
samples_and([:is_clean, :is_valid])     # Rows where ALL flags are true
samples_not(:is_bad)                    # Rows where :is_bad is false
samples_or_not([:is_extreme, :is_bad])  # Rows where NO flag is true
```
"""
samples_or(columns::Vector{Symbol}) = x -> isempty(columns) ? fill(false, nrow(x)) : reduce((a,b)->a.|b, (x[!, col] for col in columns))
samples_and(columns::Vector{Symbol}) = x -> isempty(columns) ? fill(true, nrow(x)) : reduce((a,b)->a.&b, (x[!, col] for col in columns))
samples_not(column::Symbol) = x -> .!(x[!, column])
samples_or_not(columns::Vector{Symbol}) = x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a,b)->a.|b, (x[!, col] for col in columns)))
samples_and_not(columns::Vector{Symbol}) = x -> isempty(columns) ? fill(true, nrow(x)) : .!(reduce((a,b)->a.&b, (x[!, col] for col in columns)))

# Helper functions for time interval selection (Interval type)
# These return Interval values (tuples, ranges, or nothing), not predicates
"""
    times()
    times(t::Real)
    times(start::Real, stop::Real)
    times(interval::Tuple{Real,Real})

Create time interval values for the `interval_selection` keyword argument.
Unlike other predicate generators, `times()` returns an `AbstractSelection` value
(a `TimeSelection` or `AllSelection`), not a predicate function.

# Examples
```julia
times()                             # Select all time points (returns AllSelection)
times(0.3)                          # Select single time point (snaps to nearest sample)
times(-0.2, 0.5)                    # Select time window from -200ms to 500ms
times((-0.2, 0.5))                  # Same, from a Tuple
```
"""
times() = AllSelection()
times(t::Real) = TimeSelection(Float64(t), Float64(t))
times(start::Real, stop::Real) = TimeSelection(Float64(start), Float64(stop))
times(interval::Tuple{Real,Real}) = TimeSelection(Float64(interval[1]), Float64(interval[2]))




"""
    epochs()
    epochs(index::Int)
    epochs(indices::Int...)
    epochs(indices::Union{Vector{Int}, UnitRange})
    epochs(predicate::Function)

Predicate generators for epoch selection in subsetting operations.
The `predicate` variant receives individual epoch DataFrames and should return `Bool`.

# Examples
```julia
epochs()                              # Select all epochs (default)
epochs(1:5)                           # Select first 5 epochs by numerical index
epochs(2, 4, 6)                       # Select epochs 2, 4, and 6
epochs(df -> nrow(df) > 100)          # Keep epochs with >100 samples
epochs(df -> df.interval[1] < 1.0)    # Keep epochs by metadata value
```
"""
epochs() = x -> fill(true, length(x))  # Default: select all epochs given
epochs(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> [i in epoch_numbers for i = 1:length(x)]
epochs(epoch_number::Int) = x -> [i == epoch_number for i = 1:length(x)]
epochs(epoch_numbers::Int...) = epochs(collect(epoch_numbers))
epochs(predicate::Function) = x -> map(predicate, x)  # Allow custom function predicates over dataframes

"""
    epochs_not(index::Int)
    epochs_not(indices::Int...)
    epochs_not(indices::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified epochs. Operates inversely to `epochs`.

# Examples
```julia
epochs_not(1)                         # Exclude epoch 1
epochs_not(1:10)                      # Exclude first 10 epochs
```
"""
epochs_not(epoch_numbers::Union{Vector{Int},UnitRange}) = x -> .!([i in epoch_numbers for i = 1:length(x)])
epochs_not(epoch_number::Int) = x -> .!([i == epoch_number for i = 1:length(x)])
epochs_not(epoch_numbers::Int...) = epochs_not(collect(epoch_numbers))

"""
    participants()
    participants(id::Int)
    participants(ids::Int...)
    participants(ids::Union{Vector{Int}, UnitRange})
    participants(predicate::Function)

Predicate generators for participant ID selection in batch operations.

# Examples
```julia
participants()                        # Select all participants (default)
participants(1, 2, 3)                 # Select participants 1, 2, and 3
participants(1:10)                    # Select participants 1 through 10
```
"""
participants() = x -> fill(true, length(x))  # Default: select all participants given
participants(participant_ids::Union{Vector{Int},UnitRange}) = x -> [id in participant_ids for id in x]
participants(participant_id::Int) = x -> x .== participant_id
participants(participant_ids::Int...) = participants(collect(participant_ids))
participants(predicate::Function) = predicate  # Allow custom function predicates
"""
    participants_not(id::Int)
    participants_not(ids::Int...)
    participants_not(ids::Union{Vector{Int}, UnitRange})

Predicate generators that **exclude** the specified participant IDs.

# Examples
```julia
participants_not(5)                   # Exclude participant 5
participants_not(1:3)                 # Exclude participants 1, 2, and 3
```
"""
participants_not(participant_ids::Union{Vector{Int},UnitRange}) = x -> .!([id in participant_ids for id in x])
participants_not(participant_id::Int) = x -> .!(x .== participant_id)
participants_not(participant_ids::Int...) = participants_not(collect(participant_ids))

"""
    conditions()
    conditions(index::Int)
    conditions(indices::Int...)
    conditions(indices::Union{Vector{Int}, UnitRange})
    conditions(name::String)
    conditions(names::Vector{String})
    conditions(predicate::Function)

Predicate generators for condition selection by index or name.
Used with the `condition_selection` keyword argument on `Vector{ErpData}`,
`Vector{EpochData}`, and similar multi-condition containers.

# Examples
```julia
conditions()                          # Select all conditions (default)
conditions(1, 3)                      # Select conditions 1 and 3 by index
conditions("congruent")               # Select condition by name
conditions(["congruent", "neutral"])   # Select multiple conditions by name
```
"""
conditions() = x -> fill(true, length(x))  # Default: select all conditions given
conditions(condition_indices::Union{Vector{Int},UnitRange}) = x -> [i in condition_indices for i = 1:length(x)]
conditions(condition_index::Int) = x -> [i == condition_index for i = 1:length(x)]
conditions(condition_indices::Int...) = conditions(collect(condition_indices))
conditions(condition_names::Vector{String}) = x -> [condition_name(dat) in condition_names for dat in x]
conditions(cond_name::String) = x -> [condition_name(dat) == cond_name for dat in x]
conditions(predicate::Function) = predicate  # Allow custom function predicates
"""
    conditions_not(index::Int)
    conditions_not(indices::Int...)
    conditions_not(indices::Union{Vector{Int}, UnitRange})
    conditions_not(name::String)
    conditions_not(names::String...)
    conditions_not(names::Vector{String})

Predicate generators that **exclude** the specified conditions.

# Examples
```julia
conditions_not(1)                     # Exclude condition 1
conditions_not("incongruent")         # Exclude condition by name
```
"""
conditions_not() = x -> fill(true, length(x))
conditions_not(condition_indices::Union{Vector{Int},UnitRange}) = x -> .!([i in condition_indices for i = 1:length(x)])
conditions_not(condition_index::Int) = x -> .!([i == condition_index for i = 1:length(x)])
conditions_not(condition_indices::Int...) = conditions_not(collect(condition_indices))
conditions_not(condition_names::Vector{String}) = x -> .!([condition_name(dat) in condition_names for dat in x])
conditions_not(cond_name::String) = x -> .!([condition_name(dat) == cond_name for dat in x])
conditions_not(condition_names::String...) = conditions_not(collect(condition_names))

"""Preserve user-specified channel name ordering if all names are in the selection."""
function _handle_channel_names_order(user_order::Vector{Symbol}, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
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

    # Check if this is channels([...]) or channels_not([...])
    existing_in_order = [ch for ch in user_order if ch in selectable_cols]
    existing_in_selected = [ch for ch in existing_in_order if ch in selected]

    # If all existing channels are in selected, it's channels([...]) - preserve order
    if !isempty(existing_in_order) && length(existing_in_selected) == length(existing_in_order)
        result = Symbol[]
        for ch in user_order
            if ch ∉ selectable_cols || ch ∈ result
                continue
            end
            push!(result, ch)
        end
        return result
    end
    return nothing  # Use default order
end

"""Preserve user-specified channel number ordering if all indices are in the selection."""
function _handle_channel_numbers_order(user_order_numbers, selectable_cols::Vector{Symbol}, selected::Vector{Symbol})
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
    valid_indices = [i for i in user_order_numbers if 1 <= i <= length(selectable_cols)]
    selected_indices = [i for i in valid_indices if selectable_cols[i] in selected]

    # If all valid indices are in selected, it's channels([...]) - preserve order
    if !isempty(valid_indices) && length(selected_indices) == length(valid_indices)
        result = Symbol[]
        for i in user_order_numbers
            if i < 1 || i > length(selectable_cols) || selectable_cols[i] ∈ result
                continue
            end
            push!(result, selectable_cols[i])
        end
        return result
    end
    return nothing  # Use default order
end

"""Apply a channel predicate to data and return the matching column names (with optional metadata/extra)."""
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
        user_order_names = getfield(channel_selection, :channel_names)
        ordered = _handle_channel_names_order(user_order_names, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    elseif hasfield(selection_type, :channel_numbers)
        user_order_numbers = getfield(channel_selection, :channel_numbers)
        ordered = _handle_channel_numbers_order(user_order_numbers, selectable_cols, selected)
        !isnothing(ordered) && (selected = ordered)
    end

    return vcat(metadata_cols, selected)
end


"""Apply a component predicate and return matching indices."""
function get_selected_components(ica_result::InfoIca, component_selection::Function)
    all_components = 1:length(ica_result.ica_label)
    return all_components[component_selection(all_components)]
end

"""Apply a sample predicate and return matching row indices."""
function get_selected_samples(dat::SingleDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data))
end
function get_selected_samples(dat::SingleDataFrameEeg, ::Union{Nothing, AllSelection})
    return 1:nrow(dat.data)
end

# Helper to select samples based on a predicate
function get_selected_samples(dat::MultiDataFrameEeg, sample_selection::Function)
    return findall(sample_selection(dat.data[1])) # assume all data have the same samples
end
function get_selected_samples(dat::MultiDataFrameEeg, ::Union{Nothing, AllSelection})
    return 1:nrow(dat.data[1])
end

# Helper to select samples from a DataFrame
function get_selected_samples(dat::DataFrame, sample_selection::Function)
    return findall(sample_selection(dat))
end
function get_selected_samples(dat::DataFrame, ::Union{Nothing, AllSelection})
    return 1:nrow(dat)
end

# Helper to convert TimeSelection and Tuple into sample functions
function get_selected_samples(dat::Union{SingleDataFrameEeg, MultiDataFrameEeg, DataFrame}, sel::Union{TimeSelection, Tuple{Real, Real}})
    func = _interval_to_samples(sel)
    return get_selected_samples(dat, func)
end

"""Apply an epoch predicate and return matching epoch indices."""
function get_selected_epochs(dat::MultiDataFrameEeg, epoch_selection::Function)
    # The dat.data accessor maps to data_power for TimeFreqEpochData
    return findall(epoch_selection(dat.data))
end

"""Apply a condition predicate and return matching dataset indices."""
function get_selected_conditions(datasets::Vector{<:EegFunData}, condition_selection::Function)
    return findall(condition_selection(datasets))
end




"""
    epoch_to_continuous(dat::MultiDataFrameEeg, epoch_idx::Int) -> ContinuousData

Extract a single epoch from multi-epoch EEG data and return it as ContinuousData.

# Arguments
- `dat::MultiDataFrameEeg`: The multi-DataFrame EEG data
- `epoch_idx::Int`: Index of the epoch to extract (1-based)

# Returns
- `ContinuousData`: Single DataFrame containing only the specified epoch

# Examples
```julia
# Extract epoch 3 as continuous data
single_dat = epoch_to_continuous(dat, 3)

# Now you can use single DataFrame functions
summary = channel_summary(single_dat)
```
"""
function epoch_to_continuous(dat::T, epoch_idx::Int)::ContinuousData where {T<:MultiDataFrameEeg}
    # Validate epoch index
    if epoch_idx < 1 || epoch_idx > length(dat.data)
        @minimal_error "Epoch index $epoch_idx out of range (1:$(length(dat.data)))"
    end
    return ContinuousData(dat.file, dat.data[epoch_idx], dat.layout, dat.sample_rate, dat.analysis_info)
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

"""Rename DataFrame columns using a two-phase swap-safe approach."""
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
        duplicate_names = filter(x -> count(==(x), final_names) > 1, unique(final_names))
        @minimal_error "Cannot rename channels to duplicate names: $(join(duplicate_names, ", "))"
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
        tmp = only(last.(filter(p -> first(p) == old_name, temp_renames)))
        rename!(df, tmp => new_name)
    end
end

function _rename_data_columns!(dat::SingleDataFrameEeg, rename_dict::Dict{Symbol,Symbol}, existing_channels::Set{Symbol})
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
    _create_layout_from_labels(labels::Vector{Symbol})::Layout

Creates a fake layout with channel labels and zero positions for quick visualization.
Useful when no proper electrode layout is available.
"""
function _create_layout_from_labels(labels::Vector{Symbol})::Layout
    n_channels = length(labels)
    df = DataFrame(label = labels, inc = zeros(n_channels), azi = zeros(n_channels))
    return Layout(df, nothing, nothing, nothing)
end

"""Create a time+sample+trigger+channel DataFrame from raw BiosemiData."""
function _create_eegfun_dataframe(dat::BiosemiDataFormat.BiosemiData)::DataFrame
    df = DataFrame(time = dat.time, sample = 1:length(dat.time), trigger = _clean_triggers(dat.triggers.raw))

    channels = Symbol.(dat.header.channel_labels[1:(end-1)])
    for (i, ch) in enumerate(channels)
        # Allocate Float64 vector directly into the DataFrame column
        df[!, ch] = Float64.(@view dat.data[:, i])
    end

    return df
end


function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame from Biosemi data"
    file_name = filename(dat)
    df = _create_eegfun_dataframe(dat)
    return ContinuousData(file_name, df, layout, dat.header.sample_rate[1], AnalysisInfo())
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BiosemiDataFormat.BiosemiData)
    channel_labels = Symbol.(dat.header.channel_labels[1:(end-1)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw EdfData."""
function _create_eegfun_dataframe(dat::EuropeanDataFormat.EdfData)::DataFrame

    n_samples = size(dat.data, 1)
    n_channels = size(dat.data, 2)
    sample_rate = dat.header.sample_rate[1]

    time = collect(0:(n_samples-1)) ./ sample_rate
    sample = 1:n_samples

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    if !isnothing(dat.triggers) && length(dat.triggers.idx) > 0
        for i = 1:length(dat.triggers.idx)
            sample_idx = dat.triggers.idx[i]
            if 1 <= sample_idx <= n_samples
                trigger[sample_idx] = dat.triggers.val[i]
                trigger_info[sample_idx] = string(dat.triggers.val[i])
            end
        end
    end

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(Float64.(dat.data), Symbol.(dat.header.channel_labels[1:size(dat.data, 2)])),
    )
    return df
end

"""
    download_eegfun_datasets()

Automatically download and cache built-in tutorial datasets using `DataDeps.jl`.
Returns the absolute file path to the directory containing the downloaded datasets.

# Examples
```julia
data_dir = download_eegfun_datasets()
file_path = joinpath(data_dir, "participant1.bdf")
dat = read_raw_data(file_path)
```
"""
function download_eegfun_datasets()
    return datadep"TutorialDataSets"
end

function create_eegfun_data(dat::EuropeanDataFormat.EdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.edf)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    return ContinuousData(file_name, df, layout, dat.header.sample_rate[1], AnalysisInfo())
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::EuropeanDataFormat.EdfData)
    channel_labels = Symbol.(dat.header.channel_labels[1:size(dat.data, 2)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw BrainVision data."""
function _create_eegfun_dataframe(dat::BrainVisionDataFormat.BrainVisionData)::DataFrame

    # Check if data is available
    if isnothing(dat.data)
        @minimal_error "BrainVision data is empty (data field is nothing)"
    end

    if isnothing(dat.header)
        @minimal_error "BrainVision header is empty (header field is nothing)"
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
    channel_labels = dat.header.label[1:n_channels]

    # Verify channel count matches
    if length(channel_labels) != n_channels
        @minimal_error "Number of channel labels ($(length(channel_labels))) does not match number of channels ($n_channels)"
    end

    # Create trigger and marker strings columns from markers
    if isnothing(dat.markers) || isempty(dat.markers)
        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)
    else
        trigger, trigger_info = _extract_triggers_from_markers(dat.markers, n_samples)
    end

    # Create the DataFrame 
    df = hcat(DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info), DataFrame(dat.data, channel_labels))

    return df
end

"""
    create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData

Creates a ContinuousData object from a BrainVisionDataFormat data structure and a layout.

# Arguments
- `dat::BrainVisionDataFormat.BrainVisionData`: The BrainVisionDataFormat data structure containing EEG data.
- `layout::Layout`: The layout object containing electrode information.

# Returns
- `ContinuousData`: ContinuousData object containing the EEG data and layout information.

# Examples
```julia
# Create ContinuousData from BrainVisionDataFormat data and layout
eeg_data = create_eegfun_data(brainvision_data, layout)

# Quick visualization without layout (for databrowser only)
eeg_data = create_eegfun_data(brainvision_data)
```
"""
function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.vhdr/*.vmrk/*.eeg)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    return ContinuousData(file_name, df, layout, dat.header.Fs, AnalysisInfo())
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::BrainVisionDataFormat.BrainVisionData)
    channel_labels = Symbol.(dat.header.label[1:size(dat.data, 2)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
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
    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    @debug "Processing $(length(markers)) markers for $n_samples samples"

    # Debug: Show first few marker sample indices
    if length(markers) > 0
        first_few_samples = [marker.sample for marker in markers[1:min(5, length(markers))]]
        @debug "First few marker sample indices: $first_few_samples"
    end

    # First pass: extract all unique trigger values (including empty strings for system markers)
    # Check if markers are 0-based or 1-based by looking at the first marker
    is_zero_based = length(markers) > 0 && markers[1].sample == 0
    @debug "Detected $(is_zero_based ? "0-based" : "1-based") indexing for marker samples"

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
            @minimal_warning "Marker sample $sample_idx out of bounds (1:$n_samples), skipping"
        end
    end

    @debug "Found $valid_markers valid markers with $(length(unique_values)) unique values: $(collect(unique_values))"

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
            trigger[sample_idx] = value_to_trigger[marker.value]
            trigger_info[sample_idx] = marker.value
        end
    end

    non_zero_triggers = count(x -> x != 0, trigger)
    non_empty_strings = count(x -> x != "", trigger_info)
    @debug "Created trigger with $non_zero_triggers non-zero values and $non_empty_strings non-empty marker strings"

    return trigger, trigger_info
end

_is_stim_channel(ch) = ch.kind == 3 || startswith(uppercase(ch.ch_name), "STI") || uppercase(ch.ch_name) == "STATUS"

"""Create a time+sample+trigger+channel DataFrame from raw FifData."""
function _create_eegfun_dataframe(dat::FunctionalImageFormat.FifData)::DataFrame
    _, n_samples = size(dat.data)
    sample_rate = Float64(dat.header.sfreq)

    time = collect(0:(n_samples-1)) ./ sample_rate
    sample = 1:n_samples

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    if !isnothing(dat.triggers) && length(dat.triggers.idx) > 0
        for i = 1:length(dat.triggers.idx)
            sample_idx = dat.triggers.idx[i]
            if 1 <= sample_idx <= n_samples
                trigger[sample_idx] = dat.triggers.val[i]
                trigger_info[sample_idx] = string(dat.triggers.val[i])
            end
        end
    end

    # Filter out stimulus channels (kind == 3 or named Status/STI)
    data_ch_indices = Int[]
    scales = Float64[]
    channel_labels = Symbol[]

    for (i, ch) in enumerate(dat.header.ch_info)
        if !_is_stim_channel(ch)
            push!(data_ch_indices, i)
            push!(channel_labels, Symbol(ch.ch_name))
            if ch.unit == 107 # Volts -> Microvolts
                push!(scales, 1e6)
            else
                push!(scales, 1.0)
            end
        end
    end

    scaled_data = Float64.(dat.data[data_ch_indices, :]') .* scales'

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(scaled_data, channel_labels),
    )
    return df
end

function create_eegfun_data(dat::FunctionalImageFormat.FifData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.fif)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)
    return ContinuousData(file_name, df, layout, Int(round(dat.header.sfreq)), AnalysisInfo())
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifData)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end

"""Create a time+sample+trigger+channel DataFrame from raw XdfData."""
function _create_eegfun_dataframe(dat::ExtensibleDataFormat.XdfData)::DataFrame
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    eeg_stream = eeg_streams[1]

    n_samples = size(eeg_stream.time_series, 1)
    n_channels = eeg_stream.header.channel_count

    srate =
        eeg_stream.header.nominal_srate > 0 ? eeg_stream.header.nominal_srate :
        (eeg_stream.header.effective_srate > 0 ? eeg_stream.header.effective_srate : 1.0)

    # Enforce uniform time array for ContinuousData (required for DSP filters)
    time = collect((0:(n_samples-1)) ./ srate)
    sample = 1:n_samples

    channel_labels = String[]
    if haskey(eeg_stream.header.info, "desc") &&
       haskey(eeg_stream.header.info["desc"], "channels") &&
       haskey(eeg_stream.header.info["desc"]["channels"], "channel")
        channels_info = eeg_stream.header.info["desc"]["channels"]["channel"]
        if channels_info isa Vector
            for ch in channels_info
                push!(channel_labels, get(ch, "label", "Unknown"))
            end
        elseif channels_info isa Dict
            push!(channel_labels, get(channels_info, "label", "Unknown"))
        end
    end

    if length(channel_labels) != n_channels
        channel_labels = ["ch$i" for i = 1:n_channels]
    end

    trigger = zeros(Int, n_samples)
    trigger_info = fill("", n_samples)

    marker_streams = [s for s in values(dat.streams) if s.header.type == "Markers"]
    if !isempty(marker_streams)
        marker_stream = marker_streams[1]
        for (i, t) in enumerate(marker_stream.timestamps)
            idx = find_closest_time_index(eeg_stream.timestamps, t)
            marker_val = marker_stream.time_series[i, 1]
            int_val = tryparse(Int, string(marker_val))
            trigger[idx] = int_val !== nothing ? int_val : 1
            trigger_info[idx] = string(marker_val)
        end
    end

    df = hcat(
        DataFrame(time = time, sample = sample, trigger = trigger, trigger_info = trigger_info),
        DataFrame(Float64.(eeg_stream.time_series), Symbol.(channel_labels[1:n_channels])),
    )
    return df
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData, layout::Layout)::ContinuousData
    @info "Creating EEG DataFrame (*.xdf)"
    file_name = basename_without_ext(dat.filename)
    df = _create_eegfun_dataframe(dat)

    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    sample_rate = !isempty(eeg_streams) ? eeg_streams[1].header.nominal_srate : 1.0

    return ContinuousData(file_name, df, layout, Int(round(sample_rate)), AnalysisInfo())
end

function create_eegfun_data(dat::ExtensibleDataFormat.XdfData)
    eeg_streams = [s for s in values(dat.streams) if s.header.type == "EEG"]
    if isempty(eeg_streams)
        @minimal_error "No EEG stream found in XDF data"
    end
    eeg_stream = eeg_streams[1]
    n_channels = eeg_stream.header.channel_count

    channel_labels = String[]
    if haskey(eeg_stream.header.info, "desc") &&
       haskey(eeg_stream.header.info["desc"], "channels") &&
       haskey(eeg_stream.header.info["desc"]["channels"], "channel")
        channels_info = eeg_stream.header.info["desc"]["channels"]["channel"]
        if channels_info isa Vector
            for ch in channels_info
                push!(channel_labels, get(ch, "label", "Unknown"))
            end
        elseif channels_info isa Dict
            push!(channel_labels, get(channels_info, "label", "Unknown"))
        end
    end
    if length(channel_labels) != n_channels
        channel_labels = ["ch$i" for i = 1:n_channels]
    end

    layout = _create_layout_from_labels(Symbol.(channel_labels))
    return create_eegfun_data(dat, layout)
end

function _create_eegfun_epoch_dataframes(dat::FunctionalImageFormat.FifEpochs)::Vector{DataFrame}
    _, n_samples, n_epochs = size(dat.data)
    sample_rate = Float64(dat.header.sfreq)

    # Filter out stimulus channels
    data_ch_indices = Int[]
    scales = Float64[]
    channel_labels = Symbol[]

    stim_idx = nothing
    stim_info = nothing

    for (i, ch) in enumerate(dat.header.ch_info)
        if _is_stim_channel(ch) && stim_idx === nothing
            stim_idx = i
            stim_info = ch
        end
        if !_is_stim_channel(ch)
            push!(data_ch_indices, i)
            push!(channel_labels, Symbol(ch.ch_name))
            if ch.unit == 107 # Volts -> Microvolts
                push!(scales, 1e6)
            else
                push!(scales, 1.0)
            end
        end
    end

    dfs = Vector{DataFrame}(undef, n_epochs)

    for e = 1:n_epochs
        first_samp = 0
        if !isempty(dat.first_sample)
            first_samp = length(dat.first_sample) >= e ? dat.first_sample[e] : dat.first_sample[1]
        end
        time = collect(first_samp:(first_samp+n_samples-1)) ./ sample_rate
        sample = 1:n_samples
        epoch_col = fill(e, n_samples)

        trigger = zeros(Int, n_samples)
        trigger_info = fill("", n_samples)

        if stim_idx !== nothing
            raw_float = @views dat.data[stim_idx, :, e] ./ (stim_info.range * stim_info.cal)
            raw = Int32.(round.(raw_float))

            if raw[1] > 0
                trigger[1] = raw[1]
                trigger_info[1] = string(raw[1])
            end
            for i = 2:n_samples
                if raw[i] != raw[i-1] && raw[i] > 0
                    trigger[i] = raw[i]
                    trigger_info[i] = string(raw[i])
                end
            end
        end

        scaled_data = Float64.(dat.data[data_ch_indices, :, e]') .* scales'

        dfs[e] = hcat(
            DataFrame(time = time, sample = sample, epoch = epoch_col, trigger = trigger, trigger_info = trigger_info),
            DataFrame(scaled_data, channel_labels),
        )
    end

    return dfs
end

function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs, layout::Layout)::Union{EpochData,ErpData}
    file_name = basename_without_ext(dat.filename)
    dfs = _create_eegfun_epoch_dataframes(dat)

    condition_name = isempty(dat.comments) ? "FIF Epochs" : dat.comments[1]

    if length(dfs) == 1 && !isempty(dat.nave) && dat.nave[1] > 1
        @info "Creating EEG DataFrame (*.fif evoked)"
        return ErpData(file_name, 1, condition_name, dfs[1], layout, Int(round(dat.header.sfreq)), AnalysisInfo(), Int(dat.nave[1]))
    else
        @info "Creating EEG DataFrame (*.fif epochs)"
        return EpochData(file_name, 1, condition_name, dfs, layout, Int(round(dat.header.sfreq)), AnalysisInfo())
    end
end

# Optional Layout variant for quick databrowser visualization
function create_eegfun_data(dat::FunctionalImageFormat.FifEpochs)
    channel_labels = Symbol.([ch.ch_name for ch in dat.header.ch_info if !_is_stim_channel(ch)])
    layout = _create_layout_from_labels(channel_labels)
    return create_eegfun_data(dat, layout)
end
