# === EEG DATA UTILITIES ===
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
        isnothing(last_channel_idx) && @minimal_error "Last channel label not found in data"
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


condition_number(dat::EpochData)::Int = dat.data[1].condition[1]
condition_name(dat::EpochData)::String = dat.data[1].condition_name[1]
file_name(dat::EpochData)::String = dat.data[1].file[1]


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
    limits = filter(!isnothing, limits)
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
function subset_dataframes(dataframes::Vector{DataFrame}, selected_epochs::Vector{Int}, selected_channels::Vector{Symbol}, selected_samples::Vector{Int})::Vector{DataFrame}
    return subset_dataframe.(dataframes[selected_epochs], Ref(selected_channels), Ref(selected_samples))
end


"""
    epochs_table(epochs::Vector{EpochData}; io::Union{IO, Nothing} = stdout)

Display a pretty table showing epoch information and return the DataFrame.
"""
function epochs_table(epochs::Vector{EpochData}; io::Union{IO, Nothing} = stdout)::DataFrame
    isempty(epochs) && throw(ArgumentError("epochs vector cannot be empty"))
    
    df = _build_base_epochs_df(epochs)
    df.n_epochs = [n_epochs(epoch) for epoch in epochs]
    
    _print_epochs_table(df, io, ["Condition", "Condition Name", "N Epochs"], [:r, :l, :r])
    return df
end

"""
    epochs_table(epochs_original, epochs_cleaned; io::Union{IO, Nothing} = stdout)

Display comparison table between original and cleaned epochs and return DataFrame.
"""
function epochs_table(epochs_original::Vector{EpochData}, epochs_cleaned::Vector{EpochData}; io::Union{IO, Nothing} = stdout)::DataFrame
    length(epochs_original) != length(epochs_cleaned) && 
        throw(ArgumentError("epochs_original and epochs_cleaned must have same length"))
    
    df = _build_base_epochs_df(epochs_original)
    df.n_epochs_original = [n_epochs(epoch) for epoch in epochs_original]
    df.n_epochs_cleaned = [n_epochs(epoch) for epoch in epochs_cleaned]
    df.percentage = round.((df.n_epochs_cleaned ./ df.n_epochs_original) .* 100; digits=1)
    
    _print_epochs_table(df, io, [:l, :r, :l, :r, :r, :r])
    return df
end

# Helper functions to reduce repetition
function _build_base_epochs_df(epochs::Vector{EpochData})::DataFrame
    return DataFrame(
        file = [filename(epoch) for epoch in epochs],
        condition = [condition_number(epoch) for epoch in epochs],
        condition_name = [condition_name(epoch) for epoch in epochs]
    )
end

function _print_epochs_table(df::DataFrame, io::Union{IO, Nothing}, alignments::Vector{Symbol})
    if io !== nothing
        pretty_table(io, df; alignment = alignments, crop = :none, show_subheader = false)
    end
end

"""
    log_epochs_table(message::String, epochs...; kwargs...)

Log an epochs table with message and return the DataFrame.
Combines logging and table creation in one clean call.
"""
function log_epochs_table(message::String, epochs...; kwargs...)
    io_buffer = IOBuffer()
    df = epochs_table(epochs...; io = io_buffer, kwargs...)
    @info "$message\n$(String(take!(io_buffer)))"
    return df
end

"""
    log_pretty_table(message::String, df::DataFrame; kwargs...)

Log a pretty table with message. For general DataFrame logging.
Sets show_row_number=false and show_subheader=false by default for cleaner logs.
"""
function log_pretty_table(message::String, df::DataFrame; kwargs...)
    io_buffer = IOBuffer()
    pretty_table(io_buffer, df; show_row_number=false, show_subheader=false, kwargs...)
    @info "$message\n$(String(take!(io_buffer)))"
    return nothing
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
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    return selected_epochs, selected_channels, selected_samples, layout_subset
end

"""
    _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info) -> ContinuousData

Internal helper to create ContinuousData from subset DataFrame.
"""
function _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info)
    return ContinuousData(data_subset, layout, sample_rate, analysis_info)
end

"""
    _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info, n_epochs::Int) -> ErpData

Internal helper to create ErpData from subset DataFrame.
"""
function _create_subset(data_subset::DataFrame, layout, sample_rate::Int, analysis_info, n_epochs::Int)
    return ErpData(data_subset, layout, sample_rate, analysis_info, n_epochs)
end

"""
    _create_subset(data_subset::Vector{DataFrame}, layout, sample_rate::Int, analysis_info) -> EpochData

Internal helper to create EpochData from subset DataFrames.
"""
function _create_subset(data_subset::Vector{DataFrame}, layout, sample_rate::Int, analysis_info)
    return EpochData(data_subset, layout, sample_rate, analysis_info)
end

# === SUBSET IMPLEMENTATIONS ===

function subset(
    dat::SingleDataFrameEeg;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::SingleDataFrameEeg
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info)
end

function subset(
    dat::ErpData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    include_extra::Bool = false,
)::ErpData
    selected_channels, selected_samples, layout_subset = _subset_common(dat, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframe(dat.data, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info, dat.n_epochs)
end

function subset(
    dat::EpochData;
    channel_selection::Function = channels(),
    sample_selection::Function = samples(),
    epoch_selection::Function = epochs(),
    include_extra::Bool = false,
)::EpochData
    selected_epochs, selected_channels, selected_samples, layout_subset = _subset_common(dat, epoch_selection, channel_selection, sample_selection, include_extra)
    dat_subset = subset_dataframes(dat.data, selected_epochs, selected_channels, selected_samples)
    return _create_subset(dat_subset, layout_subset, dat.sample_rate, dat.analysis_info)
end
