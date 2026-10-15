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
# NOTE: Assumes all EegData objects in the collection have identical channel labels.
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
has_channels(dat::EegData, chans::Vector{Symbol})::Bool = all(ch -> ch in dat, chans)

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
    to_data_frame(dat::Vector{<:MultiDataFrameEeg}) -> DataFrame

Convert a multi-dataframe EEG object (or vector of conditions) to a single DataFrame
by concatenating all epochs. Works for `EpochData`, `TimeFreqEpochData`, and any
other `MultiDataFrameEeg` subtype.
"""
function to_data_frame(dat::MultiDataFrameEeg)
    isempty(dat.data) && return DataFrame()
    return vcat(dat.data...)
end
function to_data_frame(dat::Vector{<:MultiDataFrameEeg})
    isempty(dat) && return DataFrame()
    isempty(dat[1].data) && return DataFrame()
    return vcat([vcat(dat[idx].data[:]...) for idx in eachindex(dat)]...)
end
