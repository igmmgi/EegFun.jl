# Abstract types
abstract type EegData end
abstract type SingleDataFrameEeg <: EegData end
abstract type MultiDataFrameEeg <: EegData end

"""
    AnalysisInfo

Basic information about data preprocessing.

# Fields
- `reference::Symbol`: Reference type used (e.g., :avg, :mastoid, :none)
- `hp_filter::Float64`: High-pass filter cutoff in Hz (0.0 if none)
- `lp_filter::Float64`: Low-pass filter cutoff in Hz (0.0 if none)
"""
@kwdef mutable struct AnalysisInfo
    reference::Symbol = :none
    hp_filter::Float64 = 0.0
    lp_filter::Float64 = 0.0
end

# Concrete types with basic fields
mutable struct ContinuousData <: SingleDataFrameEeg
    data::DataFrame
    layout::DataFrame
    sample_rate::Int64
    analysis_info::AnalysisInfo
end
 
mutable struct ErpData <: SingleDataFrameEeg
    data::DataFrame
    layout::DataFrame
    sample_rate::Int64
    analysis_info::AnalysisInfo
    n_epochs::Int64
end
 
mutable struct EpochData <: MultiDataFrameEeg
    data::Vector{DataFrame}
    layout::DataFrame
    sample_rate::Int64
    analysis_info::AnalysisInfo
end


struct IntervalIdx
    interval_start::Int
    interval_end::Int
end

struct IntervalTime
    interval_start::Float64
    interval_end::Float64
end

mutable struct CoordXY
    x::Any
    y::Any
end

mutable struct Coord
    coord::Array{CoordXY}
end

# Basic information functions right with the types
channels(dat::EegData) = dat.layout.label
extra_channels(dat::EegData) = setdiff(propertynames(data(dat)), [channels(dat); :time; :sample; :triggers])
times(dat::SingleDataFrameEeg) = dat.data.time
times(dat::MultiDataFrameEeg) = first(dat.data).time  # assume all epochs are the same
sample_rate(dat::ContinuousData) = dat.sample_rate
sample_rate(dat::EegData) = dat.sample_rate
sample_rate(dat::DataFrame) = Int(1 / mean(diff(dat.time)))
reference(dat::EegData) = dat.analysis_info.reference
reference(dat::AnalysisInfo) = dat.reference
filter_info(dat::AnalysisInfo) = [dat.hp_filter, dat.lp_filter]
data(dat::SingleDataFrameEeg) = dat.data # single data frame
data(dat::MultiDataFrameEeg) = to_data_frame(dat) # single data frame with all epochs

# Timepoints, channel, epoch, and duration 
n_samples(dat::SingleDataFrameEeg) = nrow(dat.data)
n_samples(dat::MultiDataFrameEeg) = nrow(first(dat.data))
n_channels(dat::EegData) = length(channels(dat))
n_epochs(dat::SingleDataFrameEeg) = 1
n_epochs(dat::MultiDataFrameEeg) = length(dat.data)
duration(dat::EegData) = last(times(dat)) - first(times(dat))

n_average(dat::ErpData) = dat.n_epochs

# channel information
has_channels(dat::EegData, chans::Vector{Symbol}) = all(in(channels(dat)), chans)
common_channels(dat1::EegData, dat2::EegData) = intersect(channels(dat1), channels(dat2))


function Base.show(io::IO, dat::EegData)
    println(io, "Type: $(typeof(dat))")
    println(io, "Size: $(n_epochs(dat)) (epoch) x $(n_samples(dat)) (rows) x $(n_channels(dat)) (columns)")
    println(io, "Labels: ", print_vector(channels(dat)))
    println(io, "Duration: ", duration(dat), " S")
    println(io, "Sample Rate: ", sample_rate(dat))
end

function Base.show(io::IO, dat::AnalysisInfo)
    print(io, "Reference: ", reference(dat), ", ")
    print(io, "HP Filter: $(filter_info(dat)[1]), LP Filter: $(filter_info(dat)[2])")
end

########### ICA ##############
mutable struct IcaPrms
    l_rate::Float64
    max_iter::Int
    w_change::Float64
    anneal_deg::Float64
    anneal_step::Float64
    blowup::Float64
    blowup_fac::Float64
    max_weight::Float64
    restart_factor::Float64
    degconst::Float64
    default_stop::Float64
end

struct InfoIca
    unmixing::Matrix{Float64}
    mixing::Matrix{Float64}
    sphere::Matrix{Float64}
    variance::Vector{Float64}
    scale::Float64
    mean::Vector{Float64}
    ica_label::Vector{Symbol}
    data_label::Vector{Symbol}
end

"""
    Neighbours

Stores neighbor information for an electrode in layout-based operations.

# Fields
- `electrodes::Vector{Symbol}`: List of neighboring electrode labels
- `distances::Vector{Float64}`: Distances to each neighbor
- `weights::Vector{Float64}`: Interpolation weights for each neighbor
"""
struct Neighbours
    electrodes::Vector{Symbol}
    distances::Vector{Float64}
    weights::Vector{Float64}
end

