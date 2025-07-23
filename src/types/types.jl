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

mutable struct Layout
    data::DataFrame
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

"""
    EpochCondition

Defines parameters for extracting epochs for a specific experimental condition.

# Fields
- `name::String`: Descriptive condition name
- `trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}`: Trigger sequences to match (e.g., [[1, 2, 3]], [[1, :any, 3]], [[1:5], [10:15]])
- `reference_index::Int`: Which trigger position is t=0 (1-based, default: 1)
- `timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}`: Which trigger pairs to apply min/max intervals to (optional, default: nothing)
- `min_interval::Union{Nothing,Float64}`: Minimum time between specified trigger pairs (optional, default: nothing)
- `max_interval::Union{Nothing,Float64}`: Maximum time between specified trigger pairs (optional, default: nothing)
- `after::Union{Nothing,Int}`: Only search for sequences after this trigger value (optional, default: nothing)
- `before::Union{Nothing,Int}`: Only search for sequences before this trigger value (optional, default: nothing)
"""
@kwdef struct EpochCondition
    name::String
    trigger_sequences::Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}
    reference_index::Int = 1
    timing_pairs::Union{Nothing,Vector{Tuple{Int,Int}}} = nothing
    min_interval::Union{Nothing,Float64} = nothing
    max_interval::Union{Nothing,Float64} = nothing
    after::Union{Nothing,Int} = nothing
    before::Union{Nothing,Int} = nothing
end




# Basic information functions right with the types
channels(dat::EegData) = dat.layout.label
all_channels(dat::EegData) = propertynames(dat.data)
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
samples(dat::SingleDataFrameEeg) = dat.data.sample
samples(dat::MultiDataFrameEeg, epoch::Int) = dat.data[epoch].sample
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
    removed_activations::OrderedDict{Int, Matrix{Float64}}
end

# Custom display method for InfoIca
function Base.show(io::IO, ica::InfoIca)
    n_components = length(ica.ica_label)
    n_channels = length(ica.data_label)
    
    println(io, "InfoIca Result")
    println(io, "├─ Components: $n_components")
    println(io, "├─ Channels: $n_channels")
    println(io, "├─ Scale: $(round(ica.scale, digits=3))")
    println(io, "├─ Top 5 variance explained:")
    
    # Show top 5 components with their variance
    for i in 1:min(5, n_components)
        var_pct = round(ica.variance[i] * 100, digits=2)
        println(io, "│  $(ica.ica_label[i]): $(var_pct)%")
    end
    
    if n_components > 5
        println(io, "│  ... and $(n_components - 5) more components")
    end
    
    println(io, "├─ Matrix sizes:")
    println(io, "│  ├─ Unmixing: $(size(ica.unmixing))")
    println(io, "│  ├─ Mixing: $(size(ica.mixing))")
    println(io, "│  └─ Sphere: $(size(ica.sphere))")
    
    println(io, "└─ Channel labels: $(join(ica.data_label[1:min(5, n_channels)], ", "))$(n_channels > 5 ? " ..." : "")")
end

# Compact display for arrays
function Base.show(io::IO, ::MIME"text/plain", ica::InfoIca)
    show(io, ica)
end

# Custom copy method for InfoIca
function Base.copy(ica::InfoIca)::InfoIca
    return InfoIca(
        copy(ica.unmixing),
        copy(ica.mixing),
        copy(ica.sphere),
        copy(ica.variance),
        ica.scale,
        copy(ica.mean),
        copy(ica.ica_label),
        copy(ica.data_label),
        copy(ica.removed_activations)
    )
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

