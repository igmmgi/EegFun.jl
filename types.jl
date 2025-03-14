# # Abstract types
# abstract type EegData end
# abstract type SingleFrameEeg <: EegData end
# abstract type MultiFrameEeg <: EegData end
# 
# # Concrete types with basic fields
# mutable struct ContinuousData <: SingleFrameEeg
#     data::DataFrame
#     layout::Layout
#     history::Dict{String, Any}
# end
# 
# mutable struct ErpData <: SingleFrameEeg
#     data::DataFrame
#     layout::Layout
#     history::Dict{String, Any}
# end
# 
# mutable struct EpochData <: MultiFrameEeg
#     data::Vector{DataFrame}
#     layout::Layout
#     history::Dict{String, Any}
# end
# 
# # Outer constructors - provide convenient initialization with empty history
# ContinuousData(data, layout) = ContinuousData(data, layout, Dict{String,Any}())
# ErpData(data, layout) = ErpData(data, layout, Dict{String,Any}())
# EpochData(data, layout) = EpochData(data, layout, Dict{String,Any}())



# Define an abstract type for EEG Data
abstract type EegData end

# Extend existing types to inherit from EEGData
mutable struct ContinuousData <: EegData
    data::DataFrame
    layout::DataFrame
    sample_rate::Int
end

mutable struct ErpData <: EegData
    data::DataFrame
    layout::DataFrame
    sample_rate::Int
end

mutable struct EpochData <: EegData
    data::Vector{DataFrame}
    layout::DataFrame
    sample_rate::Int
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


function Base.show(io::IO, dat::ContinuousData)
    println(io, "Data: Rows ($(nrow(dat.data))) x Columns ($(ncol(dat.data))) ")
    println(io, "Labels: ", join(names(dat.data), ", "))
    println(io, "Sample Rate: ", dat.sample_rate)
end

function Base.show(io::IO, dat::EpochData)
    println(io, "Number of trials: ", length(dat.data))
    println(io, "Timepoints ($(nrow(dat.data[1]))) x Channels ($(ncol(dat.data[1]))) ")
    println(io, "Channel Labels: ", join(names(dat.data[1]), ", "))
    println(io, "Sample Rate: ", dat.sample_rate)
end

function Base.show(io::IO, dat::ErpData)
    println(io, "Data: Rows ($(nrow(dat.data))) x Columns ($(ncol(dat.data))) ")
    println(io, "Labels: ", join(names(dat.data), ", "))
    println(io, "Sample Rate: ", dat.sample_rate)
end

# Basic information functions right with the types
channels(dat::Union{ContinuousData,ErpData,EpochData}) = dat.layout.label
times(dat::Union{ContinuousData,ErpData}) = dat.times
times(dat::EpochData) = first(dat.data).times
sample_rate(dat::Union{ContinuousData,ErpData,EpochData}) = dat.sample_rate
data(dat::Union{ContinuousData,ErpData}) = dat.data
data(dat::EpochData) = dat.data


# Data structure interfaces
"""Number of time points"""
n_times(dat::SingleFrameEeg) = length(times(dat))
n_times(dat::MultiFrameEeg) = length(times(dat))

"""Number of channels"""
n_channels(dat::EegData) = length(channels(dat))

"""Number of epochs (1 for SingleFrameEeg)"""
n_epochs(dat::SingleFrameEeg) = 1
n_epochs(dat::MultiFrameEeg) = length(data(dat))

# Validation interfaces
"""Check if channels exist"""
has_channels(dat::EegData, chans::Vector{Symbol}) = all(in(channels(dat)), chans)

"""Get channels present in both datasets"""
common_channels(dat1::EegData, dat2::EegData) = intersect(channels(dat1), channels(dat2))

# Time interfaces
"""Duration in seconds"""
duration(dat::EegData) = last(times(dat)) - first(times(dat))

"""Time resolution in seconds"""
time_resolution(dat::EegData) = 1/sample_rate(dat)

# Layout interfaces
"""Get channel coordinates"""
channel_coords(dat::EegData) = (dat.layout.x, dat.layout.y)
channel_radius(dat::EegData) = dat.layout.radius

# Data access interfaces
"""Get data for specific channels"""
function channel_data(dat::SingleFrameEeg, chans::Vector{Symbol})
    validate_channels(dat, chans)
    return dat.data[!, chans]
end

"""Get data at specific time point"""
function time_point_data(dat::SingleFrameEeg, time::Real)
    idx = findnearest(times(dat), time)
    return dat.data[idx, channels(dat)]
end


# # Add history field to abstract type interface
# history(dat::EegData) = dat.history
# 
# # Helper function to add processing steps
# function add_history!(dat::EegData, operation::String)
#     timestamp = Dates.now()
#     push!(dat.history, "[$timestamp] $operation")
#     return nothing
# end
# 
# # Use in processing functions
# function rereference!(dat::EegData, reference::Symbol)
#     # Do rereferencing...
#     
#     # Add to history
#     if reference == :avg
#         add_history!(dat, "Rereferenced to average reference")
#     elseif reference == :mastoid
#         add_history!(dat, "Rereferenced to linked mastoids (M1, M2)")
#     else
#         add_history!(dat, "Rereferenced to channels: $(join(reference, ", "))")
#     end
# end
# 
# # Show history in detailed display
# function Base.show(io::IO, ::MIME"text/plain", dat::EegData)
#     show(io, dat)  # Show basic info
#     println(io, "\nProcessing history:")
#     for entry in history(dat)
#         println(io, "  ", entry)
#     end
# end



# function Base.show(io::IO, dat::EegData)
#     print(io, "$(typeof(dat)) with:")
#     print(io, "\n  $(n_channels(dat)) channels")
#     print(io, "\n  $(sample_rate(dat)) Hz sampling rate")
#     print(io, "\n  $(duration(dat)) seconds duration")
# end
# 
# # Additional info for SingleFrameEeg
# function Base.show(io::IO, dat::SingleFrameEeg)
#     show(io, dat::EegData)  # Show common info first
#     print(io, "\n  $(n_times(dat)) time points")
#     print(io, "\n  Reference: $(dat.reference)")
# end
# 
# # Additional info for MultiFrameEeg
# function Base.show(io::IO, dat::MultiFrameEeg)
#     show(io, dat::EegData)  # Show common info first
#     print(io, "\n  $(n_epochs(dat)) epochs")
#     print(io, "\n  $(n_times(dat)) time points per epoch")
#     print(io, "\n  Reference: $(dat.reference)")
# end










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
    ica_label::Vector{String}
    data_label::Vector{String}
end
data(dat::Union{ContinuousData,ErpData}) = dat.data
data(dat::EpochData) = dat.data



# First, define your types and basic functions
# abstract type EegData end
# # ... other type definitions ...
# 
# # Define a function to add history tracking to a method
# function add_history_tracking(f::Function)
#     # Get all methods of the function
#     for m in methods(f)
#         # Check if method operates on EegData
#         if any(T <: EegData for T in m.sig.parameters)
#             # Create wrapper with history tracking
#             @eval function $(f.name)(args...)
#                 timestamp = Dates.now()
#                 result = $(Base.invokelatest)($(f), args...)
#                 if result isa EegData
#                     call_str = string($(f.name), "(", join(typeof.(args), ", "), ")")
#                     add_history!(result, "[$timestamp] $call_str")
#                 end
#                 return result
#             end
#         end
#     end
# end
# 
# # At end of your package, add tracking to all relevant functions
# function __init__()
#     # Get all functions in current scope
#     for name in names(@__MODULE__, all=true)
#         if isdefined(@__MODULE__, name)
#             f = getfield(@__MODULE__, name)
#             # Check if it's a function that modifies EegData (ends with !)
#             if f isa Function && endswith(string(name), "!")
#                 add_history_tracking(f)
#             end
#         end
#     end
# end


# # Basic show method (without history)
# function Base.show(io::IO, dat::EegData)
#     print(io, "$(typeof(dat)) with:")
#     print(io, "\n  $(n_channels(dat)) channels")
#     print(io, "\n  $(sample_rate(dat)) Hz sampling rate")
#     print(io, "\n  $(duration(dat)) seconds duration")
# end
# 
# # Detailed show method
# function Base.show(io::IO, ::MIME"text/plain", dat::EegData)
#     show(io, dat)
# end
# 
# # Convenience function for showing history
# """
#     show_history(dat::EegData)
# 
# Display data information including processing history.
# """
# function show_history(dat::EegData)
#     show(stdout, dat)
#     if !isempty(history(dat))
#         println("\nProcessing history:")
#         for entry in history(dat)
#             println("  ", entry)
#         end
#     end
# end
# 
# # Usage:
# julia> data  # normal display
# ContinuousData with:
#   64 channels
#   512 Hz sampling rate
#   300.0 seconds duration
# 
# julia> show_history(data)  # with history
# ContinuousData with:
#   64 channels
#   512 Hz sampling rate
#   300.0 seconds duration
# Processing history:
#   [2024-01-20T10:30:15] rereference!(EegData, Symbol)
#   [2024-01-20T10:30:16] filter!(EegData, Tuple{Float64,Float64})