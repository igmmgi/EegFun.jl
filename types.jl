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

