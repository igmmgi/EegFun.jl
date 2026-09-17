# Tables.jl Interface for EegFun.jl types

import Tables

# 1. Declare that all EegData types are Tables
Tables.istable(::Type{<:EegData}) = true
Tables.rowaccess(::Type{<:EegData}) = true
Tables.rows(dat::EegData) = Tables.rows(Tables.columns(dat))
Tables.schema(dat::EegData) = Tables.schema(Tables.columns(dat))

# 1.5 DataAPI / DataFrames Base Interface Forwarding
# We forward standard DataFrame functions to the internal `.data` so users don't have to extract it.
import DataFrames: nrow, ncol, describe
import Base: names

nrow(dat::SingleDataFrameEeg) = nrow(dat.data)
ncol(dat::SingleDataFrameEeg) = ncol(dat.data)
names(dat::SingleDataFrameEeg) = names(dat.data)
describe(dat::SingleDataFrameEeg; kwargs...) = describe(dat.data; kwargs...)

nrow(dat::MultiDataFrameEeg) = isempty(dat.data) ? 0 : sum(nrow, dat.data)
ncol(dat::MultiDataFrameEeg) = isempty(dat.data) ? 0 : ncol(first(dat.data))
names(dat::MultiDataFrameEeg) = isempty(dat.data) ? String[] : names(first(dat.data))
describe(dat::MultiDataFrameEeg; kwargs...) = isempty(dat.data) ? DataFrame() : describe(vcat(dat.data...; source = nothing); kwargs...)

# 2. Expose the column-based interface (since we wrap DataFrames)
Tables.columnaccess(::Type{<:EegData}) = true

# 3. Inject metadata using Idiomatic Julia Multiple Dispatch
function _inject_metadata!(res::DataFrame, dat::ContinuousData, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    return res
end

function _inject_metadata!(res::DataFrame, dat::Union{ErpData,EpochData}, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    res.condition = fill(dat.condition, n_rows)
    res.condition_name = fill(dat.condition_name, n_rows)
    if dat isa ErpData
        res.n_epochs = fill(dat.n_epochs, n_rows)
    end
    return res
end

function _inject_metadata!(res::DataFrame, dat::Union{TimeFreqData,TimeFreqEpochData,SpectrumData}, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    res.condition = fill(dat.condition, n_rows)
    res.condition_name = fill(dat.condition_name, n_rows)
    res.method = fill(dat.method, n_rows)
    return res
end

# Fallback for any future EegData type that doesn't have specific metadata
function _inject_metadata!(res::DataFrame, dat::EegData, n_rows::Int)
    return res
end

# 4. Implement column access for SingleDataFrameEeg (ContinuousData, ErpData, etc.)
function Tables.columns(dat::SingleDataFrameEeg)
    df = dat.data

    # Create a shallow copy so we can insert metadata columns without modifying the original
    res = copy(df, copycols = false)

    # Inject standard metadata via dispatch
    _inject_metadata!(res, dat, nrow(df))

    return Tables.columns(res)
end

# 5. Implement column access for MultiDataFrameEeg (EpochData, etc.)
function Tables.columns(dat::MultiDataFrameEeg)
    # Concatenate all DataFrames in the vector (epoch column is already present internally)
    df = vcat(dat.data...; source = nothing)

    # Create a shallow copy for metadata injection
    res = copy(df, copycols = false)

    # Inject standard metadata via dispatch
    _inject_metadata!(res, dat, nrow(df))

    return Tables.columns(res)
end

# 6. Support vectors of EegData (e.g., from extract_epochs or average_epochs)
Tables.istable(::Type{<:AbstractVector{<:EegData}}) = true
Tables.columnaccess(::Type{<:AbstractVector{<:EegData}}) = true
Tables.rowaccess(::Type{<:AbstractVector{<:EegData}}) = true
Tables.rows(dats::AbstractVector{<:EegData}) = Tables.rows(Tables.columns(dats))
Tables.schema(dats::AbstractVector{<:EegData}) = Tables.schema(Tables.columns(dats))

function Tables.columns(dats::AbstractVector{<:EegData})
    # Convert each object to a DataFrame (which goes through our Tables interface above)
    # and then efficiently concatenate them all together.
    dfs = [DataFrame(d) for d in dats]
    return Tables.columns(vcat(dfs...; cols = :union))
end
