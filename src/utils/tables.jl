import Tables
import DataFrames: nrow, ncol, describe, eachcol, eachrow
import Base: names, first, last

Tables.istable(::Type{<:EegData}) = true
Tables.columnaccess(::Type{<:EegData}) = true
Tables.rowaccess(::Type{<:EegData}) = true
Tables.rows(dat::EegData) = Tables.rows(Tables.columns(dat))
Tables.schema(dat::EegData) = Tables.schema(Tables.columns(dat))

nrow(dat::SingleDataFrameEeg) = nrow(dat.data)
ncol(dat::SingleDataFrameEeg) = ncol(dat.data)
names(dat::SingleDataFrameEeg) = names(dat.data)
names(dat::SingleDataFrameEeg, cols) = names(dat.data, cols)
describe(dat::SingleDataFrameEeg; kwargs...) = describe(dat.data; kwargs...)
eachcol(dat::SingleDataFrameEeg) = eachcol(dat.data)
eachrow(dat::SingleDataFrameEeg) = eachrow(dat.data)
first(dat::SingleDataFrameEeg, n::Integer) = first(dat.data, n)
last(dat::SingleDataFrameEeg, n::Integer) = last(dat.data, n)
first(dat::SingleDataFrameEeg) = first(dat.data)
last(dat::SingleDataFrameEeg) = last(dat.data)

nrow(dat::MultiDataFrameEeg) = isempty(dat.data) ? 0 : sum(nrow, dat.data)
ncol(dat::MultiDataFrameEeg) = isempty(dat.data) ? 0 : ncol(first(dat.data))
names(dat::MultiDataFrameEeg) = isempty(dat.data) ? String[] : names(first(dat.data))
names(dat::MultiDataFrameEeg, cols) = isempty(dat.data) ? String[] : names(first(dat.data), cols)
describe(dat::MultiDataFrameEeg; kwargs...) =
    isempty(dat.data) ? DataFrame() : describe(reduce(vcat, dat.data; source = nothing); kwargs...)
first(dat::MultiDataFrameEeg, n::Integer) = first(dat.data, n)
last(dat::MultiDataFrameEeg, n::Integer) = last(dat.data, n)

function _inject_metadata!(res::DataFrame, dat::ContinuousData, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    return res
end

function _inject_metadata!(res::DataFrame, dat::Union{ErpData,EpochData}, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    res.condition = fill(dat.condition, n_rows)
    res.condition_name = fill(dat.condition_name, n_rows)
    dat isa ErpData && (res.n_epochs = fill(dat.n_epochs, n_rows))
    return res
end

function _inject_metadata!(res::DataFrame, dat::Union{TimeFreqData,TimeFreqEpochData,SpectrumData}, n_rows::Int)
    res.file = fill(dat.file, n_rows)
    res.condition = fill(dat.condition, n_rows)
    res.condition_name = fill(dat.condition_name, n_rows)
    res.method = fill(dat.method, n_rows)
    return res
end

_inject_metadata!(res::DataFrame, ::EegData, ::Int) = res

function Tables.columns(dat::SingleDataFrameEeg)
    df = dat.data
    res = copy(df; copycols = false)
    _inject_metadata!(res, dat, nrow(df))
    return Tables.columns(res)
end

function Tables.columns(dat::MultiDataFrameEeg)
    if isempty(dat.data)
        res = DataFrame()
        for ch in dat.layout.data.label
            res[!, ch] = Float64[]
        end
        _inject_metadata!(res, dat, 0)
        return Tables.CopiedColumns(Tables.columns(res))
    end
    df = reduce(vcat, dat.data; source = nothing)
    _inject_metadata!(df, dat, nrow(df))
    return Tables.CopiedColumns(Tables.columns(df))
end

Tables.istable(::Type{<:AbstractVector{<:EegData}}) = true
Tables.columnaccess(::Type{<:AbstractVector{<:EegData}}) = true
Tables.rowaccess(::Type{<:AbstractVector{<:EegData}}) = true
Tables.rows(dats::AbstractVector{<:EegData}) = Tables.rows(Tables.columns(dats))
Tables.schema(dats::AbstractVector{<:EegData}) = Tables.schema(Tables.columns(dats))

function Tables.columns(dats::AbstractVector{<:EegData})
    isempty(dats) && return Tables.CopiedColumns(Tables.columns(DataFrame()))
    dfs = [DataFrame(d) for d in dats]
    return Tables.CopiedColumns(Tables.columns(reduce(vcat, dfs; cols = :union)))
end
