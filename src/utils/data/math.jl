function find_closest_time_index(times::AbstractVector{<:Real}, target::Real)
    isempty(times) && @minimal_error("Cannot find closest time index in an empty time vector")
    idx = searchsortedfirst(times, target)
    idx == 1 && return 1
    idx > length(times) && return length(times)
    return abs(times[idx] - target) < abs(times[idx-1] - target) ? idx : idx - 1
end

function _datarange(x::AbstractVector)::Float64
    isempty(x) && return 0.0
    mn, mx = extrema(x)
    return Float64(mx - mn)
end

function _add_col!(out::Vector{Float64}, col::Vector{Float64}, N::Int)
    @inbounds @simd for i = 1:N
        out[i] += col[i]
    end
end

function _add_col!(out::Vector{Float64}, col::AbstractVector, N::Int)
    @inbounds @simd for i = 1:N
        out[i] += Float64(col[i])
    end
end

function _colmeans(df::DataFrame, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = nrow(df)
    out = zeros(Float64, N)
    for c in cols
        _add_col!(out, df[!, c], N)
    end
    out .*= (1.0 / length(cols))
    return out
end

function _colmeans(mat::Matrix)::Vector{Float64}
    N, C = size(mat)
    C == 0 && return Float64[]
    out = zeros(Float64, N)
    for c = 1:C
        @inbounds @simd for i = 1:N
            out[i] += mat[i, c]
        end
    end
    out .*= (1.0 / C)
    return out
end

function _colmeans(mat::Matrix, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = size(mat, 1)
    out = zeros(Float64, N)
    for c in cols
        @inbounds @simd for i = 1:N
            out[i] += mat[i, c]
        end
    end
    out .*= (1.0 / length(cols))
    return out
end

function _data_limits_x(dat::DataFrame; col::Symbol = :time)
    isempty(dat) && return nothing
    return extrema(dat[!, col])
end

function _nan_safe_extrema(v::AbstractVector{Float64})
    mn, mx = Inf, -Inf
    found = false
    @inbounds for x in v
        if !isnan(x)
            mn = min(mn, x)
            mx = max(mx, x)
            found = true
        end
    end
    return found ? (mn, mx) : nothing
end

function _nan_safe_extrema(v::AbstractVector)
    mn, mx = Inf, -Inf
    found = false
    @inbounds for x in v
        if !ismissing(x) && !(x isa AbstractFloat && isnan(x))
            val = Float64(x)
            mn = min(mn, val)
            mx = max(mx, val)
            found = true
        end
    end
    return found ? (mn, mx) : nothing
end

function _nan_safe_extrema_rows(col::Vector{Float64}, rows::Vector{Int})
    mn, mx = Inf, -Inf
    found = false
    @inbounds for idx in rows
        val = col[idx]
        if !isnan(val)
            mn = min(mn, val)
            mx = max(mx, val)
            found = true
        end
    end
    return found ? (mn, mx) : nothing
end

function _nan_safe_extrema_rows(col::AbstractVector, rows::Vector{Int})
    mn, mx = Inf, -Inf
    found = false
    @inbounds for idx in rows
        x = col[idx]
        if !ismissing(x) && !(x isa AbstractFloat && isnan(x))
            val = Float64(x)
            mn = min(mn, val)
            mx = max(mx, val)
            found = true
        end
    end
    return found ? (mn, mx) : nothing
end

function _data_limits_y(dat::DataFrame, col::Symbol)
    isempty(dat) && return nothing
    ext = _nan_safe_extrema(dat[!, col])
    return isnothing(ext) ? nothing : [ext[1], ext[2]]
end

function _data_limits_y(dat::DataFrame, cols::Vector{Symbol})
    (isempty(dat) || isempty(cols)) && return nothing
    global_min, global_max = Inf, -Inf
    found = false
    @inbounds for col in cols
        ext = _nan_safe_extrema(dat[!, col])
        if !isnothing(ext)
            global_min = min(global_min, ext[1])
            global_max = max(global_max, ext[2])
            found = true
        end
    end
    return found ? [global_min, global_max] : nothing
end

function _filter_ylimit_cols(raw_cols::Vector{Symbol}, lay_labels::AbstractVector{Symbol}, df0::DataFrame)::Vector{Symbol}
    return Symbol[c for c in raw_cols if (c in lay_labels) || (eltype(df0[!, c]) <: Real && !(eltype(df0[!, c]) <: Bool))]
end

function _selected_ylimit_cols(dat::EegData, df0::DataFrame, channel_selection, include_extra::Bool)::Vector{Symbol}
    raw_cols = get_selected_channels(dat, channel_selection; include_meta = false, include_extra = include_extra)
    !include_extra && return raw_cols
    return _filter_ylimit_cols(raw_cols, dat.layout.data[!, :label], df0)
end

function _ylimits(
    dat::Union{ErpData,ContinuousData};
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}
    df = dat.data
    isempty(df) && return (0.0, 1.0)
    cols = _selected_ylimit_cols(dat, df, channel_selection, include_extra)
    sel_samples = get_selected_samples(df, _combine_interval_sample(interval_selection, sample_selection))
    (isempty(cols) || isempty(sel_samples)) && return (0.0, 1.0)

    mn, mx = Inf, -Inf
    found = false
    for c in cols
        ext = _nan_safe_extrema_rows(df[!, c], sel_samples)
        if !isnothing(ext)
            mn = min(mn, ext[1])
            mx = max(mx, ext[2])
            found = true
        end
    end
    return found ? (mn, mx) : (0.0, 1.0)
end

function _ylimits(
    dat::EpochData;
    channel_selection = channels(),
    sample_selection = samples(),
    interval_selection::Interval = times(),
    include_extra::Bool = false,
)::Tuple{Float64,Float64}
    isempty(dat.data) && return (0.0, 1.0)
    df1 = dat.data[1]
    isempty(df1) && return (0.0, 1.0)
    cols = _selected_ylimit_cols(dat, df1, channel_selection, include_extra)
    sel_samples = get_selected_samples(df1, _combine_interval_sample(interval_selection, sample_selection))
    (isempty(cols) || isempty(sel_samples)) && return (0.0, 1.0)

    mn, mx = Inf, -Inf
    found = false
    for df in dat.data
        for c in cols
            ext = _nan_safe_extrema_rows(df[!, c], sel_samples)
            if !isnothing(ext)
                mn = min(mn, ext[1])
                mx = max(mx, ext[2])
                found = true
            end
        end
    end
    return found ? (mn, mx) : (0.0, 1.0)
end

function _ylimits(datasets::Vector{<:Union{ErpData,ContinuousData,EpochData}}; kwargs...)::Tuple{Float64,Float64}
    isempty(datasets) && return (0.0, 1.0)
    mn, mx = Inf, -Inf
    for d in datasets
        lo, hi = _ylimits(d; kwargs...)
        mn = min(mn, lo)
        mx = max(mx, hi)
    end
    return (mn, mx)
end
