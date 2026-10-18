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
function _datarange(x::AbstractVector)::Float64
    isempty(x) && return 0.0
    mn, mx = extrema(x)
    return Float64(mx - mn)
end


"""
    _colmeans(df::DataFrame, cols) -> Vector{Float64}
    _colmeans(df::Matrix) -> Vector{Float64}
    _colmeans(df::Matrix, cols) -> Vector{Float64}

Calculate the row-wise mean across specified columns (channel-average waveform).

# Arguments
- `df::DataFrame`: The DataFrame containing the data.
- `cols`: The columns across which to average for each row. Can be a vector of column names or indices.

# Returns
- `Vector{Float64}`: A vector of length `nrow(df)` containing the mean across the specified columns at each sample.
"""
function _colmeans(df::DataFrame, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = nrow(df)
    n_cols = length(cols)
    out = zeros(Float64, N)
    for c in cols
        col_vec = df[!, c]
        if col_vec isa Vector{Float64}
            col_f64 = col_vec::Vector{Float64}
            @inbounds @simd for i = 1:N
                out[i] += col_f64[i]
            end
        else
            @inbounds @simd for i = 1:N
                out[i] += Float64(col_vec[i])
            end
        end
    end
    return out ./ n_cols
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
    return out ./ C
end

function _colmeans(mat::Matrix, cols)::Vector{Float64}
    isempty(cols) && return Float64[]
    N = size(mat, 1)
    n_cols = length(cols)
    out = zeros(Float64, N)
    for c in cols
        @inbounds @simd for i = 1:N
            out[i] += mat[i, c]
        end
    end
    return out ./ n_cols
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
    (isempty(dat) || isempty(cols)) && return nothing
    global_min, global_max = Inf, -Inf
    @inbounds for col in cols
        mn, mx = extrema(dat[!, col])
        global_min = min(global_min, mn)
        global_max = max(global_max, mx)
    end
    return [global_min, global_max]
end




# === DEFAULT Y-RANGE HELPERS ===
"""
    _ylimits(dat::Union{ErpData, ContinuousData}; channel_selection::Function = channels(), sample_selection::Function = samples(), interval_selection::Interval = times(), include_extra::Bool = false)

Compute y-limits (min, max) for single-dataframe EEG data (ERP or Continuous) using channel and time selection predicates.
"""
function _ylimits(
    dat::Union{ErpData, ContinuousData};
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
    isnothing(lims) && return (0.0, 1.0)
    return (lims[1], lims[2])
end

"""
    _ylimits(dat::EpochData; channel_selection::Function = channels(), sample_selection::Function = samples(), interval_selection::Interval = times(), include_extra::Bool = false)

Compute y-limits (min, max) for epoch data using channel and time selection predicates.
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
