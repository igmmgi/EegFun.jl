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
    mn, mx = extrema(x)
    return mx - mn
end


"""
    _colmeans(df::DataFrame, cols) -> Vector{Float64}
    _colmeans(df::Matrix) -> Vector{Float64}
    _colmeans(df::Matrix, cols) -> Vector{Float64}

Calculate the mean of specified columns in a DataFrame.

# Arguments
- `df::DataFrame`: The DataFrame containing the data.
- `cols`: The columns for which to calculate the mean. This can be a vector of column names or indices.

# Returns
- `Vector{Float64}`: A vector containing the mean of each specified column.
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
    isempty(dat) && return nothing
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
    yrange(dat::ErpData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1)

Compute a padded y-range for ERP data using channel selection predicate.
Returns (min, max) padded by `buffer` proportion.
"""
function _ylimits(
    dat::ErpData;
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
    return (lims[1], lims[2])
end

"""
    yrange(dat::EpochData; channel_selection::Function = channels(), sample_selection::Function = samples(), include_extra::Bool = false, buffer::Float64 = 0.1, average_only::Bool = false)

Compute a padded y-range for epoch data using channel selection predicate.
If `average_only=true`, uses the averaged waveform across epochs.
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




"""
    _subset_dataframes(dataframes::Vector{DataFrame}, selected_epochs::Vector{Int}, selected_channels::Vector{Symbol}, selected_samples::Vector{Int}) -> Vector{DataFrame}

Create subsets of multiple DataFrames by selecting specific epochs, channels, and samples.

# Arguments
- `dataframes::Vector{DataFrame}`: Vector of DataFrames to subset
- `selected_epochs::Vector{Int}`: Indices of epochs (DataFrames) to include
- `selected_channels::Vector{Symbol}`: Column names to include in each DataFrame
- `selected_samples::Vector{Int}`: Row indices to include in each DataFrame

# Returns
- `Vector{DataFrame}`: Vector of subset DataFrames
"""
function _subset_dataframes(
    dataframes::Vector{DataFrame},
    selected_epochs::Vector{Int},
    selected_channels::Vector{Symbol},
    selected_samples::Vector{Int},
)::Vector{DataFrame}
    return _subset_dataframe.(dataframes[selected_epochs], Ref(selected_channels), Ref(selected_samples))
end
