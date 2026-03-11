"""
    _tf_build_row_index(df) -> Dict{Tuple{Float64,Float64}, Int}

Build a lookup from (rounded_freq, rounded_time) → row index for a TF DataFrame.
Uses `round(x, digits=6)` for robust floating-point matching.
"""
function _tf_build_row_index(df::DataFrame)
    idx = Dict{Tuple{Float64,Float64},Int}()
    sizehint!(idx, nrow(df))
    for i = 1:nrow(df)
        key = (round(df.freq[i], digits = 6), round(df.time[i], digits = 6))
        if !haskey(idx, key)
            idx[key] = i
        end
    end
    return idx
end

"""
    _tf_df_to_matrix(df, channel, frequencies, time_points[, row_index]) -> Matrix{Float64}

Extract a `[n_freqs × n_time]` power matrix for one channel from a TF DataFrame.

Uses a Dict-based row index for O(1) lookups. Optionally accepts a precomputed `row_index`
(from `_tf_build_row_index`) to avoid rebuilding it when extracting multiple channels
from the same DataFrame.

Missing cells are filled with `NaN`.
"""
function _tf_df_to_matrix(
    df::DataFrame,
    channel::Symbol,
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df),
)
    n_freqs = length(frequencies)
    n_time = length(time_points)
    mat = fill(NaN, n_freqs, n_time)

    for (ti, t) in enumerate(time_points)
        rt = round(t, digits = 6)
        for (fi, f) in enumerate(frequencies)
            row = get(row_index, (round(f, digits = 6), rt), nothing)
            if !isnothing(row)
                mat[fi, ti] = df[row, channel]
            end
        end
    end

    return mat
end

"""
    _tf_matrix_to_df!(df, channel, mat, frequencies, time_points[, row_index])

Write a `[n_freqs × n_time]` matrix back into the DataFrame column for `channel`.

Uses the same Dict-based row index as `_tf_df_to_matrix` for consistency.
Only writes cells that exist in the DataFrame (skips missing rows).
"""
function _tf_matrix_to_df!(
    df::DataFrame,
    channel::Symbol,
    mat::Matrix{Float64},
    frequencies::Vector{Float64},
    time_points::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df),
)
    for (ti, t) in enumerate(time_points)
        rt = round(t, digits = 6)
        for (fi, f) in enumerate(frequencies)
            row = get(row_index, (round(f, digits = 6), rt), nothing)
            if !isnothing(row)
                df[row, channel] = mat[fi, ti]
            end
        end
    end

    return nothing
end
