"""
    _tf_build_row_index(df) -> Dict{Tuple{Float64,Float64}, Int}

Build a lookup from (rounded_freq, rounded_time) → row index for a TF DataFrame.
Uses `round(x, digits=6)` for robust floating-point matching.
"""
function _tf_build_row_index(df::DataFrame)
    idx = Dict{Tuple{Float64,Float64},Int}()
    n = nrow(df)
    sizehint!(idx, n)
    freq_col = df.freq
    time_col = df.time
    for i = 1:n
        key = (round(freq_col[i], digits = 6), round(time_col[i], digits = 6))
        if !haskey(idx, key)
            idx[key] = i
        end
    end
    return idx
end

"""
    _tf_df_to_matrix!(mat, df, channel, rounded_freqs, rounded_times, row_index)

In-place extraction of a `[n_freqs × n_time]` power matrix for one channel from a TF DataFrame.
Writes into pre-allocated `mat` which should be initialized (e.g. with `fill!(mat, NaN)`).
"""
function _tf_df_to_matrix!(
    mat::AbstractMatrix{Float64},
    df::DataFrame,
    channel::Symbol,
    rounded_freqs::Vector{Float64},
    rounded_times::Vector{Float64},
    row_index::Dict{Tuple{Float64,Float64},Int},
)
    col_data = df[!, channel]::Vector{Float64}
    @inbounds for (ti, rt) in enumerate(rounded_times)
        for (fi, rf) in enumerate(rounded_freqs)
            row = get(row_index, (rf, rt), nothing)
            if !isnothing(row)
                mat[fi, ti] = col_data[row]
            end
        end
    end
    return mat
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
    row_index::Dict{Tuple{Float64,Float64},Int} = _tf_build_row_index(df);
    rounded_freqs::Vector{Float64} = round.(frequencies, digits = 6),
    rounded_times::Vector{Float64} = round.(time_points, digits = 6),
)
    mat = fill(NaN, length(frequencies), length(time_points))
    return _tf_df_to_matrix!(mat, df, channel, rounded_freqs, rounded_times, row_index)
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
    col = df[!, channel]::Vector{Float64}
    for (ti, t) in enumerate(time_points)
        rt = round(t, digits = 6)
        for (fi, f) in enumerate(frequencies)
            row = get(row_index, (round(f, digits = 6), rt), nothing)
            if !isnothing(row)
                col[row] = mat[fi, ti]
            end
        end
    end

    return nothing
end

"""
    _filter_edges!(eegpower, eegconv, num_frex, time_indices, window_lengths_per_freq, n_samples_per_epoch)

Set edge regions to NaN where the analysis window extends beyond the data boundaries.
Applied to both power and (optionally) complex convolution results.
"""
function _filter_edges!(
    eegpower::AbstractArray,
    eegconv::Union{AbstractArray,Nothing},
    num_frex::Int,
    time_indices::AbstractVector{Int},
    window_lengths_per_freq::Union{Vector{Int},Vector{Float64}},
    n_samples_per_epoch::Int,
)
    for fi = 1:num_frex

        half_nsamplefreqoi = window_lengths_per_freq[fi] / 2.0
        min_valid_threshold = half_nsamplefreqoi
        max_valid_threshold = n_samples_per_epoch - half_nsamplefreqoi

        for ti in eachindex(time_indices)
            sample_idx = time_indices[ti]
            if !(sample_idx >= min_valid_threshold && sample_idx < max_valid_threshold)
                if ndims(eegpower) == 3  # return_trials = true, layout is [trial, freq, time]
                    @views eegpower[:, fi, ti] .= NaN
                    if !isnothing(eegconv)
                        @views eegconv[:, fi, ti] .= NaN * im
                    end
                else  # return_trials = false
                    eegpower[fi, ti] = NaN
                    if !isnothing(eegconv)
                        eegconv[fi, ti] = NaN * im
                    end
                end
            end
        end
    end
end

"""
    _validate_tf_pad(pad)

Throw an error unless `pad` is `nothing`, `:pre`, `:post` or `:both`.
"""
function _validate_tf_pad(pad)
    if !isnothing(pad) && pad ∉ [:pre, :post, :both]
        @minimal_error("`pad` must be `nothing`, `:pre`, `:post`, or `:both`, got :$pad")
    end
    return nothing
end

"""
    _tf_gpu_active(use_gpu::Bool) -> Bool

Return whether the GPU should be used. Warns and returns `false` if a GPU was
requested but CUDA is not available or functional.
"""
function _tf_gpu_active(use_gpu::Bool)
    gpu_active = use_gpu && is_gpu_available()
    if use_gpu && !gpu_active
        @minimal_warning "GPU requested but CUDA is not available or functional. Falling back to CPU."
    end
    return gpu_active
end

"""
    _subset_tf_data(dat::EpochData, channel_selection, interval_selection) -> EpochData

Subset `dat` by channel and interval selection, erroring if no data remain.
"""
function _subset_tf_data(dat::EpochData, channel_selection::Function, interval_selection::Interval)
    dat = subset(dat; channel_selection = channel_selection, interval_selection = interval_selection)
    isempty(dat.data) && @minimal_error("No data remaining after subsetting")
    return dat
end

"""
    _validate_tf_frequencies(frequencies)

Throw an error unless `frequencies` is non-empty and strictly positive.
"""
function _validate_tf_frequencies(frequencies)
    if length(frequencies) == 0
        @minimal_error("`frequencies` must contain at least one frequency")
    end
    if any(f -> f <= 0, frequencies)
        @minimal_error("All frequencies in `frequencies` must be positive")
    end
    return nothing
end

"""
    _validate_tf_window_samples(n_window_samples)

Throw an error if any window (given in samples, scalar or vector) is shorter than 2 samples.
"""
function _validate_tf_window_samples(n_window_samples)
    min_samples = minimum(n_window_samples)
    if min_samples < 2
        @minimal_error("Window length is too short. Minimum window size is 2 samples, got $min_samples samples.")
    end
    return nothing
end

"""
    _tf_time_grid(dat::EpochData, time_steps::Real) -> (time_indices, times_out)

Determine the output time grid (sample indices and time values) for a step size of
`time_steps` seconds. Errors if no valid time points are found.
"""
function _tf_time_grid(dat::EpochData, time_steps::Real)
    times = time_vector(dat)
    time_min = minimum(times)
    time_max = maximum(times)
    time_steps_range = time_min:time_steps:time_max
    time_indices, times_out = find_times(times, time_steps_range)
    if isempty(time_indices)
        @minimal_error("No valid time points found with step size $time_steps in range ($time_min to $time_max seconds)")
    end
    return time_indices, times_out
end
