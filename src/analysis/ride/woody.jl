"""
    find_nearest_peak(signal::AbstractVector, center::Int)

Find the strict local maximum in `signal` closest to index `center`, returning `(peak_index, found_peak)`.
"""
function find_nearest_peak(signal::AbstractVector, center::Int)
    n = length(signal)
    best_dist = typemax(Int)
    best_idx = clamp(center, 1, max(1, n))
    found_peak = false

    @inbounds for i = 2:(n-1)
        if signal[i] > signal[i-1] && signal[i] > signal[i+1]
            found_peak = true
            dist = abs(i - center)
            if dist < best_dist
                best_dist = dist
                best_idx = i
            end
        end
    end

    return best_idx, found_peak
end

@inline function _center_vector!(dest::Vector{Float64}, src::AbstractVector)
    n = length(src)
    s = 0.0
    @inbounds @simd for i = 1:n
        v = Float64(src[i])
        dest[i] = v
        s += v
    end
    m = s / n
    @inbounds @simd for i = 1:n
        dest[i] -= m
    end
    return dest
end

@inline function _accumulate_xcov_window!(
    win_xcov::Vector{Float64},
    x::Vector{Float64},
    y::AbstractVector{Float64},
    lag_min::Int,
    lag_max::Int,
)
    n = length(x)
    @inbounds for idx = 1:(lag_max-lag_min+1)
        τ = lag_min + idx - 1
        s = 0.0
        if τ >= 0
            @simd for k = 1:(n-τ)
                s += x[k+τ] * y[k]
            end
        else
            a = -τ
            @simd for k = 1:(n-a)
                s += x[k] * y[k+a]
            end
        end
        win_xcov[idx] += s
    end
    return nothing
end

"""
    woody_filter(data, window_radius; template=nothing, initial_latencies=nothing,
                 max_lag=nothing, rng=Random.default_rng()) -> Vector{Int}

Estimate single-trial sample shifts via multi-channel cross-covariance against `template`
(or leave-one-out trial averages when `template === nothing`), centered around the median shift.
"""
function woody_filter(
    data::AbstractArray{T,3},
    window_radius::Int;
    template::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    initial_latencies::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    max_lag::Union{Nothing,Int} = nothing,
    rng::AbstractRNG = Random.default_rng(),
) where {T<:Real}
    n_samples, n_channels, n_trials = size(data)
    n_trials == 0 && return Int[]
    if n_samples <= 2 || n_channels == 0 || (n_trials == 1 && isnothing(template))
        return zeros(Int, n_trials)
    end

    init_lats = isnothing(initial_latencies) ? zeros(Int, n_trials) : collect(Int, initial_latencies)
    length(init_lats) == n_trials || @minimal_error("initial_latencies length ($(length(init_lats))) does not match n_trials ($n_trials)")

    max_init_shift = isempty(init_lats) ? 0 : maximum(abs, init_lats)
    default_lag = min(n_samples - 1, max(div(n_samples, 2), max_init_shift + max(0, window_radius)))
    maxlag = isnothing(max_lag) ? default_lag : clamp(max_lag, 0, n_samples - 1)
    eff_radius = clamp(window_radius, 0, maxlag)

    n_lags = 2 * maxlag + 1
    center_idx = maxlag + 1
    latencies = copy(init_lats)
    found_peaks = Vector{Bool}(undef, n_trials)

    # Precompute column-centered template (or centered grand sum for leave-one-out)
    tmpl_cent = Matrix{Float64}(undef, n_samples, n_channels)
    has_fixed_template = !isnothing(template)
    if has_fixed_template
        @inbounds for c = 1:n_channels
            s_col = 0.0
            @simd for s = 1:n_samples
                v = Float64(template[s, c])
                tmpl_cent[s, c] = v
                s_col += v
            end
            m_col = s_col / n_samples
            @simd for s = 1:n_samples
                tmpl_cent[s, c] -= m_col
            end
        end
    else
        fill!(tmpl_cent, 0.0)
        @inbounds for t = 1:n_trials, c = 1:n_channels
            @simd for s = 1:n_samples
                tmpl_cent[s, c] += Float64(data[s, c, t])
            end
        end
        @inbounds for c = 1:n_channels
            s_col = 0.0
            @simd for s = 1:n_samples
                s_col += tmpl_cent[s, c]
            end
            m_col = s_col / n_samples
            @simd for s = 1:n_samples
                tmpl_cent[s, c] -= m_col
            end
        end
    end

    inv_loo = (!has_fixed_template && n_trials > 1) ? 1.0 / (n_trials - 1) : 1.0
    n_win = 2 * eff_radius + 1

    _foreach_chunk(n_trials) do trial_range
        x_cent = Vector{Float64}(undef, n_samples)
        y_cent = has_fixed_template ? Float64[] : Vector{Float64}(undef, n_samples)
        win_xcov = Vector{Float64}(undef, n_win)

        for t in trial_range
            search_center = clamp(center_idx + init_lats[t], 1 + eff_radius, n_lags - eff_radius)
            search_start = search_center - eff_radius
            search_end = search_center + eff_radius
            lag_min = search_start - center_idx
            lag_max = search_end - center_idx

            fill!(win_xcov, 0.0)
            for c = 1:n_channels
                _center_vector!(x_cent, @view(data[:, c, t]))
                if has_fixed_template
                    _accumulate_xcov_window!(win_xcov, x_cent, @view(tmpl_cent[:, c]), lag_min, lag_max)
                else
                    if n_trials > 1
                        @inbounds @simd for s = 1:n_samples
                            y_cent[s] = (tmpl_cent[s, c] - x_cent[s]) * inv_loo
                        end
                        _accumulate_xcov_window!(win_xcov, x_cent, y_cent, lag_min, lag_max)
                    else
                        _accumulate_xcov_window!(win_xcov, x_cent, @view(tmpl_cent[:, c]), lag_min, lag_max)
                    end
                end
            end

            peak_idx, found = find_nearest_peak(win_xcov, eff_radius + 1)
            found_peaks[t] = found
            latencies[t] = found ? (search_start + peak_idx - 1 - center_idx) : init_lats[t]
        end
    end

    n_valid = count(found_peaks)
    if n_valid < n_trials
        if n_valid > 1
            valid_latencies = [latencies[t] for t = 1:n_trials if found_peaks[t]]
            std_lat = std(valid_latencies)
            if isfinite(std_lat) && std_lat > 0
                for t = 1:n_trials
                    if !found_peaks[t]
                        jitter = round(Int, randn(rng) * std_lat)
                        latencies[t] = clamp(latencies[t] + jitter, -maxlag, maxlag)
                    end
                end
            end
        elseif n_valid == 1
            only_lat = latencies[findfirst(found_peaks)::Int]
            for t = 1:n_trials
                !found_peaks[t] && (latencies[t] = only_lat)
            end
        end
    end

    lat_median = round(Int, median(latencies))
    latencies .-= lat_median
    return latencies
end
