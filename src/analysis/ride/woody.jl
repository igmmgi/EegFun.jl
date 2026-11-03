"""
    RIDE Woody Filter Implementation

Native Julia implementation of the Woody filter for single-trial latency estimation.
Uses multi-threaded cross-correlation for massive performance improvements.
"""

using Statistics
using DSP
using Random

"""
    find_nearest_peak(signal::AbstractVector, center::Int)

Find the local maximum in `signal` that is closest to index `center`.
Returns `(peak_index, found_peak_boolean)`.
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

"""
    woody_filter(
        data::AbstractArray{T, 3}, 
        window_radius::Int;
        template::Union{Nothing, AbstractMatrix{T}} = nothing,
        initial_latencies::Union{Nothing, AbstractVector{Int}} = nothing,
        max_lag::Union{Nothing, Int} = nothing,
        rng::AbstractRNG = Random.default_rng(),
    ) where {T<:Real}

Estimate single-trial latencies using cross-covariance against a template.
If `template` is nothing, a leave-one-out average is used.

# Arguments
- `data`: 3D array of `(samples, channels, trials)`
- `window_radius`: How many samples from the center to search for a peak (equivalent to MATLAB `dur`)
- `template`: Optional fixed template `(samples, channels)`
- `initial_latencies`: Optional initial guesses for each trial (defaults to 0 shift)
- `max_lag`: Maximum lag for cross-correlation. Defaults to `size(data,1) ÷ 2` (expanded up to `size(data,1) - 1` if needed for `initial_latencies` + `window_radius`).
- `rng`: Random number generator for fallback jitter when a trial has no interior cross-covariance peak.

# Returns
- `latencies::Vector{Int}`: The estimated sample shifts for each trial relative to the median.
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
    length(init_lats) != n_trials && @minimal_error("initial_latencies length ($(length(init_lats))) does not match n_trials ($n_trials)")

    max_init_shift = isempty(init_lats) ? 0 : maximum(abs, init_lats)
    default_lag = min(n_samples - 1, max(div(n_samples, 2), max_init_shift + max(0, window_radius)))
    maxlag = isnothing(max_lag) ? default_lag : clamp(max_lag, 0, n_samples - 1)
    eff_radius = clamp(window_radius, 0, maxlag)

    lags = (-maxlag):maxlag
    n_lags = length(lags)
    center_idx = maxlag + 1

    latencies = copy(init_lats)
    found_peaks = falses(n_trials)

    # Calculate the grand average if no template is provided
    # We will subtract each trial from this sum to get the leave-one-out average
    grand_sum = zeros(T, n_samples, n_channels)
    if isnothing(template)
        @inbounds for t = 1:n_trials
            for c = 1:n_channels
                for s = 1:n_samples
                    grand_sum[s, c] += data[s, c, t]
                end
            end
        end
    end

    # Pre-allocate per-trial cross-covariance in parallel across trials
    Threads.@threads for t = 1:n_trials

        # 1. Get the template for this trial
        trial_template = zeros(T, n_samples, n_channels)
        if !isnothing(template)
            copyto!(trial_template, template)
        elseif n_trials > 1
            # Leave-one-out average
            @inbounds for c = 1:n_channels
                for s = 1:n_samples
                    trial_template[s, c] = (grand_sum[s, c] - data[s, c, t]) / (n_trials - 1)
                end
            end
        else
            copyto!(trial_template, grand_sum)
        end

        # 2. Compute cross-covariance across channels and average them
        mean_xcov = zeros(T, n_lags)

        for c = 1:n_channels
            # Extract trial data and template for this channel
            x = view(data, :, c, t)
            y = view(trial_template, :, c)

            # Mean center
            mx = mean(x)
            my = mean(y)
            x_centered = x .- mx
            y_centered = y .- my

            # Fast cross-correlation using DSP.xcorr
            # Note: DSP.xcorr returns 2N-1 points. We want the center 2*maxlag+1 points
            full_xcov = DSP.xcorr(x_centered, y_centered)

            # The center of full_xcov is at index n_samples
            # We extract from n_samples - maxlag to n_samples + maxlag
            start_idx = n_samples - maxlag

            @inbounds for l = 1:n_lags
                mean_xcov[l] += full_xcov[start_idx+l-1] / n_channels
            end
        end

        # 3. Extract the search window centered around initial latency guess
        search_center = center_idx + init_lats[t]
        search_center = clamp(search_center, 1 + eff_radius, n_lags - eff_radius)

        search_start = search_center - eff_radius
        search_end = search_center + eff_radius

        search_window = view(mean_xcov, search_start:search_end)

        # Find nearest peak relative to the center of the search window
        local_center = eff_radius + 1
        peak_idx, found = find_nearest_peak(search_window, local_center)

        found_peaks[t] = found

        # Convert local window peak index back to full lag index, then to latency
        if found
            full_peak_idx = search_start + peak_idx - 1
            latencies[t] = full_peak_idx - center_idx
        else
            latencies[t] = init_lats[t]
        end
    end

    # Deal with trials where no peak was found
    # Randomly assign them based on the distribution of found peaks (guarding <= 1 valid peak)
    valid_latencies = latencies[found_peaks]
    if length(valid_latencies) < n_trials
        if length(valid_latencies) > 1
            std_lat = std(valid_latencies)
            if isfinite(std_lat) && std_lat > 0
                missing_indices = findall(.!found_peaks)
                for idx in missing_indices
                    random_jitter = round(Int, randn(rng) * std_lat)
                    latencies[idx] = clamp(latencies[idx] + random_jitter, -maxlag, maxlag)
                end
            end
        elseif length(valid_latencies) == 1
            # Only 1 peak found: std would be NaN; assign missing trials to the found latency without NaN jitter
            missing_indices = findall(.!found_peaks)
            for idx in missing_indices
                latencies[idx] = valid_latencies[1]
            end
        end
    end

    # RIDE converts latencies to relative values by subtracting the median
    lat_median = round(Int, median(latencies))
    latencies .-= lat_median

    return latencies
end
