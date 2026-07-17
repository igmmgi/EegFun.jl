"""
    RIDE Woody Filter Implementation

Native Julia implementation of the Woody filter for single-trial latency estimation.
Uses multi-threaded cross-correlation for massive performance improvements.
"""

using Statistics
using DSP

"""
    find_nearest_peak(signal::AbstractVector, center::Int)

Find the local maximum in `signal` that is closest to index `center`.
Returns `(peak_index, found_peak_boolean)`.
"""
function find_nearest_peak(signal::AbstractVector, center::Int)
    n = length(signal)
    best_dist = typemax(Int)
    best_idx = center
    found_peak = false
    
    @inbounds for i in 2:(n-1)
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
        max_lag::Union{Nothing, Int} = nothing
    ) where {T<:Real}

Estimate single-trial latencies using cross-covariance against a template.
If `template` is nothing, a leave-one-out average is used.

# Arguments
- `data`: 3D array of `(samples, channels, trials)`
- `window_radius`: How many samples from the center to search for a peak (equivalent to MATLAB `dur`)
- `template`: Optional fixed template `(samples, channels)`
- `initial_latencies`: Optional initial guesses for each trial (defaults to 0 shift)
- `max_lag`: Maximum lag for cross-correlation. Defaults to `size(data,1) ÷ 2`.

# Returns
- `latencies::Vector{Int}`: The estimated sample shifts for each trial relative to the median.
"""
function woody_filter(
    data::AbstractArray{T, 3}, 
    window_radius::Int;
    template::Union{Nothing, AbstractMatrix{T}} = nothing,
    initial_latencies::Union{Nothing, AbstractVector{Int}} = nothing,
    max_lag::Union{Nothing, Int} = nothing
) where {T<:Real}

    n_samples, n_channels, n_trials = size(data)
    
    maxlag = isnothing(max_lag) ? div(n_samples, 2) : max_lag
    lags = -maxlag:maxlag
    n_lags = length(lags)
    center_idx = maxlag + 1
    
    if isnothing(initial_latencies)
        initial_latencies = zeros(Int, n_trials)
    end
    
    latencies = copy(initial_latencies)
    found_peaks = falses(n_trials)
    
    # Calculate the grand average if no template is provided
    # We will subtract each trial from this sum to get the leave-one-out average
    grand_sum = zeros(T, n_samples, n_channels)
    if isnothing(template)
        @inbounds for t in 1:n_trials
            for c in 1:n_channels
                for s in 1:n_samples
                    grand_sum[s, c] += data[s, c, t]
                end
            end
        end
    end
    
    # Pre-allocate thread-local buffers to prevent GC pressure in the multithreaded loop
    # We need a buffer for the cross-covariance of one channel and the mean across channels
    
    Threads.@threads for t in 1:n_trials
        
        # 1. Get the template for this trial
        trial_template = zeros(T, n_samples, n_channels)
        if !isnothing(template)
            copyto!(trial_template, template)
        else
            # Leave-one-out average
            @inbounds for c in 1:n_channels
                for s in 1:n_samples
                    trial_template[s, c] = (grand_sum[s, c] - data[s, c, t]) / (n_trials - 1)
                end
            end
        end
        
        # 2. Compute cross-covariance across channels and average them
        mean_xcov = zeros(T, n_lags)
        
        for c in 1:n_channels
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
            
            @inbounds for l in 1:n_lags
                mean_xcov[l] += full_xcov[start_idx + l - 1] / n_channels
            end
        end
        
        # 3. Extract the search window
        # The search window is centered around the initial latency guess
        # initial_latencies[t] is the expected shift.
        # In MATLAB: nearest_latency(temp5, initial_latency + center)
        
        search_center = center_idx + initial_latencies[t]
        search_center = clamp(search_center, 1 + window_radius, n_lags - window_radius)
        
        search_start = search_center - window_radius
        search_end = search_center + window_radius
        
        search_window = view(mean_xcov, search_start:search_end)
        
        # Find nearest peak relative to the center of the search window
        # Center of search_window is at index (window_radius + 1)
        local_center = window_radius + 1
        peak_idx, found = find_nearest_peak(search_window, local_center)
        
        found_peaks[t] = found
        
        # Convert local window peak index back to full lag index, then to latency
        if found
            full_peak_idx = search_start + peak_idx - 1
            latencies[t] = full_peak_idx - center_idx
        else
            latencies[t] = initial_latencies[t]
        end
    end
    
    # Deal with trials where no peak was found
    # MATLAB randomly assigns them based on the distribution of found peaks
    valid_latencies = latencies[found_peaks]
    if !isempty(valid_latencies) && length(valid_latencies) < n_trials
        std_lat = std(valid_latencies)
        missing_indices = findall(.!found_peaks)
        for idx in missing_indices
            random_jitter = round(Int, randn() * std_lat)
            latencies[idx] = latencies[idx] + random_jitter
        end
    end
    
    # RIDE usually converts latencies to relative values by subtracting the median
    lat_median = round(Int, median(latencies))
    latencies .-= lat_median
    
    return latencies
end
