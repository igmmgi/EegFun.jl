"""
    generate_signal(
        n_trials::Int,
        time_window::Vector{Float64},
        sample_rate::Real,
        signal_freqs::Vector{<:Real},
        signal_amps::Vector{<:Real},
        signal_times::Vector{Vector{Float64}},
        noise::Real
    ) -> Matrix{Float64}

Generate multiple trials of signals with different frequency components in specific time windows.

# Arguments
- `n_trials`: Number of trials to generate
- `time_window`: [start_time, end_time] window
- `sample_rate`: Sampling rate in Hz
- `signal_freqs`: Vector of frequencies for each component
- `signal_amps`: Vector of amplitudes for each component
- `signal_times`: Vector of [start_time, end_time] for each component
- `noise`: Amplitude of random noise

# Returns
- Matrix of size (n_samples × n_trials) containing generated signals
"""
function generate_signal(
    n_trials::Int,
    time_window::Vector{Float64},
    sample_rate::Real,
    signal_freqs::Vector{<:Real},
    signal_amps::Vector{<:Real},
    signal_times::Vector{Vector{Float64}},
    noise::Real,
)
    # Pre-compute time vector
    dt = 1 / sample_rate
    x_time = range(time_window[1], time_window[2] - dt, step = dt)
    n_samples = length(x_time)

    # Pre-allocate output matrix with noise
    signal_trials = zeros(n_samples, n_trials)
    if noise > 0
        @views for trial = 1:n_trials
            signal_trials[:, trial] .= randn(n_samples) .* noise
        end
    end

    # Pre-allocate temporary arrays
    temp_signal = zeros(n_samples)

    # Add frequency components to all trials
    @inbounds for trial = 1:n_trials
        # Reset temporary signal array
        fill!(temp_signal, 0.0)

        # Add each frequency component
        for (freq, amp, times) in zip(signal_freqs, signal_amps, signal_times)
            # Find time window indices
            time_mask = (x_time .>= times[1]) .& (x_time .< times[2])
            window_times = @view x_time[time_mask]

            # Generate signal for this component
            @. temp_signal[time_mask] += amp * sin(2π * freq * window_times)
        end

        # Add to output
        signal_trials[:, trial] .+= temp_signal
    end

    return x_time, signal_trials
end


function average_over_trials(x::Array{<:AbstractFloat,3})
    return dropdims(mean(x, dims = 3), dims = 3)  # Returns a freq × time matrix
end

# baseline operations
function apply_tf_baseline_db(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return 10 .* log10.(tf_data ./ baseline_power)
end

function apply_tf_baseline_perchange(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return 100 .* (tf_data .- baseline_power) ./ baseline_power
end

function apply_tf_baseline_relchange(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(@view(tf_data[:, base_idx]), dims = 2)
    return tf_data ./ baseline_power
end

function apply_tf_baseline_ztransforms(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = @view tf_data[:, base_idx]
    baseline_mean = mean(baseline_power, dims = 2)
    baseline_std = std(baseline_power, dims = 2)
    return (tf_data .- baseline_mean) ./ baseline_std
end


function tf_morlet(signal, times, sample_rate, freqs, cycles; log_freqs = true, tois = nothing)

    # if vector
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal

    # data dimensions
    n_samples, n_trials = size(signal)
    n_freq = length(freqs)

    # time of interest
    tois_idx = 1:n_samples
    if !isnothing(tois)
        tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in tois]
    end

    # frequency steps
    if log_freqs
        freqs = exp.(range(log(freqs[1]), log(freqs[end]), length = n_freq))
    end

    # variable or fixed cycles/frequency
    if length(cycles) == 1
        cycles = cycles[1] ./ (2 * pi .* freqs) # variable cycles/frequency
    else
        cycles = exp.(range(log(cycles[1]), log(cycles[end]), length = n_freq)) ./ (2 * pi .* freqs)
    end

    # pre-compute signal FFT once maintining trial dimension
    # use lowest freq for max padding
    n_convolution = n_samples + ceil(Int, 6 * (maximum(cycles) / (2π * freqs[1])) * sample_rate)
    n_conv_pow2 = nextpow(2, n_convolution)
    signal_fft = fft([signal; zeros(n_conv_pow2 - n_samples, n_trials)], 1)  # FFT along time dimension

    tf_trials = zeros(n_freq, length(tois_idx), n_trials)
    @inbounds for idx_freq in eachindex(freqs)
        # frequency-dependent window length
        window_length = min(n_convolution, ceil(Int, 6 * cycles[idx_freq] * sample_rate))  # Cap window length
        window_length += iseven(window_length)  # ensure odd length

        #  wavelet
        wavelet_time = range(-window_length / 2, window_length / 2, length = window_length) ./ sample_rate
        wavelet =
            sqrt(1 ./ (cycles[idx_freq] .* sqrt(pi))) .* exp.(2im * pi * freqs[idx_freq] .* wavelet_time) .*
            exp.(-wavelet_time .^ 2 ./ (2 * (cycles[idx_freq] .^ 2)))

        # zero pad wavelet to match signal FFT length
        wavelet_fft = fft([wavelet; zeros(n_conv_pow2 - length(wavelet))])

        # convolution across trials
        eegconv = ifft(wavelet_fft .* signal_fft, 1)  # IFFT along time dimension

        # trim to original size and handle edge effects
        half_window = window_length ÷ 2
        valid_idx = (half_window+1):(half_window+n_samples)

        # store result
        @views tf_trials[idx_freq, :, :] = abs2.(eegconv[valid_idx, :])[tois_idx, :]

    end

    return tf_trials, times[tois_idx], freqs

end


function tf_hanning(signal, times, sample_rate, frequencies, time_steps; window_length = nothing, cycles = nothing)

    if isnothing(window_length) && isnothing(cycles)
        throw(ArgumentError("Must specify either window_length or cycles"))
    end

    # if vector
    signal = ndims(signal) == 1 ? reshape(signal, :, 1) : signal

    # convert time_steps to index values
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_steps]

    # Pre-allocate output and temporary arrays
    n_frex = length(frequencies)
    n_timepoints = length(time_steps)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)

    # Pre-compute maximum window length
    max_timewinidx = if !isnothing(cycles)
        round(Int, cycles / minimum(frequencies) * sample_rate)
    else
        round(Int, window_length * sample_rate)
    end
    max_timewinidx += iseven(max_timewinidx)

    # Pre-allocate reusable arrays
    tmpdat = zeros(max_timewinidx, n_trials)
    fdat = zeros(ComplexF64, max_timewinidx, n_trials)

    # Pre-compute window lengths and Hanning windows for each frequency
    timewin_indices = Vector{Int}(undef, n_frex)
    hann_windows = Vector{Vector{Float64}}(undef, n_frex)
    frex_indices = Vector{Int}(undef, n_frex)

    for (fi, freq) in enumerate(frequencies)
        timewinidx = if !isnothing(cycles)
            round(Int, cycles / freq * sample_rate)
        else
            round(Int, window_length * sample_rate)
        end
        timewinidx += iseven(timewinidx)
        timewin_indices[fi] = timewinidx

        # Pre-compute Hanning window
        hann_windows[fi] = 0.5 * (1 .- cos.(2π * (0:timewinidx-1) / (timewinidx - 1)))
        frex_indices[fi] = round(Int, freq * timewinidx / sample_rate) + 1
    end

    # Main computation loop
    @inbounds for fi = 1:n_frex
        timewinidx = timewin_indices[fi]
        hann_win = hann_windows[fi]
        frex_idx = frex_indices[fi]
        half_win = fld(timewinidx, 2)

        for (timepointi, center_idx) in enumerate(tois_idx)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win

            # Handle edge cases and data extraction
            fill!(view(tmpdat, 1:timewinidx, :), 0)
            if start_idx < 1 || end_idx > size(signal, 1)
                pad_left = max(1 - start_idx, 0)
                valid_start = max(start_idx, 1)
                valid_end = min(end_idx, size(signal, 1))
                valid_length = valid_end - valid_start + 1
                @views tmpdat[(pad_left+1):(pad_left+valid_length), :] .= signal[valid_start:valid_end, :]
            else
                @views tmpdat[1:timewinidx, :] .= signal[start_idx:end_idx, :]
            end

            # Apply Hanning taper and compute FFT
            @views tmpdat[1:timewinidx, :] .*= hann_win
            @views fdat[1:timewinidx, :] .= tmpdat[1:timewinidx, :]
            fft!(view(fdat, 1:timewinidx, :), 1)  # FFT along time dimension

            # Store power
            @views tf_trials[fi, timepointi, :] .= abs2.(fdat[frex_idx, :])
        end
    end

    return tf_trials, times[tois_idx], frequencies

end

# very basic tf plot
function plot_tf(tf_dat, time, freq; xlim=nothing, ylim=nothing, log_yscale=false, interpolate=false, colormap = :jet, colorrange=nothing, colorbar=false, colorbar_label = "")
  fig = Figure()
  # axes
  if log_yscale
    ax = Axis(fig[1, 1], yscale = log10)
  else
    ax = Axis(fig[1, 1])
  end
  # limits
  if isnothing(xlim)
    xlim = [time[1], time[end]]
  end
  if isnothing(ylim)
    ylim = [freq[1], freq[end]]
  end
  if isnothing(colorrange)
    colorrange = [minimum(tf_dat), maximum(tf_dat)]
  end
  hm = heatmap!(ax, times_out, freqs_out, transpose(tf_dat), colormap = colormap, interpolate = interpolate, colorrange=colorrange)
  xlims!(ax, xlim[1], xlim[end])
  ylims!(ax, ylim[1], ylim[end])
  ax.yticks = round.(freq);
  ax.ylabel = "Frequency [Hz]"
  ax.xlabel = "Time [S]"
  if colorbar
    Colorbar(fig[:, end+1], hm, label = colorbar_label, labelrotation=deg2rad(270))
  end
  display(fig)
  return fig, ax
end
