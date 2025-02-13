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

# Generate 100 trials with multiple frequency components
# times, signal = generate_signal(
#     100,                             # n_trials
#     [-2.0, 3.0],                     # time window
#     256.0,                           # sampling rate
#     [10.0, 15.0, 25.0, 40.0],        # frequencies
#     [1.0, 2.0, 1.0, 1.0],            # amplitudes
#     [
#         [0.0, 0.25],
#         [0.25, 0.5],      # time windows
#         [0.5, 1.0],
#         [1.0, 1.5],
#     ],
#     0,                              # noise level
# )
# 
# lines(times, vec(mean(signal, dims = 2)))


function tf_morlet(signal, times, sample_rate, freqs, cycles; log_freqs = true, tois = nothing)

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
        wavelet_time = range(-window_length/2, window_length/2, length=window_length) ./ sample_rate
        wavelet = sqrt(1 ./ (cycles[idx_freq] .* sqrt(pi))) .* 
                 exp.(2im * pi * freqs[idx_freq] .* wavelet_time) .*
                 exp.(-wavelet_time .^ 2 ./ (2 * (cycles[idx_freq] .^ 2)))
        
        # zero pad wavelet to match signal FFT length
        wavelet_fft = fft([wavelet; zeros(n_conv_pow2 - length(wavelet))])
        
        # convolution across trials
        eegconv = ifft(wavelet_fft .* signal_fft, 1)  # IFFT along time dimension
        
        # trim to original size and handle edge effects
        half_window = window_length ÷ 2
        valid_idx = (half_window + 1):(half_window + n_samples)
        
        # store result
        @views tf_trials[idx_freq, :, :] = abs2.(eegconv[valid_idx, :])[tois_idx, :]
    end

    return tf_trials, times[tois_idx], freqs

end



function average_over_trials(x::Array{<:AbstractFloat, 3})
    return dropdims(mean(tf_data, dims=3), dims=3)  # Returns a freq × time matrix
end


# Figure 13.11 Panel A
signal = DataFrame(CSV.File("data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1:(1/sample_rate):1.5-(1/sample_rate)) 
#tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80, [10])
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80, [3, 10]; log_freqs = false, tois = -1:0.01:1)
tf = average_trials(tf_trials)
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))

fig = Figure()
ax = Axis(fig[1, 1], yscale = log10)
# contourf!(ax, times_out, freqs_out, transpose(tf), levels = -5:0.2:5, colormap = :jet)
heatmap!(ax, times_out, freqs_out, transpose(tf), colormap = :jet, interpolate = true)
xlims!(ax, -0.2, 1.0)
ylims!(ax, 2, 80)
ax.yticks = round.(freqs_out)
display(fig)



function tf_hanning(signal, times, sample_rate, frequencies, time_steps; 
                   window_length=nothing, cycles=nothing)
    if isnothing(window_length) && isnothing(cycles)
        throw(ArgumentError("Must specify either window_length or cycles"))
    end
    
    # Convert time_steps to indices
    tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_steps]
    
    # Pre-allocate output and temporary arrays
    n_frex = length(frequencies)
    n_timepoints = length(time_steps)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)
    
    # Pre-compute maximum window length for pre-allocation
    max_timewinidx = if !isnothing(cycles)
        round(Int, cycles / minimum(frequencies) * sample_rate)
    else
        round(Int, window_length * sample_rate)
    end
    max_timewinidx += iseven(max_timewinidx)
    
    # Pre-allocate reusable arrays outside all loops
    tmpdat = zeros(max_timewinidx, n_trials)
    fdat = zeros(ComplexF64, max_timewinidx, n_trials)
    
    for (fi, freq) in enumerate(frequencies)
        # Calculate frequency-dependent or fixed window length
        timewinidx = if !isnothing(cycles)
            round(Int, cycles / freq * sample_rate)
        else
            round(Int, window_length * sample_rate)
        end
        timewinidx += iseven(timewinidx)
        
        # Pre-compute window and frequency index for this frequency
        hann_win = 0.5 * (1 .- cos.(2π * (0:timewinidx-1) / (timewinidx - 1)))
        frex_idx = round(Int, freq * timewinidx / sample_rate) + 1
        
        for (timepointi, center_idx) in enumerate(tois_idx)
            # Calculate window indices
            half_win = fld(timewinidx, 2)
            start_idx = center_idx - half_win
            end_idx = center_idx + half_win
            
            # Handle edge cases and data extraction in one step
            fill!(view(tmpdat, 1:timewinidx, :), 0)
            if start_idx < 1 || end_idx > size(signal, 1)
                pad_left = max(1 - start_idx, 0)
                valid_start = max(start_idx, 1)
                valid_end = min(end_idx, size(signal, 1))
                valid_length = valid_end - valid_start + 1
                tmpdat[(pad_left + 1):(pad_left + valid_length), :] = view(signal, valid_start:valid_end, :)
            else
                tmpdat[1:timewinidx, :] = view(signal, start_idx:end_idx, :)
            end
            
            # Apply Hanning taper and compute FFT in-place
            tmpdat[1:timewinidx, :] .*= hann_win
            fdat[1:timewinidx, :] = tmpdat[1:timewinidx, :]
            fft!(view(fdat, 1:timewinidx, :), 1)  # FFT along time dimension
            
            # Store power for all trials at once
            tf_trials[fi, timepointi, :] = abs2.(fdat[frex_idx, :]) ./ timewinidx^2
        end
    end
    
    return tf_trials, times[tois_idx], frequencies
end


signal = DataFrame(CSV.File("data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1:(1/sample_rate):1.5-(1/sample_rate)) 
# tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:80, -0.5:0.2:1.5, window_length = 0.2)
tf_trials, times_out, freqs_out = tf_hanning(signal, times, sample_rate, 1:1:80, -0.5:0.2:1.5, cycles = 7)
tf = average_trials(tf_trials)
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
fig = Figure()
ax = Axis(fig[1, 1])#, yscale = log10)
# heatmap!(ax, times_out, freqs_out, transpose(tf), colormap = :jet, interpolate = true)
contourf!(ax, times_out, freqs_out, transpose(tf), levels = -6:0.2:6, colormap = :jet)
xlims!(ax, -0.5, 1.5)
ylims!(ax, 2, 80) 
ax.yticks = round.(freqs_out)
display(fig)


