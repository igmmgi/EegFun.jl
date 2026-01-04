"""
    tf_morlet(dat::EpochData, channel::Symbol; 
              min_freq::Real=2, max_freq::Real=80, num_frex::Int=30,
              cycles::Union{Real,Tuple{Real,Real}}=(3, 10),
              tois::Union{Nothing,AbstractVector}=nothing)

Time-frequency analysis using Morlet wavelets (Cohen Chapter 13).

# Arguments
- `dat::EpochData`: Epoched EEG data
- `channel::Symbol`: Channel to analyze

# Keyword Arguments
- `min_freq::Real=2`: Minimum frequency in Hz
- `max_freq::Real=80`: Maximum frequency in Hz
- `num_frex::Int=30`: Number of frequencies (log-spaced)
- `cycles::Union{Real,Tuple{Real,Real}}=(3, 10)`: Number of cycles. Can be:
  - Single number: fixed cycles for all frequencies
  - Tuple (min, max): log-spaced cycles from min to max
- `tois::Union{Nothing,AbstractVector}=nothing`: Times of interest (seconds). 
  If `nothing`, uses all time points.

# Returns
- `TimeFreqData`: Time-frequency data object containing power (and phase, currently empty)

# Example
```julia
tf_data = tf_morlet(epochs, :Cz; min_freq=2, max_freq=80, num_frex=30)
```
"""
function tf_morlet(
    dat::EpochData,
    channel::Symbol;
    min_freq::Real = 2,
    max_freq::Real = 80,
    num_frex::Int = 30,
    cycles::Union{Real,Tuple{Real,Real}} = (3, 10),
    tois::Union{Nothing,AbstractVector} = nothing,
)
    # Validate channel exists
    ch_labels = channel_labels(dat)
    channel in ch_labels || error("Channel $channel not found in data. Available channels: $ch_labels")

    # Get sample rate
    sr = Float64(dat.sample_rate)

    # Get time vector from first epoch (all epochs should have same time points)
    times = dat.data[1][!, :time]
    n_samples = length(times)

    # Get number of trials/epochs
    n_trials = length(dat.data)

    # Extract signal data: reshape to (n_samples * n_trials, 1) for convolution
    # MATLAB: reshape(EEG.data(chan,:,:), 1, EEG.pnts*EEG.trials)
    signal = zeros(Float64, n_samples * n_trials)
    for (trial_idx, epoch_df) in enumerate(dat.data)
        signal[(trial_idx-1)*n_samples+1:trial_idx*n_samples] = epoch_df[!, channel]
    end

    # Define frequencies (log-spaced, matching MATLAB logspace)
    freqs = exp.(range(log(min_freq), log(max_freq), length = num_frex))

    # Define wavelet time window (MATLAB: -1:1/EEG.srate:1)
    wavelet_time = -1:(1/sr):1
    n_wavelet = length(wavelet_time)
    half_of_wavelet_size = (n_wavelet - 1) / 2

    # Define cycles/sigma
    # MATLAB: s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex)
    if cycles isa Tuple
        cycles_vec = exp.(range(log(cycles[1]), log(cycles[2]), length = num_frex))
    else
        cycles_vec = fill(Float64(cycles), num_frex)
    end
    sigma = cycles_vec ./ (2 * pi .* freqs)

    # Define convolution parameters
    n_data = n_samples * n_trials
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)

    # Get FFT of data (MATLAB: fft(reshape(...), n_conv_pow2))
    signal_padded = zeros(ComplexF64, n_conv_pow2)
    signal_padded[1:n_data] = signal
    eegfft = fft(signal_padded)

    # Initialize output power matrix (frequencies × time)
    eegpower = zeros(Float64, num_frex, n_samples)

    # Loop through frequencies
    for fi = 1:num_frex
        # Create wavelet (MATLAB: sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))))
        A = sqrt(1 / (sigma[fi] * sqrt(pi)))
        wavelet = A .* exp.(2im * pi * freqs[fi] .* wavelet_time) .* exp.(-wavelet_time .^ 2 ./ (2 * (sigma[fi]^2)))

        # FFT of wavelet
        wavelet_padded = zeros(ComplexF64, n_conv_pow2)
        wavelet_padded[1:n_wavelet] = wavelet
        wavelet_fft = fft(wavelet_padded)

        # Convolution (MATLAB: ifft(wavelet.*eegfft))
        eegconv = ifft(wavelet_fft .* eegfft)

        # Extract valid part (MATLAB: eegconv(1:n_convolution) then eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size))
        eegconv = eegconv[1:n_convolution]
        eegconv = eegconv[Int(half_of_wavelet_size)+1:end-Int(half_of_wavelet_size)]

        # Reshape to (n_samples × n_trials) and average power over trials
        # MATLAB: mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2)
        eegconv_reshaped = reshape(eegconv, n_samples, n_trials)
        eegpower[fi, :] = vec(mean(abs2.(eegconv_reshaped), dims = 2))
    end

    # Handle times of interest
    if !isnothing(tois)
        tois_idx = [argmin(abs.(times .- t)) for t in tois]
        times_out = times[tois_idx]
        eegpower = eegpower[:, tois_idx]
    else
        times_out = times
    end

    # Convert power matrix to DataFrame format (long format: time, freq, channel)
    # eegpower is (num_frex × n_times) where eegpower[fi, ti] = power at frequency fi, time ti
    # Baseline code expects: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, freq2_t2, ..., freqN_t2, ...]
    # Baseline calculates: row_idx = (ti - 1) * n_freqs + fi
    # So we need: for each time, all frequencies
    n_times = length(times_out)
    
    # Create DataFrame with all combinations of (time, freq)
    # Order: for each time, all freqs (time1: freq1, freq2, ..., freqN, time2: freq1, freq2, ...)
    power_df = DataFrame()
    power_df.time = repeat(times_out, inner = num_frex)  # Each time repeated for all freqs
    power_df.freq = repeat(freqs, outer = n_times)        # All freqs repeated for each time
    
    # Reshape power matrix to match expected order: [freq1_t1, freq2_t1, ..., freqN_t1, freq1_t2, ...]
    # eegpower is (num_frex × n_times), we need column-major order
    # vec(eegpower) gives: [freq1_time1, freq2_time1, ..., freqN_time1, freq1_time2, ...]
    # This matches the expected order!
    power_values = vec(eegpower)
    power_df[!, channel] = power_values

    # Create empty phase DataFrame with same structure
    phase_df = DataFrame()
    phase_df.time = copy(power_df.time)
    phase_df.freq = copy(power_df.freq)
    phase_df[!, channel] = fill(NaN, nrow(power_df))  # Phase not computed yet

    # Create TimeFreqData object
    return TimeFreqData(
        dat.file,
        dat.condition,
        dat.condition_name,
        power_df,
        phase_df,
        dat.layout,
        dat.sample_rate,
        :wavelet,
        nothing,  # baseline
        dat.analysis_info,
    )
end