"""
    generate_signal(n_trials, time_window, sample_rate, signal_freqs, signal_amps, signal_times, noise)

Generate synthetic signals with known frequency components for testing TF analysis.

# Arguments
- `n_trials`: Number of trials to generate
- `time_window`: [start_time, end_time] in seconds
- `sample_rate`: Sampling rate in Hz
- `signal_freqs`: Vector of frequencies for each component
- `signal_amps`: Vector of amplitudes for each component
- `signal_times`: Vector of [start_time, end_time] for each component
- `noise`: Amplitude of random noise

# Returns
- `times`: Time vector
- `signal`: Matrix (samples × trials)
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

    # Find time 0 in the time array (same for all trials and frequencies)
    # signal_times are relative to time 0 (event/stimulus time), not time_window[1]
    time_zero_idx = argmin(abs.(x_time))
    time_zero = x_time[time_zero_idx]

    # Add frequency components to all trials
    @inbounds for trial = 1:n_trials
        # Reset temporary signal array
        fill!(temp_signal, 0.0)

        # Add each frequency component
        for (freq, amp, times) in zip(signal_freqs, signal_amps, signal_times)
            # Convert signal_times from relative to time 0 to absolute times
            abs_times = [time_zero + times[1], time_zero + times[2]]
            
            # Find time window indices (now using absolute times)
            time_mask = (x_time .>= abs_times[1]) .& (x_time .< abs_times[2])
            window_times = @view x_time[time_mask]

            # Generate signal for this component
            @. temp_signal[time_mask] += amp * sin(2π * freq * window_times)
        end

        # Add to output
        signal_trials[:, trial] .+= temp_signal
    end

    return x_time, signal_trials
end

"""
    _signal_to_epochs(times, signal, channel_name::Symbol, sample_rate::Int;
                     file::String="synthetic", condition::Int=1, condition_name::String="test") -> EpochData

Internal function: Convert signal data from `generate_signal` to `EpochData` format.

Used internally for testing. Not part of the public API.
"""
function _signal_to_epochs(
    times,
    signal,
    channel_name::Symbol,
    sample_rate::Int;
    file::String = "synthetic",
    condition::Int = 1,
    condition_name::String = "test",
)
    # Ensure signal is 2D (samples × trials)
    if ndims(signal) == 1
        signal = reshape(signal, :, 1)
    end

    n_samples, n_trials = size(signal)

    # Convert times to vector if needed
    times_vec = collect(Float64, times)
    @assert length(times_vec) == n_samples "Time vector length ($(length(times_vec))) must match signal samples ($n_samples)"

    # Create Vector{DataFrame} - one DataFrame per trial
    trial_dfs = Vector{DataFrame}(undef, n_trials)

    for trial_idx = 1:n_trials
        # Create DataFrame for this trial with time and channel columns
        trial_dfs[trial_idx] = DataFrame(:time => times_vec, channel_name => signal[:, trial_idx], :epoch => trial_idx)
    end

    # Create minimal layout for the channel
    layout_df = DataFrame(:label => [channel_name])
    layout = Layout(layout_df, nothing, nothing)

    # Create AnalysisInfo
    analysis_info = AnalysisInfo()

    # Create and return EpochData
    return EpochData(file, condition, condition_name, trial_dfs, layout, sample_rate, analysis_info)
end
