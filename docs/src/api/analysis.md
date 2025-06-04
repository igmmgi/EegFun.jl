# Analysis

This section documents functions for EEG data analysis including ICA, time-frequency analysis, and statistical computations.

## Independent Component Analysis

### `infomax_ica(dat_ica, data_labels; n_components, params)`
Perform Independent Component Analysis using the Infomax algorithm.

**Arguments:**
- `dat_ica::Matrix{Float64}`: EEG data matrix
- `data_labels`: Channel labels
- `n_components::Union{Nothing,Int}`: Number of components (optional)
- `params::IcaPrms`: ICA parameters

### `create_work_arrays(n_components::Int, block_size::Int)`
Create working arrays for ICA computation.

**Returns:** WorkArrays structure containing temporary matrices

## Time-Frequency Analysis

### `generate_signal(n_trials, time_window, sample_rate, signal_freqs, signal_amps, signal_times, noise)`
Generate multiple trials of signals with different frequency components in specific time windows.

**Arguments:**
- `n_trials::Int`: Number of trials to generate
- `time_window::Vector{Float64}`: [start_time, end_time] window
- `sample_rate::Real`: Sampling rate in Hz
- `signal_freqs::Vector{<:Real}`: Vector of frequencies for each component
- `signal_amps::Vector{<:Real}`: Vector of amplitudes for each component
- `signal_times::Vector{Vector{Float64}}`: Vector of [start_time, end_time] for each component
- `noise::Real`: Amplitude of random noise

**Returns:** Matrix of size (n_samples × n_trials) containing generated signals

### `average_over_trials(x::Array{<:AbstractFloat,3})`
Average time-frequency data over trials.

### `apply_tf_baseline_db(tf_data, times, baseline_window)`
Apply baseline correction to time-frequency data using decibel conversion.

### `apply_tf_baseline_perchange(tf_data, times, baseline_window)`
Apply baseline correction to time-frequency data using percent change.

### `plot_tf(tf_dat, time, freq; kwargs...)`
Plot time-frequency data as a heatmap.

## Statistical Analysis

### `compute_probability!(probaMap, data, bins)`
Compute probability distribution for data.

**Arguments:**
- `probaMap::Vector{Float64}`: Output probability map
- `data::AbstractVector{Float64}`: The data vector to compute probabilities for
- `bins::Int`: Number of bins to use for probability computation

### `trim_extremes(x::Vector{Float64})`
Remove extreme values (top and bottom 10%) from data.

### `get_mean_amplitude(erp_data::ErpData, time_window)`
Calculate mean amplitude for each electrode within a specified time window.

**Returns:** DataFrame with electrode labels and corresponding mean amplitudes

## Event Analysis

### `_trigger_time_count(time, triggers)`
Internal function to process trigger data and count occurrences.

**Returns:** Trigger times, values, and count dictionary

### `plot_events(trigger_times, trigger_values, trigger_count)`
Plot trigger events as a scatter plot with vertical lines 