# Examples

This section provides practical examples of using eegfun for common EEG analysis tasks.

## Loading and Basic Preprocessing

```julia
using eegfun

# Load EEG data
data = load_eeg("subject01.bdf")

# Basic preprocessing pipeline
data = filter_data(data, 0.1, 40.0)  # Bandpass filter
data = reference_data(data, :average)  # Average reference
data = remove_bad_channels(data)       # Remove bad channels
```

## Event-Related Potentials

```julia
# Load event information
events = load_events("events.txt")

# Create epochs around events
epochs = epoch_data(data, events, (-0.2, 0.8))

# Compute ERPs
erp = compute_erp(epochs)

# Plot results
fig, ax = plot_erp(erp, ["Fz", "Cz", "Pz"])
```

## Time-Frequency Analysis

```julia
# Generate test signal
time_window = [-1.0, 2.0]
signal_freqs = [10.0, 20.0]
signal_amps = [1.0, 0.5]
signal_times = [[-0.5, 0.5], [0.0, 1.0]]

time, signals = generate_signal(
    100,  # 100 trials
    time_window,
    1000.0,  # 1000 Hz sampling rate
    signal_freqs,
    signal_amps,
    signal_times,
    0.1  # noise level
)

# Perform time-frequency analysis
tf_result = compute_tf(signals, time, freqs=1:50)

# Plot results
fig, ax = plot_tf(tf_result, time, 1:50)
```

## Independent Component Analysis

```julia
# Prepare data for ICA
ica_data = prepare_ica_data(data)

# Run ICA
ica_result = infomax_ica(ica_data, channel_labels)

# Plot component spectra
fig, ax = plot_components_spectra(ica_result, data, [1, 2, 3])
```

## Topographic Plotting

```julia
# Load electrode layout
layout = read_layout("biosemi64.csv")

# Create topographic plot of mean amplitude
mean_data = mean(data.data[!, 1:64], dims=1)
fig, ax = plot_topoplot(mean_data, layout)
```

## Automated Processing Pipeline

```julia
# Define preprocessing configuration
config = """
[preprocessing]
highpass = 0.1
lowpass = 40.0
reference = "average"
remove_bad_channels = true

[epochs]
prestim = -0.2
poststim = 0.8
baseline = [-0.2, 0.0]
"""

# Run automated preprocessing
preprocess_eeg_data("config.toml")
``` 