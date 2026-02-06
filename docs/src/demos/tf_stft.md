# TF STFT

This demo demonstrates time-frequency analysis using the Short-Time Fourier Transform (STFT), a classic approach to spectrograms.

This demo demonstrates time-frequency analysis using the Short-Time Fourier Transform (STFT), a classic approach to spectrograms.

### What is STFT?

STFT applies the Fast Fourier Transform (FFT) to windowed segments of the signal:

1. **Window the signal**: Apply time-limited window (e.g., Hann, Tukey)
2. **Compute FFT**: Calculate frequency spectrum for that segment
3. **Slide window**: Move forward in time and repeat
4. **Build spectrogram**: Stack spectra to create time-frequency representation

### Key Parameters

**Window length** (`window_length`):

Controls time-frequency resolution:

```julia
tf_stft(epochs, window_length = 0.5)  # 500 ms windows
```

- **Longer windows** (e.g., 500 ms): Better frequency resolution, worse time resolution
- **Shorter windows** (e.g., 200 ms): Better time resolution, worse frequency resolution

**Cycles** (alternative):

Instead of specifying window length, specify cycles at each frequency:

```julia
tf_stft(epochs, cycles = 7)
```

Creates frequency-adaptive windows similar to Morlet wavelets.

**Time steps**:

Controls temporal resolution of output:

```julia
tf_stft(epochs, time_steps = 0.005)  # 5 ms steps
```

Smaller = smoother time course, larger = faster computation.

### Frequency Spacing

**Linear** (`frequencies = 2:1:80`):

- Equal Hz spacing
- Natural for narrow bands

**Logarithmic** (`frequencies = logrange(2, 80, length = 100)`):

- Proportional spacing
- Better for wide ranges
- Use `ylogscale = true` in plots

### Baseline Correction

```julia
plot_time_frequency(
    tf_data,
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3)
)
```

Shows relative power changes from pre-stimulus baseline.

### Edge Effects

**Filter edges** (`filter_edges = true`):

- Removes edge samples where windowing causes artifacts
- Recommended for cleaner results

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz)
2. **Fixed window length**: Constant resolution analysis
3. **Cycles mode**: Adaptive resolution
4. **Linear vs log spacing**: Different frequency sampling
5. **Baseline correction**: Event-related power changes
6. **Edge filtering**: Removing edge artifacts

### Further Reading

Cohen, M. X. (2014). *Analyzing Neural Time Series Data: Theory and Practice*. Chapter 15: Short-Time FFT


## Code Examples

::: details Show Code

```julia
using EegFun

#######################################################################
@info EegFun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal
sample_rate = 1000.0
times, signal = EegFun.generate_signal(
    400,                                    # n_trials
    [-1.0, 3.0],                            # time_window
    sample_rate,                            # sample_rate
    [5.0, 25, 35.0],                        # frequencies
    [5.0, 5.0, 5.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time windows for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = EegFun.signal_to_data(times, signal, :Channel1, sample_rate)
EegFun.plot_epochs(epochs_synthetic, channel_selection = EegFun.channels([:Channel1]))

spectrum = EegFun.freq_spectrum(epochs_synthetic, max_freq = 80.0)
EegFun.plot_frequency_spectrum(spectrum, channel_selection = EegFun.channels([:Channel1]))

# tf_stft_fixed
tf_data = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_stft(epochs_synthetic, frequencies = logrange(1, 40, length = 30), window_length = 0.5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

tf_data = EegFun.tf_stft(epochs_synthetic, frequencies = logrange(1, 40, length = 30), window_length = 0.5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

#######################################################################
@info EegFun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################
# This is some data that was presented in Cohen: Analyzin Neural Time Series Data
data_cohen = EegFun.read_data("./data/files/tf_test_epochs.jld2");

# Figure 13.11 A)
tf_data =
    EegFun.tf_stft(data_cohen, frequencies = logrange(2, 80, length = 100), window_length = 0.5, time_steps = 0.001, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.11 A)
tf_data = EegFun.tf_stft(data_cohen, frequencies = 2:1:80, cycles = 7, time_steps = 0.005, filter_edges = false)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)

# Figure 13.14 A)
tf_data =
    EegFun.tf_stft(data_cohen, frequencies = logrange(2, 80, length = 30), window_length = 0.5, time_steps = 0.005, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 B)
tf_data =
    EegFun.tf_stft(data_cohen, frequencies = logrange(2, 80, length = 30), window_length = 0.5, time_steps = 0.005, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 C)
tf_data = EegFun.tf_stft(data_cohen, frequencies = logrange(2, 80, length = 30), cycles = 5, time_steps = 0.005, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)
```

:::

## See Also

- [API Reference](../reference/index.md)
