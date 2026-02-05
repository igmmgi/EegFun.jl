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

### Window Types

The STFT uses **Tukey windows** by default:
- Smooth edges reduce spectral leakage
- Similar to Hann windows
- Good compromise between resolution and leakage

### Fixed vs Adaptive Resolution

**Window length mode** (fixed resolution):
- Same time-frequency resolution at all frequencies
- Simple and predictable
- Good for narrow frequency ranges

**Cycles mode** (adaptive resolution):
- Resolution adapts with frequency  
- Longer windows at low frequencies
- Shorter windows at high frequencies
- Better for wide frequency ranges

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

### When to Use STFT

**Best for**:
- Quick exploratory analyses
- When constant resolution is acceptable
- Broad-band power changes
- Familiar, widely understood method

**Consider alternatives**:
- **Morlet wavelets**: Need adaptive resolution across frequencies
- **Multitaper**: Need best spectral quality and variance reduction

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz)
2. **Fixed window length**: Constant resolution analysis
3. **Cycles mode**: Adaptive resolution
4. **Linear vs log spacing**: Different frequency sampling
5. **Baseline correction**: Event-related power changes
6. **Edge filtering**: Removing edge artifacts

### Advantages

- **Simple and fast**: Straightforward implementation
- **Widely understood**: Classic method with extensive literature  
- **Flexible**: Can use fixed or adaptive windows
- **Predictable**: Fixed resolution makes interpretation easier

### Trade-offs

**Fixed resolution**:
- Unlike wavelets which adapt, STFT in window_length mode has same resolution everywhere
- May be suboptimal for wide frequency ranges

**Less flexible than wavelets**:
- Morlet wavelets offer more precise control over time-frequency trade-off
- But STFT is simpler and faster

### Practical Tips

**Window length**:
- **200-300 ms**: Good time precision
- **300-500 ms**: Balanced resolution
- **500-1000 ms**: Good frequency precision

**Use cycles mode for wide ranges**:
- Bridges gap between STFT and Morlet wavelets
- Gets adaptive resolution without additional complexity


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
