# TF Multitaper

This demo demonstrates time-frequency analysis using the multitaper method for high-quality spectral estimation.

This demo demonstrates time-frequency analysis using the multitaper method for high-quality spectral estimation.

### What is the Multitaper Method?

The multitaper method improves spectral estimates by averaging across multiple orthogonal tapers (windows):

- **Slepian tapers**: Mathematically optimal windows that minimize spectral leakage
- **Multiple estimates**: Each taper provides an independent spectral estimate
- **Averaging**: Reduces variance while preserving spectral features
- **Time-bandwidth product (NW)**: Controls frequency smoothing

### Key Parameters

**Cycles**:

Similar to Morlet wavelets, controls the time window length:
```julia
tf_multitaper(epochs, cycles = 5)
```
Higher cycles = longer windows = better frequency resolution

**Time-bandwidth product**:

Internally determined by the number of cycles and controls:
- Frequency smoothing bandwidth
- Number of tapers used
- Variance reduction

### Advantages Over Single-Taper Methods

**Variance reduction**:
- Multiple independent estimates average out noise
- Smoother, more stable spectra
- Better for noisy data

**Minimal spectral leakage**:
- Slepian tapers have optimal concentration properties
- Less contamination from nearby frequencies

**Gold standard**:
- Considered best practice for power spectral estimation
- Widely used in neuroscience

### Frequency Spacing

**Linear spacing** (`frequencies = 1:1:40`):
- Equal Hz steps across range
- Good for narrow frequency bands

**Logarithmic spacing** (`frequencies = logrange(1, 40, length = 30)`):
- Proportional steps (constant ratios)
- Better for wide ranges (2-80 Hz)
- Use with `ylogscale = true` for visualization

### Baseline Correction

```julia
plot_time_frequency(
    tf_data,
    baseline_window = (-0.5, -0.2),
    baseline_method = :db
)
```

Isolates event-related changes relative to pre-stimulus baseline.

### When to Use Multitaper

**Best for**:
- High-quality power spectra
- Narrow-band oscillatory analyses
- Stationary signals
- When variance reduction is critical

**Consider alternatives**:
- Morlet wavelets: Better time precision at high frequencies
- STFT: Simpler, faster for exploratory work

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz)
2. **Fixed cycles**: Constant resolution analysis
3. **Linear vs log spacing**: Different frequency sampling strategies
4. **Baseline correction**: Event-related power changes
5. **Edge filtering**: Removing edge artifacts

### Practical Tips

**Choose cycles based on needs**:
- **3-5 cycles**: Better temporal resolution
- **7-10 cycles**: Better frequency resolution
- **5 cycles**: Good balance for most EEG work

**Frequency sampling**:
- Dense sampling (every 0.5 Hz) for detailed spectra
- Coarse sampling (every 2-5 Hz) for faster computation

**Time steps**:
- Control temporal smoothness of the output
- Smaller steps = smoother time course but slower computation


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
tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

#######################################################################
@info EegFun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################
# This is some data that was presented in Cohen: Analyzin Neural Time Series Data
data_cohen = EegFun.read_data("./data/files/tf_test_epochs.jld2");

# Figure 13.11 A)
tf_data = EegFun.tf_multitaper(data_cohen, frequencies = logrange(2, 80, length = 30), cycles = 5, time_steps = 0.05, filter_edges = true)
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
