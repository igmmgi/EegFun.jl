# TF Morlet

This demo demonstrates time-frequency analysis using Morlet wavelets to decompose EEG signals into time-varying frequency content.

This demo demonstrates time-frequency analysis using Morlet wavelets to decompose EEG signals into time-varying frequency content.

### What are Morlet Wavelets?

Morlet wavelets are oscillating waves with Gaussian envelopes that provide localized time-frequency analysis:

- **Wavelet**: Small wave that is both time-limited and frequency-specific
- **Gaussian envelope**: Smooth windowing that minimizes spectral leakage
- **Time-frequency trade-off**: Cannot have perfect precision in both simultaneously

### Key Parameter: Number of Cycles

The `cycles` parameter controls the time-frequency resolution trade-off:

**Fixed cycles** (e.g., `cycles = 5`):

- Same spectral bandwidth at all frequencies
- **Low frequencies**: Broad temporal windows (poor time resolution)
- **High frequencies**: Narrow temporal windows (good time resolution)

**Variable cycles** (e.g., `cycles = (3, 10)`):

- Adaptive resolution across frequencies
- **Low frequencies**: Fewer cycles = better temporal precision
- **High frequencies**: More cycles = better spectral precision
- Balances resolution optimally across the spectrum

### Frequency Spacing

**Linear spacing** (`frequencies = 1:1:40`):

- Equal Hz steps
- Good for narrow frequency ranges
- Easier interpretation

**Logarithmic spacing** (`frequencies = logrange(1, 40, length = 30)`):

- Proportional frequency steps
- Better for wide ranges (e.g., 2-80 Hz)
- Use with `ylogscale = true` in plots

### Baseline Correction

Apply baseline correction to isolate event-related changes:

```julia
baseline_window = (-0.5, -0.2)  # Pre-stimulus baseline
baseline_method = :db           # Decibel conversion (10*log10(activity/baseline))
```

**Baseline methods**:

- **:db**: Decibel scale (relative power change)
- **:percent**: Percent change from baseline
- **:zscore**: Z-score normalization

### Edge Effects

The `filter_edges` and `pad` parameters control edge artifacts:

```julia
tf_morlet(epochs, pad = :both)         # Mirror padding
tf_morlet(epochs, filter_edges = true) # Remove edge samples
```

### Workflow Summary

This demo shows:

1. **Synthetic validation**: Known frequencies (5, 25, 35 Hz) for sanity checking
2. **Fixed cycles**: Constant resolution across frequencies
3. **Variable cycles**: Adaptive resolution with `cycles = (3, 10)`
4. **Linear vs log spacing**: Different frequency sampling strategies
5. **Baseline correction**: Isolating event-related power changes


### Further Reading

Cohen, M. X. (2014). *Analyzing Neural Time Series Data: Theory and Practice*. Chapter 12: Morlet Wavelets and Wavelet Convolution


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

# tf_morlet
tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, pad = :both)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

# tf_morlet
tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:0.5:40, cycles = 10)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 3)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 10)
EegFun.plot_time_frequency(tf_data, ylogscale = true)


#######################################################################
@info EegFun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################
# This is some data that was presented in Cohen: Analyzin Neural Time Series Data
data_cohen = EegFun.read_data("./data/files/tf_test_epochs.jld2");

# Figure 13.11 A)
tf_data = EegFun.tf_morlet(reconstructed_data, frequencies = range(2, 80, length = 80), cycles = 3, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)

# Figure 13.11 A)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = 2:1:80, cycles = (3, 10))
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)

# Figure 13.14 A)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 30), cycles = 3)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 B)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 30), cycles = 10)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 C)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 30), cycles = (3, 10))
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
