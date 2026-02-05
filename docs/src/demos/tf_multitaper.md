# TF Multitaper

## Overview

## Overview

This demo demonstrates time-frequency analysis using the multitaper method.

### Multitaper Method

Reduce variance in spectral estimates using multiple tapers:
- **Tapers**: Data windows with optimal properties (Slepian sequences)
- **Multiple estimates**: Average across tapers for reduced variance
- **Time-bandwidth product**: Controls frequency smoothing

### Key Parameters

**Time-bandwidth product (NW):**
- Controls spectral concentration
- Typical values: 2-4
- Higher = more frequency smoothing, less variance

**Number of tapers:**
- Usually 2×NW - 1
- More tapers = smoother spectra

### Advantages

- Superior variance reduction vs. single-taper methods
- Minimal spectral leakage
- Optimal for narrow-band analyses
- Gold standard for power spectral estimation

### Applications

- High-quality power spectra
- Coherence analysis
- Stationary oscillatory activity
- Resting-state analyses


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
