# TF Morlet

Time-frequency analysis using Morlet wavelets with various cycle configurations.

## Overview

Demonstrates Time-frequency analysis using Morlet wavelets with various cycle configurations.

## Source Code

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
