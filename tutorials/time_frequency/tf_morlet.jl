# Demo: Time-Frequency Analysis - Morlet Wavelets
# Shows wavelet-based time-frequency decomposition with synthetic signals
# and real data, demonstrating different cycle counts and frequency resolution.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

#######################################################################
@info EegFun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal
sample_rate = 1000.0
times, signal = EegFun.generate_signal(
    400,                                    # n_trials
    [-1.0, 3.0],                            # time_interval
    sample_rate,                            # sample_rate
    [5.0, 25, 35.0],                        # frequencies
    [5.0, 5.0, 5.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time intervals for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = EegFun.signal_to_data(times, signal, :Channel1, sample_rate)
# EegFun.plot_epochs(epochs_synthetic, channel_selection = EegFun.channels([:Channel1]))

spectrum = EegFun.freq_spectrum(epochs_synthetic, max_freq = 80.0)
# EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels([:Channel1]))

# tf_morlet
@time tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, pad = :both, use_gpu = true)
EegFun.plot_tf(tf_data, ylogscale = false)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:0.5:40, cycles = 10)
EegFun.plot_tf(tf_data, ylogscale = false)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 3)
EegFun.plot_tf(tf_data, ylogscale = true)

tf_data = EegFun.tf_morlet(epochs_synthetic, frequencies = logrange(1, 40, length = 30), cycles = 10)
EegFun.plot_tf(tf_data, ylogscale = true)


#######################################################################
@info EegFun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################
# This is some data that was presented in Cohen: Analyzin Neural Time Series Data
data_cohen = EegFun.read_data(EegFun.example_path("data/julia/tf/tf_test_epochs.jld2"));

# Figure 13.11 A)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 80), cycles = 3, filter_edges = false)
EegFun.plot_tf(
    tf_data;
    title = "3-cycle wavelets",
    colorbar_label = "dB change from baseline",
    baseline_interval = (-0.5, -0.2),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
    interpolate = true,
    yticks = [2, 4, 10, 20, 40, 80],
    xticks = [0, 0.2, 0.4, 0.6, 0.8],
    time_unit = :ms,
)

# Figure 13.11 B)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = 2:1:80, cycles = (3, 10))
EegFun.plot_tf(
    tf_data;
    baseline_interval = (-0.5, -0.2),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)

# Figure 13.14 A)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 40), cycles = 3)
EegFun.plot_tf(
    tf_data;
    baseline_interval = (-0.5, -0.2),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 B)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 40), cycles = 10, filter_edges = false)
EegFun.plot_tf(
    tf_data;
    baseline_interval = (-0.5, -0.2),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 C)
tf_data = EegFun.tf_morlet(data_cohen, frequencies = logrange(2, 80, length = 40), cycles = (3, 10))
EegFun.plot_tf(
    tf_data;
    baseline_interval = (-0.5, -0.2),
    baseline_method = :db,
    xlim = (-0.2, 1.0),
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
    interpolate = true,
)
