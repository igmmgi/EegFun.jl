using EegFun
using GLMakie
using JLD2
using DataFrames
using BenchmarkTools

#######################################################################
@info EegFun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal
sample_rate = 1000.0
times, signal = EegFun.generate_signal(
    10,                                      # n_trials
    [-1.0, 3.0],                            # time_window
    sample_rate,                            # sample_rate
    [5.0, 25, 35.0],                        # frequencies
    [5.0, 5.0, 5.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time windows for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = EegFun.signal_to_data(times, signal, :Channel1, sample_rate)
EegFun.plot_epochs(epochs_synthetic, channel_selection = EegFun.channels([:Channel1]))

spectrum = EegFun.freq_spectrum(epochs_synthetic, max_freq=80.0)
EegFun.plot_freq_spectrum(spectrum, channel_selection = EegFun.channels([:Channel1]))

# tf_stft_fixed
tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = logrange(1, 40, length=30), cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

tf_data = EegFun.tf_multitaper(epochs_synthetic, frequencies = logrange(1, 40, length=30), cycles = 5)
EegFun.plot_time_frequency(tf_data, ylogscale = true)

#######################################################################
@info EegFun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################

data_cohen = EegFun.load_data("/home/ian/Desktop/tf_test_epochs.jld2")

# Figure 13.11 A)
tf_data = EegFun.tf_multitaper(data_cohen, frequencies = logrange(2, 80, length=30), cycles = 5, time_steps = 0.05, filter_edges = true)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)



#######################################################################
@info EegFun.section("TEST 3: Exported data from MATLAB FieldTrip")
#######################################################################

epochs = EegFun.load_csv("/home/ian/Documents/MATLAB/BioPsyLab/Data/TestData/data1/", file = "epoch_data.csv")


tf_data = EegFun.tf_multitaper(
    epochs;
    channel_selection = EegFun.channels(:Cz),
    frequencies = 1:2:30,           # cfg.foi = 1:2:30
    cycles = 7,                       # cfg.t_ftimwin = 3./cfg.foi
    frequency_smoothing = 0.3,        # cfg.tapsmofrq = 0.1 * cfg.foi
    time_steps = 0.05    # cfg.toi = -0.5:0.05:1.5
)
EegFun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.1),
    baseline_method = :relchange,
    # colorrange = (-0.6, 1),
    ylogscale = false,
    #colormap = :jet,
)







