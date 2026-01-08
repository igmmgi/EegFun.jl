using eegfun
using GLMakie
using JLD2
using DataFrames
using BenchmarkTools

#######################################################################
@info eegfun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal
sample_rate = 1000.0
times, signal = eegfun.generate_signal(
    100,                                      # n_trials
    [-1.0, 3.0],                            # time_window
    sample_rate,                            # sample_rate
    [5.0, 25, 35.0],                        # frequencies
    [5.0, 5.0, 5.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time windows for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = eegfun.signal_to_data(times, signal, :Channel1, sample_rate)
eegfun.plot_epochs(epochs_synthetic, channel_selection = eegfun.channels([:Channel1]))

spectrum = eegfun.freq_spectrum(epochs_synthetic, max_freq=80.0)
eegfun.plot_freq_spectrum(spectrum, channel_selection = eegfun.channels([:Channel1]))

# tf_stft_fixed
tf_data = eegfun.tf_stft_fixed(epochs_synthetic, lin_freqs = (1, 40, 0.5), window_length = 0.5)
eegfun.plot_time_frequency(tf_data, ylogscale = false)

tf_data = eegfun.tf_stft_fixed(epochs_synthetic, lin_freqs = (1, 40, 0.5), window_length = 0.5)
eegfun.plot_time_frequency(tf_data, log_freqs = false)

tf_data = eegfun.tf_stft_fixed(epochs_synthetic, log_freqs = (1, 40, 30), window_length = 0.5)
eegfun.plot_time_frequency(tf_data, ylogscale = true)

tf_data = eegfun.tf_stft_fixed(epochs_synthetic, log_freqs = (1, 40, 30), window_length = 0.5)
eegfun.plot_time_frequency(tf_data, ylogscale = true)

#######################################################################
@info eegfun.section("TEST 2: Cohen Data Chapter 13")
#######################################################################

data_cohen = eegfun.load_data("/home/ian/Desktop/tf_test_epochs.jld2")

# Figure 13.11 A)
tf_data = eegfun.tf_stft_fixed(data_cohen, log_freqs = (2, 80, 30), window_length = 0.5, time_steps = (-0.5, 1.0, 0.005), filter_edges = false)
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.11 A)
tf_data = eegfun.tf_morlet(data_cohen, lin_freqs = (2, 80, 1), cycles = (3, 10), time_steps = (-0.5, 1.0, 0.005))
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)

# Figure 13.14 A)
tf_data = eegfun.tf_morlet(data_cohen, log_freqs = (2, 80, 30), cycles = 3, time_steps = (-0.5, 1.0, 0.005))
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 B)
tf_data = eegfun.tf_morlet(data_cohen, log_freqs = (2, 80, 30), cycles = 10, time_steps = (-0.5, 1.0, 0.005))
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# Figure 13.14 C)
tf_data = eegfun.tf_morlet(data_cohen, log_freqs = (2, 80, 30), cycles = (3, 10), time_steps = (-0.5, 1.0, 0.005))
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

#######################################################################
@info eegfun.section("TEST 3: Exported data from MATLAB FieldTrip")
#######################################################################

data = eegfun.load_csv("/home/ian/Documents/MATLAB/BioPsyLab/Data/TestData/data1/", file = "epoch_data.csv")

tf_data = eegfun.tf_stft(data, channel_selection = eegfun.channels([:Cz]), lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01)) 
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.1),
    baseline_method = :relchange,
    # colorrange = (-0.6, 1),
    ylogscale = false,
    #colormap = :jet,
)

tf_data = eegfun.tf_stft(data, channel_selection = eegfun.channels([:Cz]), lin_freqs = (2, 30, 1), cycles = 7, time_steps = (-0.5, 1.5, 0.05)) 
eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.1),
    baseline_method = :relchange,
    # colorrange = (-0.6, 1),
    ylogscale = false,
    #colormap = :jet,
)
