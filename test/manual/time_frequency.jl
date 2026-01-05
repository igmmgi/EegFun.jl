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
    [5.0, 5.0, 4.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time windows for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = eegfun.signal_to_data(times, signal, :Channel1, sample_rate)
# eegfun.plot_erp(epochs_synthetic, channel_selection = eegfun.channels([:Channel1]))
spectrum = eegfun.freq_spectrum(epochs_synthetic, max_freq=80.0)
eegfun.plot_freq_spectrum(spectrum, channel_selection = eegfun.channels([:Channel1]))

# Generate synthetic signal with noise
sample_rate = 256.0
times, signal = eegfun.generate_signal(
    10,                                      # n_trials
    [-1.0, 2.0],                            # time_window
    sample_rate,                            # sample_rate
    [2.0, 15, 25.0],                        # frequencies
    [2.0, 3.0, 2.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time windows for each freq 
    0.5,                                    # noise amplitude
);
epochs_synthetic = eegfun.signal_to_data(times, signal, :Channel1, sample_rate)
eegfun.plot_epochs(epochs_synthetic, channel_selection = eegfun.channels([:Channel1]))


tf_data = eegfun.tf_stft_fixed(epochs_synthetic, lin_freqs = (1, 40, 0.5), window_length = 0.5)

@btime tf_data = eegfun.tf_morlet(epochs_synthetic, lin_freqs = (1, 40, 1))
@btime tf_data = eegfun.tf_stft_fixed(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5)
@btime tf_data = eegfun.tf_stft_adaptive(epochs_synthetic, lin_freqs = (1, 40, 1), cycles = 5)

@btime tf_data = eegfun.tf_multitaper(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5)


spectrum = eegfun.freq_spectrum(epochs_synthetic, max_freq=100.0)
eegfun.plot_freq_spectrum(spectrum, channel_selection = eegfun.channels([:Channel1]))


tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 80, 1), window_length = 0.5) 
@btime tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5) 


# time-frequency analysis
# 13.11
@btime tf_data = eegfun.tf_morlet(epochs_synthetic, lin_freqs = (1, 80, 1))
@btime tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 80, 1), window_length = 0.5) 

@btime tf_data = eegfun.tf_multitaper(epochs_synthetic, lin_freqs = (1, 80, 1), window_length = 0.5) 

fig1 = eegfun.plot_time_frequency(tf_data)

@btime tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (2, 80, 4), window_length = 0.5) 
@btime tf_data = eegfun.tf_multitaper(epochs_synthetic, lin_freqs = (2, 80, 2), window_length = 0.5) 
fig1 = eegfun.plot_time_frequency(tf_data)

@btime tf_data = eegfun.tf_multitaper(epochs_synthetic, lin_freqs = (2, 80, 2), window_length = 0.5) 
fig1 = eegfun.plot_time_frequency(tf_data)


# Load real data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
epochs_real = eegfun.load_data(epoch_file)[1] # take single epoch
@btime tf_data = eegfun.tf_stft_fixed(epochs_real, lin_freqs = (1, 80, 1), window_length = 0.5)
@btime tf_data = eegfun.tf_stft_adaptive(epochs_real, lin_freqs = (1, 40, 1), cycles = 5)



@btime tf_data = eegfun.tf_morlet(epochs_real, lin_freqs = (1, 40, 1), time_steps = (-0.5, 2.0, 0.01)) 
@btime tf_data = eegfun.tf_stft(epochs_real, lin_freqs = (1, 40, 1), time_steps = (-0.5, 2.0, 0.01), window_length = 0.5) 
@btime tf_data = eegfun.tf_stft(epochs_real, lin_freqs = (1, 40, 1), time_steps = (-0.5, 2.0, 0.01), cycles = 5) 


@btime tf_data = eegfun.tf_morlet(epochs_real, lin_freqs = (1, 40, 1), time_steps = (-0.5, 2.0, 0.1)) 
tf_data = eegfun.tf_morlet(epochs_real, lin_freqs = (1, 40, 5), time_steps = (-0.5, 2.0, 0.01)) 
fig1 = eegfun.plot_time_frequency( tf_data, baseline_window = (-0.5, -0.2), baseline_method = :db)

@btime tf_data = eegfun.tf_stft(epochs_real, lin_freqs = (1, 40, 1), window_length = 0.5) 

@btime tf_data = eegfun.tf_morlet(epochs_real, lin_freqs = (1, 40, 1), time_steps = (-0.5, 2.0, 0.1)) 
tf_data = eegfun.tf_morlet(epochs_real, lin_freqs = (1, 40, 5), time_steps = (-0.5, 2.0, 0.01)) 
fig1 = eegfun.plot_time_frequency( tf_data, baseline_window = (-0.5, -0.2), baseline_method = :db)










# Load real EEG data from Cohen's dataset
data_cohen = eegfun.load_data("/home/ian/Desktop/tf_test_epochs.jld2")
# tf_data = eegfun.tf_morlet(data_cohen, lin_freqs = (1, 40, 0.5), cycles = 5)
# tf_data = eegfun.tf_stft_fixed(data_cohen, lin_freqs = (1, 40, 0.5), window_length = 0.5)
# tf_data = eegfun.tf_stft_adaptive(data_cohen, lin_freqs = (1, 40, 1), cycles = 5)
tf_data = eegfun.tf_multitaper(data_cohen, lin_freqs = (1, 40, 0.5), cycles = 5, time_steps = (-0.5, 2.0, 0.02))
fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)



# fig1 = eegfun.plot_time_frequency(tf_data)



# 13.11
tf_data = eegfun.tf_morlet(data_cohen, log_freqs = (2, 80, 100), cycles = 7, time_steps = (-0.5, 1.0, 0.01))
fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

# 13.11
@btime tf_data = eegfun.tf_morlet(data_cohen, lin_freqs = (2, 80, 1), cycles = 7, time_steps = (-2.5, 3.0, 0.01))

fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)


fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = false,
    colormap = :jet,
)


# Load real data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
epochs_real = eegfun.load_data(epoch_file)[1] # take single epoch

@time tf_data = eegfun.tf_morlet(epochs_real, log_freqs = (2, 80, 100), cycles = 7, time_steps = (-0.5, 1.0, 0.01))
fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)

@time tf_data =
    eegfun.tf_morlet_optim(epochs_real, log_freqs = (2, 80, 100), cycles = 7, time_steps = (-0.5, 1.0, 0.01))
fig1, _ = eegfun.plot_time_frequency(
    tf_data;
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
    colorrange = (-3, 3),
    ylogscale = true,
    colormap = :jet,
)




# Time and frequency parameters matching Cohen's Figure 13.14 exactly
# Note: Fine time resolution (0.001s) and frequency resolution (0.5 Hz) for high-quality plot
#time_of_interest_cohen = -1:(1/data_cohen.sample_rate):1.0  # 0.001s steps (matches Cohen Figure 13.14)
time_of_interest_cohen = -1:0.01:1.0  # 0.001s steps (matches Cohen Figure 13.14)
frequency_of_interest_cohen = 2:1:80   # 0.5 Hz steps (matches Cohen Figure 13.14)

@info "Testing all time-frequency methods on Cohen's real EEG data"
@info "Data: $(eegfun.filename(data_cohen)), $(eegfun.n_epochs(data_cohen)) epochs, $(eegfun.sample_rate(data_cohen)) Hz"
@info "Channels: $(eegfun.channel_labels(data_cohen))"

# 1. Wavelet method (Morlet wavelets - Cohen Ch. 13-14)
# Parameters matching Cohen's Figure 13.14:
# - width=(3, 10): 3 cycles at lowest freq, 10 cycles at highest (logarithmic scaling)
# - baseline_window=(-0.5, -0.2): Pre-stimulus baseline
# - baseline_method=:db: Decibel conversion
# - colorrange=(-3, 3): dB range matching Cohen's figure
@info "Running wavelet analysis (this may take a while due to fine time/freq resolution)..."
tf_data1_cohen = eegfun.tf_analysis(
    data_cohen,
    frequency_of_interest_cohen,
    time_of_interest_cohen;
    method = :wavelet,
    width = 3,  # Frequency-dependent cycles: 3 at low freq, 10 at high freq
    log_freqs = false,  # Use linear frequency spacing (1:2:80) instead of logarithmic
)

fig1_cohen, _ = eegfun.plot_time_frequency(
    tf_data1_cohen,
    :Channel1;
    title = "Figure 13.14 - Time-Frequency Analysis (Cohen)",
    #baseline_window = (-0.5, -0.2),
    #baseline_method = :db,
    ylogscale = false,  # Logarithmic frequency scale (matches Cohen)
    #colorrange = (-3, 3),  # dB range matching Cohen's figure
    colormap = :jet,
)























# Time and frequency parameters for TF analysis
time_of_interest = -1:0.01:3
frequency_of_interest = 1:1:40

# 1. Wavelet method
tf_data = tf_analysis(epochs_synthetic, frequency_of_interest, time_of_interest; method = :wavelet, width = 7)
fig1, _ = eegfun.plot_time_frequency(tf_data, :Channel1; title = "tf_analysis (wavelet, 7 cycles)")



# Load real EEG data from Cohen's dataset
data_cohen = eegfun.load_data("/home/ian/Desktop/tf_test_epochs.jld2")

# Time and frequency parameters matching Cohen's Figure 13.14 exactly
# Note: Fine time resolution (0.001s) and frequency resolution (0.5 Hz) for high-quality plot
#time_of_interest_cohen = -1:(1/data_cohen.sample_rate):1.0  # 0.001s steps (matches Cohen Figure 13.14)
time_of_interest_cohen = -1:0.01:1.0  # 0.001s steps (matches Cohen Figure 13.14)
frequency_of_interest_cohen = 2:1:80   # 0.5 Hz steps (matches Cohen Figure 13.14)

@info "Testing all time-frequency methods on Cohen's real EEG data"
@info "Data: $(eegfun.filename(data_cohen)), $(eegfun.n_epochs(data_cohen)) epochs, $(eegfun.sample_rate(data_cohen)) Hz"
@info "Channels: $(eegfun.channel_labels(data_cohen))"

# 1. Wavelet method (Morlet wavelets - Cohen Ch. 13-14)
# Parameters matching Cohen's Figure 13.14:
# - width=(3, 10): 3 cycles at lowest freq, 10 cycles at highest (logarithmic scaling)
# - baseline_window=(-0.5, -0.2): Pre-stimulus baseline
# - baseline_method=:db: Decibel conversion
# - colorrange=(-3, 3): dB range matching Cohen's figure
@info "Running wavelet analysis (this may take a while due to fine time/freq resolution)..."
tf_data1_cohen = eegfun.tf_analysis(
    data_cohen,
    frequency_of_interest_cohen,
    time_of_interest_cohen;
    method = :wavelet,
    width = 3,  # Frequency-dependent cycles: 3 at low freq, 10 at high freq
    log_freqs = false,  # Use linear frequency spacing (1:2:80) instead of logarithmic
)

fig1_cohen, _ = eegfun.plot_time_frequency(
    tf_data1_cohen,
    :Channel1;
    title = "Figure 13.14 - Time-Frequency Analysis (Cohen)",
    #baseline_window = (-0.5, -0.2),
    #baseline_method = :db,
    ylogscale = false,  # Logarithmic frequency scale (matches Cohen)
    #colorrange = (-3, 3),  # dB range matching Cohen's figure
    colormap = :jet,
)



# Load real data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
epochs_real = eegfun.load_data(epoch_file)[1] # take single epoch

tf_data = eegfun.tf_analysis(epochs_real, frequency_of_interest, time_of_interest; method = :wavelet, width = 7)
fig1, _ = eegfun.plot_time_frequency(
    tf_data,
    :Cz;
    title = "tf_analysis (wavelet, 7 cycles)",
    ylogscale = false,
    baseline_window = (-0.5, -0.2),
    baseline_method = :db,
)

# 2. Hanning with fixed window length
tf_data2 = eegfun.tf_analysis(
    epochs,
    frequency_of_interest,
    time_of_interest;
    method = :hanning_fixed,
    time_window_length = 0.3,
)
fig2, _ = eegfun.plot_time_frequency(tf_data2, :Channel1; title = "tf_analysis (hanning_fixed, 0.3s window)")

# 3. Hanning with frequency-dependent window length
tf_data2b = eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :hanning_adaptive, cycles = 5)
fig2b, _ = eegfun.plot_time_frequency(tf_data2b, :Channel1; title = "tf_analysis (hanning_adaptive, 5 cycles)")

# 4. Multitaper with DPSS  
tf_data3 = eegfun.tf_analysis(
    epochs,
    frequency_of_interest,
    time_of_interest;
    method = :multitaper,
    time_window_length = 0.3,
    frequency_smoothing = 4.0,
)
fig3, _ = eegfun.plot_time_frequency(tf_data3, :Channel1; title = "tf_analysis (multitaper, DPSS)")

# 5. Superlet Unified method
tf_data4 = eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :superlet, order = 5)
fig4, _ = eegfun.plot_time_frequency(tf_data4, :Channel1; title = "tf_analysis (superlet, 5th order)")

# 6. Power spectrum (using tf_spectrum on raw signal for comparison)
power_spectrum, freqs_out =
    eegfun.tf_spectrum(signal, sample_rate, frequency_of_interest; taper = :dpss, frequency_smoothing = 2.0)
fig5 = Figure(size = (600, 400))
ax = Axis(fig5[1, 1], xlabel = "Frequency (Hz)", ylabel = "Power", title = "tf_spectrum (peaks at 2Hz, 15Hz, 25Hz)")
lines!(ax, freqs_out, power_spectrum)
display(fig5)


#######################################################################
@info eegfun.section("TEST 2: Real EEG Data")
#######################################################################

# Try to load epoched data from output_data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
epochs = eegfun.load_data(epoch_file)[1] # take single epoch

# Time and frequency parameters
time_of_interest = -0.5:0.01:2
frequency_of_interest = 1:1:40

# 1. Wavelet method
@time tf_data1 = eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :wavelet, width = 7)

# 2. Hanning with fixed window length
@time tf_data2 = eegfun.tf_analysis(
    epochs,
    frequency_of_interest,
    time_of_interest;
    method = :hanning_fixed,
    time_window_length = 0.3,
)

# 3. Hanning with frequency-dependent window length
@btime tf_data2b =
    eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :hanning_adaptive, cycles = 5)

# 4. Multitaper with DPSS  
@btime tf_data3 = eegfun.tf_analysis(
    epochs,
    frequency_of_interest,
    time_of_interest;
    method = :multitaper,
    time_window_length = 0.3,
    frequency_smoothing = 4.0,
)

# 5. Superlet Unified method
tf_data4 = eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :superlet, order = 5)




#######################################################################
@info eegfun.section("TEST 3: Baseline Correction Examples")
#######################################################################

# Baseline correction removes the mean power during a baseline period
# This is essential for visualizing event-related changes in power

selected_channel = :Cz
# Create a TF analysis for baseline examples
@btime tf_baseline_example = eegfun.tf_analysis(
    epochs,
    frequency_of_interest,
    time_of_interest;
    method = :multitaper,
    width = 7,
    channel_selection = eegfun.channels(selected_channel),
)

# Create a TF analysis for baseline examples
tf_baseline_example = eegfun.tf_analysis(epochs, frequency_of_interest, time_of_interest; method = :wavelet, width = 7)


# Example 1: Apply baseline correction using tf_baseline function
# Baseline window: -0.3 to 0.0 seconds (pre-stimulus period)
baseline_window = (-0.3, 0.0)
# Method 1: Decibel (dB) - most common, shows power relative to baseline
# Formula: 10 * log10(power / baseline_mean)
tf_db = eegfun.tf_baseline(tf_baseline_example, baseline_window; method = :db)
fig_b1_no_baseline, _ = eegfun.plot_time_frequency(
    tf_db,
    selected_channel;
    title = "Baseline: dB ($(baseline_window[1]) to $(baseline_window[2])s)",
)

fig_b1_with_baseline, _ = eegfun.plot_time_frequency(
    tf_db,
    selected_channel;
    title = "Baseline: dB ($(baseline_window[1]) to $(baseline_window[2])s)",
)

# Method 2: Percent change
# Formula: 100 * (power - baseline_mean) / baseline_mean
tf_percent = eegfun.tf_baseline(tf_baseline_example, baseline_window; method = :percent)
fig_b2, _ = eegfun.plot_time_frequency(
    tf_percent,
    selected_channel;
    title = "Baseline: Percent Change ($(baseline_window[1]) to $(baseline_window[2])s)",
)

# Method 3: Relative change (ratio)
# Formula: power / baseline_mean
tf_rel = eegfun.tf_baseline(tf_baseline_example, baseline_window; method = :relchange)
fig_b3, _ = eegfun.plot_time_frequency(
    tf_rel,
    selected_channel;
    title = "Baseline: Relative Change ($(baseline_window[1]) to $(baseline_window[2])s)",
)

# Check values look reasonable
println("\n" * "="^60)
println("Baseline Correction Value Check:")
println("="^60)
times_unique = sort(unique(tf_baseline_example.data.time))
freqs_unique = sort(unique(tf_baseline_example.data.freq))
# Get a sample time point (e.g., 0.1s after baseline ends)
sample_time_idx = findfirst(t -> t >= 0.1, times_unique)
sample_freq_idx = findfirst(f -> f >= 10, freqs_unique)  # Around 10Hz

if !isnothing(sample_time_idx) && !isnothing(sample_freq_idx)
    sample_time = times_unique[sample_time_idx]
    sample_freq = freqs_unique[sample_freq_idx]

    # Extract values at this time-frequency point
    row_idx = findfirst(row -> row.time ≈ sample_time && row.freq ≈ sample_freq, eachrow(tf_baseline_example.data))
    if !isnothing(row_idx)
        original_val = tf_baseline_example.data[row_idx, selected_channel]
        db_val = tf_db.data[row_idx, selected_channel]
        percent_val = tf_percent.data[row_idx, selected_channel]
        rel_val = tf_rel.data[row_idx, selected_channel]

        println("Sample point: $(round(sample_time, digits=2))s, $(round(sample_freq, digits=1))Hz")
        println("Original power: $(round(original_val, digits=4))")
        println("dB: $(round(db_val, digits=2)) dB")
        println("Percent: $(round(percent_val, digits=2))%")
        println("Relative: $(round(rel_val, digits=3))x")
        println("\nExpected relationships:")
        println("  dB ≈ 10 * log10(relative) = $(round(10 * log10(rel_val), digits=2))")
        println("  Percent ≈ 100 * (relative - 1) = $(round(100 * (rel_val - 1), digits=2))")
    end
end

# Example 2: Apply baseline directly in plotting (convenience method)
# This applies baseline correction on-the-fly without modifying the original data
fig_b5, _ = eegfun.plot_time_frequency(
    tf_baseline_example,
    selected_channel;
    baseline_window = baseline_window,
    baseline_method = :percent,
    title = "Plot with baseline_window parameter (dB)",
)


#######################################################################
@info eegfun.section("TEST 4: Performance Benchmarks")
#######################################################################

println("="^60)
println("Time-Frequency Performance Benchmark")
println("="^60)
println()

# Generate synthetic test data for benchmarking
println("Generating synthetic test data...")
sample_rate = 256.0
times, signal = eegfun.generate_signal(
    100,                                    # n_trials (typical for EEG)
    [-0.5, 2.5],                           # time_window
    sample_rate,                            # sample_rate
    [10.0, 20.0],                          # frequencies
    [1.0, 1.0],                            # amplitudes
    [[0.0, 1.0], [1.0, 2.0]],             # time windows for each freq
    0.5,                                     # noise amplitude
)

# Convert to EpochData
epochs_bench = eegfun._signal_to_epochs(times, signal, :TestChannel, Int(sample_rate))

# Time and frequency parameters
toi = -0.5:0.01:2.0
foi = 1:1:40

println("Benchmark configuration:")
println("  Trials: $(length(epochs_bench.data))")
println("  Samples per trial: $(nrow(epochs_bench.data[1]))")
println("  Time points: $(length(toi))")
println("  Frequencies: $(length(foi))")
println()

# Warmup
println("Warming up...")
eegfun.tf_analysis(epochs_bench, foi, toi; method = :multitaper, time_window_length = 0.5, frequency_smoothing = 4.0)

println("\n" * "="^60)
println("BENCHMARK: Multitaper (DPSS) Method")
println("="^60)
println()

# Run benchmark
println("Running benchmark (this may take a few seconds)...")
@time tf_data_bench = eegfun.tf_analysis(
    epochs_bench,
    foi,
    toi;
    method = :multitaper,
    time_window_length = 0.5,
    frequency_smoothing = 4.0,
)

println("\nDetailed benchmark:")
@btime eegfun.tf_analysis(
    $epochs_bench,
    $foi,
    $toi;
    method = :multitaper,
    time_window_length = 0.5,
    frequency_smoothing = 4.0,
)

println("\n" * "="^60)
println("Results:")
println("="^60)
println("Output size: $(size(tf_data_bench.data))")
println("Channels: $(names(tf_data_bench.data)[3:end])")
println()
println("="^60)
