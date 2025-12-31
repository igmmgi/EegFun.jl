using eegfun
using GLMakie
using JLD2
using DataFrames
using BenchmarkTools

#######################################################################
@info eegfun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal: 2Hz from -0.1-0.5s, 15Hz from 0.5-1.1s, 25Hz from 1.5-2.0s
sample_rate = 256.0
times, signal = eegfun.generate_signal(
    300,                                    # n_trials
    [-1.0, 3.0],                            # time_window
    sample_rate,                            # sample_rate
    [2.0, 15, 25.0],                        # frequencies
    [1.0, 2.0, 3.0],                        # amplitudes
    [[-0.1, 0.5], [0.5, 1.1], [1.5, 2.0]],  # time windows for each freq
    2.0                                     # noise amplitude
);

# Convert signal to EpochData format
epochs = eegfun.signal_to_epochs(times, signal, :Channel1, Int(sample_rate))

# Time and frequency parameters for TF analysis
toi = -1:0.01:3
foi = 1:1:40

# 1. Wavelet method
tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)
fig1, _ = eegfun.plot_time_frequency(tf_data1, :Channel1; title="tf_analysis (wavelet, 7 cycles)")

# 2. Multitaper with Hanning
tf_data2 = eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, taper=:hanning, t_ftimwin=0.3)
fig2, _ = eegfun.plot_time_frequency(tf_data2, :Channel1; title="tf_analysis (multitaper, Hanning)")

# 3. Multitaper with DPSS  
tf_data3 = eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, taper=:dpss, t_ftimwin=0.3, tapsmofrq=4.0)
fig3, _ = eegfun.plot_time_frequency(tf_data3, :Channel1; title="tf_analysis (multitaper, DPSS)")

# 4. Superlet Unified method
tf_data4 = eegfun.tf_analysis(epochs, foi, toi; method=:superlet, order=5)
fig4, _ = eegfun.plot_time_frequency(tf_data4, :Channel1; title="tf_analysis (superlet, 5th order)")

# 5. Power spectrum (using tf_spectrum on raw signal for comparison)
power_spectrum, freqs_out = eegfun.tf_spectrum(signal, sample_rate, foi; taper=:dpss, tapsmofrq=2.0)
fig5 = Figure(size=(600, 400))
ax = Axis(fig5[1, 1], xlabel="Frequency (Hz)", ylabel="Power", title="tf_spectrum (peaks at 2Hz, 15Hz, 25Hz)")
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
toi = -0.5:0.01:2
foi = 1:1:40

# TODO: Julia code seems really slow!!!  # takes approx 7 seconds vs. 2 seconds in matlab (check this!!)
tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:multitaper, width=7)



# Get data from a central electrode
selected_channel = [:Cz, :Fp1, :Fp2]

# Time and frequency parameters
toi = -0.5:0.01:2
foi = 1:1:40

# Use tf_analysis - it handles all the data extraction automatically!
# 1. Wavelet
@btime tf_data1 = eegfun.tf_analysis(epochs, foi, toi; 
                                method=:multitaper, taper=:hanning, t_ftimwin=0.3,
                                channel_selection=eegfun.channels(selected_channel))


@time tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)
tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)


fig1, _ = eegfun.plot_time_frequency(tf_data1, selected_channel; title="Wavelet - $selected_channel") 

# 2. Wavelet with variable cycles
tf_data2 = eegfun.tf_analysis(epochs, foi, toi; 
                                method=:wavelet, width=(3, 10),
                                channel_selection=eegfun.channels(selected_channel))
fig2, _ = eegfun.plot_time_frequency(tf_data2, selected_channel[1]; title="Wavelet Variable - $selected_channel")

# 3. Superlet
tf_data3 = eegfun.tf_analysis(epochs, foi, toi; 
                                method=:superlet, order=5,
                                channel_selection=eegfun.channels(selected_channel))
fig3, _ = eegfun.plot_time_frequency(tf_data3, selected_channel[1]; title="Superlet - $selected_channel")

# 4. Multitaper Hanning
tf_data4 = eegfun.tf_analysis(epochs, foi, toi; 
                                method=:multitaper, taper=:hanning, t_ftimwin=0.3,
                                channel_selection=eegfun.channels(selected_channel))
fig4, _ = eegfun.plot_time_frequency(tf_data4, selected_channel[1]; title="MTM Hanning - $selected_channel")

# 5. Multitaper DPSS
@btime tf_data5 = eegfun.tf_analysis(epochs, foi, toi; 
                                method=:multitaper, taper=:dpss, t_ftimwin=0.3, tapsmofrq=4.0,
                                channel_selection=eegfun.channels(selected_channel))
fig5, _ = eegfun.plot_time_frequency(tf_data5, selected_channel[1]; title="MTM DPSS - $selected_channel") 


#######################################################################
@info eegfun.section("TEST 3: Baseline Correction Examples")
#######################################################################

# Baseline correction removes the mean power during a baseline period
# This is essential for visualizing event-related changes in power

selected_channel = :Cz
# Create a TF analysis for baseline examples
@btime tf_baseline_example = eegfun.tf_analysis(epochs, foi, toi; 
                                         method=:multitaper, width=7,
                                         channel_selection=eegfun.channels(selected_channel))

# Create a TF analysis for baseline examples
tf_baseline_example = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)


# Example 1: Apply baseline correction using tf_baseline function
# Baseline window: -0.3 to 0.0 seconds (pre-stimulus period)
baseline_window = (-0.3, 0.0)
# Method 1: Decibel (dB) - most common, shows power relative to baseline
# Formula: 10 * log10(power / baseline_mean)
tf_db = eegfun.tf_baseline(tf_baseline_example, baseline_window; method=:db)
fig_b1_no_baseline, _ = eegfun.plot_time_frequency(tf_db, selected_channel; 
                                       title="Baseline: dB ($(baseline_window[1]) to $(baseline_window[2])s)")

fig_b1_with_baseline, _ = eegfun.plot_time_frequency(tf_db, selected_channel; 
                                       title="Baseline: dB ($(baseline_window[1]) to $(baseline_window[2])s)")

# Method 2: Percent change
# Formula: 100 * (power - baseline_mean) / baseline_mean
tf_percent = eegfun.tf_baseline(tf_baseline_example, baseline_window; method=:percent)
fig_b2, _ = eegfun.plot_time_frequency(tf_percent, selected_channel; 
                                       title="Baseline: Percent Change ($(baseline_window[1]) to $(baseline_window[2])s)")

# Method 3: Relative change (ratio)
# Formula: power / baseline_mean
tf_rel = eegfun.tf_baseline(tf_baseline_example, baseline_window; method=:relchange)
fig_b3, _ = eegfun.plot_time_frequency(tf_rel, selected_channel; 
                                       title="Baseline: Relative Change ($(baseline_window[1]) to $(baseline_window[2])s)")

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
fig_b5, _ = eegfun.plot_time_frequency(tf_baseline_example, selected_channel; 
                                       baseline_window=baseline_window,
                                       baseline_method=:percent,
                                       title="Plot with baseline_window parameter (dB)")
