using eegfun
using GLMakie
using JLD2
using DataFrames


# ============================================================================
# Helper function for plotting TF results
# ============================================================================

function plot_tf_result(power, times, freqs; title="Time-Frequency Analysis", 
                        colormap=:viridis, baseline_window=nothing)
    fig = Figure(size=(800, 500))
    
    # Apply baseline if specified
    if !isnothing(baseline_window)
        power = eegfun.apply_tf_baseline_db(power, times, baseline_window)
        colorbar_label = "Power (dB)"
    else
        colorbar_label = "Power"
    end
    
    ax = Axis(fig[1, 1], 
              xlabel="Time (s)", 
              ylabel="Frequency (Hz)",
              title=title)
    
    hm = heatmap!(ax, collect(times), collect(freqs), power', colormap=colormap)
    Colorbar(fig[1, 2], hm, label=colorbar_label)
    
    display(fig)
    return fig, ax
end




# ============================================================================
# Test 1: Synthetic Signal with Known Frequencies
# ============================================================================
println("\n" * "="^60)
println("TEST 1: Synthetic Signal with Known Frequencies")
println("="^60)
 
# Generate synthetic signal: 10Hz from 0-0.5s, 25Hz from 0.5-1s
sample_rate = 256.0
times, signal = eegfun.generate_signal(
    300,                                    # n_trials
    [-1.0, 3.0],                            # time_window
    sample_rate,                            # sample_rate
    [2.0, 15, 25.0],                        # frequencies
    [1.0, 2.0, 3.0],                       # amplitudes
    [[-0.1, 0.5], [0.5, 1.1], [1.5, 2.0]],  # time windows for each freq
    2.0                                     # noise amplitude
)

# Time and frequency parameters
toi = -1:0.01:3
foi = 1:1:40
  
# 1. Wavelet method
power, times_out, freqs_out = eegfun.tf_wavelet(signal, times, sample_rate, foi, toi; width=7)
fig1, _ = plot_tf_result(power, times_out, freqs_out; title="tf_wavelet (7 cycles)")

# 2. Multitaper with Hanning
power, times_out, freqs_out = eegfun.tf_multitaper(signal, times, sample_rate, foi, toi; taper=:hanning, t_ftimwin=0.3)
fig2, _ = plot_tf_result(power, times_out, freqs_out; title="tf_multitaper (Hanning)")

# 3. Multitaper with DPSS  
power, times_out, freqs_out = eegfun.tf_multitaper(signal, times, sample_rate, foi, toi; taper=:dpss, t_ftimwin=0.3, tapsmofrq=4.0)
fig3, _ = plot_tf_result(power, times_out, freqs_out; title="tf_multitaper (DPSS)")

# 4. Superlet Unified method
power, times_out, freqs_out = eegfun.tf_superlet(signal, times, sample_rate, foi, toi; order=5)
fig4, _ = plot_tf_result(power, times_out, freqs_out; title="tf_superlet (5th order)")

# 5. Power spectrum
power_spectrum, freqs_out = eegfun.tf_spectrum(signal, sample_rate, foi; taper=:dpss, tapsmofrq=2.0)
fig4 = Figure(size=(600, 400))
ax = Axis(fig4[1, 1], xlabel="Frequency (Hz)", ylabel="Power", title="tf_spectrum (peaks at 10Hz, 25Hz)")
lines!(ax, freqs_out, power_spectrum)
display(fig4)
 

# power, t, f = eegfun.tf_analysis(:wavelet, signal, times, sample_rate, foi, toi; width=5)
# fig1, _ = plot_tf_result(power, times_out, freqs_out; title="tf_wavelet (7 cycles)")
# power, t, f = eegfun.tf_analysis(:multitaper, signal, times, sample_rate, foi, toi; taper=:dpss)
# fig1, _ = plot_tf_result(power, times_out, freqs_out; title="tf_wavelet (7 cycles)")
# power, f = eegfun.tf_analysis(:spectrum, signal, nothing, sample_rate, foi, nothing)
# fig1, _ = plot_tf_result(power, times_out, freqs_out; title="tf_wavelet (7 cycles)")a
# power, f = eegfun.tf_analysis(:spectrum, signal, nothing, sample_rate, foi, nothing)
# fig1, _ = plot_tf_result(power, times_out, freqs_out; title="tf_wavelet (7 cycles)")

 

    
# ============================================================================
# Test 2: Real EEG Data
# ============================================================================

# Try to load epoched data from output_data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
epoch_file = joinpath(data_dir, epoch_files[1])
data = eegfun.load_data(epoch_file)

   
# Handle Vector{EpochData} - use first condition
epochs = data isa Vector ? data[1] : data

# Get data from a central electrode
selected_channel = :Cz

# Extract signal data from epochs (samples × trials)
# epochs.data is Vector{DataFrame}, each DataFrame is one trial
n_trials = length(epochs.data)
n_samples = nrow(epochs.data[1])
signal = zeros(n_samples, n_trials)
for (i, trial_df) in enumerate(epochs.data)
    signal[:, i] = trial_df[:, selected_channel]
end

# Get time vector from first trial's time column, sample rate from epochs
times = collect(epochs.data[1].time)
sample_rate = epochs.sample_rate

# Time and frequency parameters
toi = -0.5:0.01:2
foi = 1:1:40

# 1. Wavelet
power, t, f = eegfun.tf_wavelet(signal, times, sample_rate, foi, collect(toi); width=7)
fig1, _ = plot_tf_result(power, t, f; title="Wavelet - $selected_channel") 

# 2. Wavelet with variable cycles
power, t, f = eegfun.tf_wavelet(signal, times, sample_rate, foi, collect(toi); width=[3, 10])
fig2, _ = plot_tf_result(power, t, f; title="Wavelet Variable - $selected_channel")

# 3. Superlet
power, t, f = eegfun.tf_superlet(signal, times, sample_rate, foi, collect(toi); order=5)
fig3, _ = plot_tf_result(power, t, f; title="Superlet - $selected_channel")

# 4. Multitaper Hanning
power, t, f = eegfun.tf_multitaper(signal, times, sample_rate, foi, collect(toi); taper=:hanning, t_ftimwin=0.3)
fig4, _ = plot_tf_result(power, t, f; title="MTM Hanning - $selected_channel")

# 5. Multitaper DPSS
power, t, f = eegfun.tf_multitaper(signal, times, sample_rate, foi, collect(toi); taper=:dpss, t_ftimwin=0.3, tapsmofrq=4.0)
fig5, _ = plot_tf_result(power, t, f; title="MTM DPSS - $selected_channel") 

# 6. Spectrum (non-time-resolved)
power_spec, f_spec = eegfun.tf_spectrum(signal, sample_rate, foi; taper=:dpss, tapsmofrq=2.0)
fig6 = Figure(size=(600, 400))
ax = Axis(fig6[1, 1], xlabel="Frequency (Hz)", ylabel="Power", title="Power Spectrum - $selected_channel")
lines!(ax, f_spec, power_spec)
display(fig6)
    
