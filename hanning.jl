# signal = DataFrame(CSV.File("data_hanning.csv")).data
# # reshap into 99 trials
# signal = reshape(signal, 640, 99)
# sample_rate = 256
# t = -1:1/sample_rate:(1.5-(1/sample_rate));
# dsignal = signal #detrend(t, signal)
# fig = Figure()
# ax = Axis(fig[1, 1])
# timewin = 100  # This is in milliseconds, but needs to be converted properly
# times2save = -0.5:0.01:1
# times2saveidx = zeros(length(times2save))
# frex = 2:2:30
# for i = 1:length(times2save)
#     times2saveidx[i] = find_idx_start_end(t, [times2save[i]])[1]
# end
# # Convert timewin from ms to samples and ensure it's odd for symmetric window
# timewinidx = round(Int, timewin / (1000 / sample_rate))
# timewinidx = timewinidx + (1 - timewinidx % 2)  # Make odd
# # Center the data window around each timepoint
# tf = zeros(length(frex), length(times2save))
# for timepointi = 1:length(times2save)
#     # Calculate window indices ensuring we don't go out of bounds
#     center_idx = Int(times2saveidx[timepointi])
#     half_win = div(timewinidx, 2)
#     start_idx = max(1, center_idx - half_win)
#     end_idx = min(size(signal, 1), center_idx + half_win)
#     # Extract and window the data
#     tmpdat = dsignal[start_idx:end_idx, :]
#     win_length = size(tmpdat, 1)
#     hann_win = hanning(win_length)
#     # Apply window and compute FFT
#     taperdat = tmpdat .* hann_win
#     fdat = fft(taperdat, 1)
#     # Average power over trials
#     tf[:, timepointi] = mean(abs2.(fdat[1:length(frex), :]), dims=2)
# end
# # Plot using contourf
# contourf!(ax, times2save, frex, transpose(log10.(tf)), 
#          #levels=-1.75:0.4:1.75,
#          colormap=:viridis)




signal = DataFrame(CSV.File("dataFIC.csv")).data
# Reshape into trials
signal = reshape(signal, 900, 76)
sample_rate = 300
t = -1:1/sample_rate:(2-(1/sample_rate))
dsignal = detrend(t, signal)
# Create figure
fig = Figure()
ax = Axis(fig[1, 1])
# Set up time-frequency analysis parameters
timewin = 1  # Window length in seconds
# Adjust time range to account for edge effects
half_win_time = timewin/2
times2save = (-0.5+half_win_time):0.01:(1.5-half_win_time)  # Smaller time steps
times2saveidx = zeros(length(times2save))
for i = 1:length(times2save)
    times2saveidx[i] = find_idx_start_end(t, [times2save[i]])[1]
end
# Convert timewin to samples and ensure it's odd
timewinidx = round(Int, timewin * sample_rate)
timewinidx = timewinidx + (1 - timewinidx % 2)  # Make odd
# Fixed padding - scaled with window size
pad_samples = max(500, nextpow(2, 4 * timewinidx))  # Ensure adequate padding for longer windows
# Create Hanning window
hann_win = hanning(timewinidx)
# Define frequencies of interest
frex = 2:1:30
freq_idx = round.(Int, (frex .* pad_samples) ./ sample_rate) .+ 1
# Initialize time-frequency matrix
tf = zeros(length(frex), length(times2save))
# Compute time-frequency decomposition
for timepointi = 1:length(times2save)
    center_idx = Int(times2saveidx[timepointi])
    half_win = div(timewinidx, 2)
    start_idx = max(1, center_idx - half_win)
    end_idx = min(size(dsignal, 1), center_idx + half_win)
    # Only process if we have a full window
    if (end_idx - start_idx + 1) == timewinidx
        tmpdat = dsignal[start_idx:end_idx, :]
        taperdat = tmpdat .* hann_win
        # Zero pad the data
        padded_data = vcat(taperdat, zeros(pad_samples - timewinidx, size(taperdat, 2)))
        fdat = fft(padded_data, 1)
        tf[:, timepointi] = mean(abs2.(fdat[freq_idx, :]), dims=2)
    end
end
# Plot using heatmap
heatmap!(ax, times2save, frex, transpose(log10.(tf)), 
         colormap=:viridis)
# Add title and labels
ax.title = "Time-Frequency Analysis"
ax.xlabel = "Time (s)"
ax.ylabel = "Frequency (Hz)"



signal = DataFrame(CSV.File("dataFIC.csv")).data
# Reshape into trials
signal = reshape(signal, 900, 76)
sample_rate = 300
t = -1:1/sample_rate:(2-(1/sample_rate))
dsignal = detrend(t, signal)

# Create figure
fig = Figure()
ax = Axis(fig[1, 1])

# Set up time-frequency analysis parameters
timewin = 0.5  # Window length in seconds
# Adjust time range to account for edge effects
half_win_time = timewin/2
times2save = (-0.5+half_win_time):0.01:(1.5-half_win_time)  # Smaller time steps
times2saveidx = Int.(find_idx_start_end(t, times2save))  # Vectorized

# Convert timewin to samples and ensure it's odd
timewinidx = round(Int, timewin * sample_rate)
timewinidx = timewinidx + (1 - timewinidx % 2)  # Make odd

# Fixed padding - scaled with window size
pad_samples = max(500, nextpow(2, 4 * timewinidx))  # Ensure adequate padding for longer windows

# Create Hanning window
hann_win = hanning(timewinidx)

# Define frequencies of interest
frex = 2:1:30
freq_idx = round.(Int, (frex .* pad_samples) ./ sample_rate) .+ 1

# Initialize time-frequency matrix and pre-allocate arrays
tf = zeros(length(frex), length(times2save))
padded_data = zeros(pad_samples, size(dsignal, 2))

# Compute time-frequency decomposition
@views for timepointi = 1:length(times2save)
    center_idx = times2saveidx[timepointi]
    half_win = div(timewinidx, 2)
    start_idx = max(1, center_idx - half_win)
    end_idx = min(size(dsignal, 1), center_idx + half_win)
    
    # Only process if we have a full window
    if (end_idx - start_idx + 1) == timewinidx
        # Use views to avoid allocations
        padded_data .= 0  # Reset padding
        padded_data[1:timewinidx, :] .= dsignal[start_idx:end_idx, :] .* hann_win
        fdat = fft(padded_data, 1)
        tf[:, timepointi] = mean(abs2.(fdat[freq_idx, :]), dims=2)
    end
end

# Plot using heatmap
heatmap!(ax, times2save, frex, transpose(log10.(tf)), 
         colormap=:viridis)

# Add title and labels
ax.title = "Time-Frequency Analysis"
ax.xlabel = "Time (s)"
ax.ylabel = "Frequency (Hz)"



# TODO: test DimensionalData
using CSV, DataFrames, FFTW, DSP, CairoMakie, DimensionalData

# Load data and reshape into trials with labeled dimensions
signal = DataFrame(CSV.File("dataFIC.csv")).data
n_trials = 76
signal_2d = reshape(signal, :, n_trials)  # Auto-compute the first dimension

# Define dimensions using DimensionalData
time_dim = Dim{:time}(range(-1, step=1/300, length=size(signal_2d, 1)))
trial_dim = Dim{:trial}(1:n_trials)
ds = DimArray(signal_2d, (time_dim, trial_dim); name="signal")

ds[time=Where(x -> x <= 1.0 && x >= 0.0), trial=Where(x -> x <= 10)]

# Detrend along the time dimension
ds_detrend = detrend(ds[Ti=1:end], dims=1)  # `Ti` is the time dimension selector

# Time-frequency parameters
timewin = 0.5  # Window length (seconds)
half_win = timewin / 2
times2save = range(-0.5 + half_win, 1.5 - half_win; step=0.01)
times2save_idx = findall(t -> t ∈ times2save, ds[Ti].val)  # Find matching time indices

# Convert timewin to samples (odd)
timewin_samp = round(Int, timewin * 300)
timewin_samp += iseven(timewin_samp) ? 1 : 0

# Precompute Hanning window and padding
hann_win = hanning(timewin_samp)
pad_samples = max(500, nextpow(2, 4 * timewin_samp))

# Frequencies of interest
frex = 2:30
freq_idx = @. round(Int, (frex * pad_samples) / 300) + 1

# Initialize time-frequency DimArray
tf_data = zeros(length(frex), length(times2save))
tf = DimArray(tf_data, (Dim(frex, :freq), Dim(times2save, :time)); name="power")

# Preallocate padded data buffer
padded_data = zeros(pad_samples, size(ds_detrend, 2))

# Compute time-frequency decomposition
@views for (i, t) in enumerate(times2save_idx)
    half_win = div(timewin_samp, 2)
    win_start = max(1, t - half_win)
    win_end = min(size(ds_detrend, 1), t + half_win)

    if (win_end - win_start + 1) == timewin_samp
        padded_data .= 0
        padded_data[1:timewin_samp, :] .= ds_detrend[Ti=win_start:win_end, :] .* hann_win
        fdat = fft(padded_data, 1)
        tf[Fi=1:end, Ti=i] .= mean(abs2.(fdat[freq_idx, :]), dims=2)
    end
end

# Plot using DimensionalData's automatic axis labeling
fig = Figure()
ax = Axis(fig[1, 1])
heatmap!(ax, tf[Ti=1:end, Fi=1:end], colormap=:viridis)
ax.xlabel = "Time (s)"
ax.ylabel = "Frequency (Hz)"
fig


