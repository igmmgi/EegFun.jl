# Analyzin Neural TimeSeries Data
# Chapter 10
using LinearAlgebra
using DSP
# using BenchmarkTools
using GLMakie
using FFTW
using CSV
using DataFrames
using StatsBase

include("types.jl")
include("utils.jl")

a = [1, 2, 3]
b = [4, 5, 6]

@btime sum(a .* b)
@btime dot(a, b) # faster

# Figure 10.2
# TODO: length of my_conv result vs. DSP conv function?
function my_conv(signal, kernel)
    # step 1: pad signal with zeros at start/end size of kernel - 1
    kernel_length = length(kernel)
    signal = [zeros(kernel_length - 1); signal; zeros(kernel_length - 1)]
    # step 2: flip kernel 
    kernel = reverse(kernel)
    # step 3: compute dot product for all points in signal
    out = zeros(length(signal))
    @inbounds for idx = kernel_length:length(signal)
        @views out[idx] = dot(signal[idx-(kernel_length-1):idx], kernel)
    end
    # step 4: remove padding
    return out[kernel_length:(end-(kernel_length-1))]
end

signal = zeros(100)
signal[45:55] .= 1
kernel = LinRange(1, 0.2, 5)
fig = Figure()
ax1 = GLMakie.Axis(fig[1, 1])
lines!(ax1, signal)
xlims!(ax1, 0, 100)
ax2 = GLMakie.Axis(fig[2, 1])
lines!(ax2, kernel)
xlims!(ax2, 0, 100)
result = conv(signal, kernel) # from DSP (convolution)
# result = conv(signal, kernel)[3:end-2] # from DSP
# result = xcorr(signal, kernel) # from DSP (cross-correlation kernel NOT flipped) 
my_result = my_conv(signal, kernel)
ax3 = GLMakie.Axis(fig[3, 1])
lines!(ax3, result)
lines!(ax3, my_result)
xlims!(ax3, 0, 100)

# Figure 11.1
times = LinRange(0, 2, 2000)
freq = 10
amp = 2
signal = amp * sin.(2 * pi .* freq * times)
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, times, signal)
ax.xlabel = "Time (S)"
ax.ylabel = "Amplitude"

# Figure 11.2/11.3
times = LinRange(0, 2, 2000)
freqs = [5, 10, 15, 30]
amps = [5, 1, 2, 8]
fig = Figure()
signal_sum = zeros(length(times))
for (idx, plt) in enumerate(freqs)
    signal = amps[idx] * sin.(2 * pi .* freqs[idx] * times)
    signal_sum .+= signal
    ax = Axis(fig[idx, 1])
    ylims!(ax, -maximum(amps) * 1.1, maximum(amps) * 1.1)
    lines!(ax, times, signal)
    ax.xlabel = "Time (S)"
    ax.ylabel = "Amplitude"
end
ax = Axis(fig[length(freqs)+1, 1])
lines!(ax, times, signal_sum)
lines!(ax, times, signal_sum .+ (rand(length(signal_sum)) .- 0.5) .* 20, color = :black) # add some noise
ax.xlabel = "Time (S)"
ax.ylabel = "Amplitude"



# Figure 11.4
# naive implementation
function myfft(signal)
    signal_length = length(signal)
    fourier = ComplexF64[0.0 + 0.0im for _ = 1:signal_length]
    for fi = 1:signal_length
        # Generate the complex exponential
        sine_wave = exp.(-im * 2 * pi * (fi - 1) .* (0:(signal_length-1)) / signal_length)
        # Compute the dot product
        fourier[fi] = dot(sine_wave, signal)
    end
    return fourier
end

sample_rate = 1000
times = 0:(1/sample_rate):2
amps = [1, 0.5]
freqs = [3, 8]
hz = LinRange(0, sample_rate, length(times))
signal1 = amps[1] * sin.(2 .* pi .* freqs[1] .* times)
signal2 = amps[2] * sin.(2 .* pi .* freqs[2] .* times)
signal3 = signal1 + signal2
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, times, signal1)
xlims!(ax, 0, times[end])
ylims!(ax, -1.5, 1.5)
ax.xlabel = "Time [s]"
ax.ylabel = "Amplitude"
ax = Axis(fig[1, 2])
lines!(ax, times, signal2)
xlims!(ax, 0, times[end])
ylims!(ax, -1.5, 1.5)
ax.xlabel = "Time [s]"
ax.ylabel = "Amplitude"
ax = Axis(fig[1, 3])
lines!(ax, times, signal3)
xlims!(ax, 0, times[end])
ylims!(ax, -1.5, 1.5)
ax.xlabel = "Time [s]"
ax.ylabel = "Amplitude"
ax = Axis(fig[2, 1])
f = myfft(signal1) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)
ax.xlabel = "Frequency (Hz)"
ax.ylabel = "Amplitude"
ax = Axis(fig[2, 2])
f = myfft(signal2) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)
ax.xlabel = "Frequency (Hz)"
ax.ylabel = "Amplitude"
ax = Axis(fig[2, 3])
f = myfft(signal3) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)
ax.xlabel = "Frequency (Hz)"
ax.ylabel = "Amplitude"


# Figure 11.4
sample_rate = 500
times = 0:(1/sample_rate):2
amps = [1, 0.5]
freqs = [5, 8]
hz = LinRange(0, sample_rate, length(times))
signal1 = amps[1] * sin.(2 .* pi .* freqs[1] .* times)
signal2 = amps[2] * sin.(2 .* pi .* freqs[2] .* times)
signal3 = signal1 + signal2
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, signal1)
xlims!(ax, 0, length(times))
ylims!(ax, -1.5, 1.5)
ax = Axis(fig[1, 2])
lines!(ax, signal2)
xlims!(ax, 0, length(times))
ylims!(ax, -1.5, 1.5)
ax = Axis(fig[1, 3])
lines!(ax, signal3)
xlims!(ax, 0, length(times))
ylims!(ax, -1.5, 1.5)
ax = Axis(fig[2, 1])
f = fft(signal1) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)
ax = Axis(fig[2, 2])
f = fft(signal2) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)
ax = Axis(fig[2, 3])
f = fft(signal3) ./ length(times)
barplot!(ax, hz, abs.(f[1:length(hz)] .* 2))
xlims!(ax, 0, 10)
ylims!(ax, 0, 2)



# Figure 11.6 Extended
sample_rate = 200
nyquist_freq = sample_rate / 2
times = 0:(1/sample_rate):(1-(1/sample_rate))
data = randn(sample_rate)
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, times, data)
xlims!(ax, 0, times[end])
ylims!(ax, -3, 3)
fourier = ComplexF64[0.0 + 0.0im for _ = 1:length(data)]
frequencies = LinRange(0, sample_rate, length(times));
ax = Axis(fig[2, 1])
xlims!(ax, 0, times[end])
for ii = 1:length(data)
    sine_wave = exp.(-im * 2 * pi * (ii - 1) .* times)
    lines!(ax, times, real.(sine_wave))
    # Compute the dot product
    fourier[ii] = dot(sine_wave, data)
end
xlims!(ax, 0, times[end])
fourier = fourier ./ length(data)
ax = Axis(fig[3, 1])
barplot!(ax, frequencies, abs.(fourier[1:length(frequencies)] .* 2))
xlims!(ax, 0 - 0.5, frequencies[end] + 0.5)
# reconstruct data
reconstructed_data = zeros(length(data))
ax = Axis(fig[4, 1])
xlims!(ax, 0, times[end])
for ii = 1:length(data)
    sine_wave = fourier[ii] * exp.(-im .* 2 .* pi * (ii - 1) .* times)
    plot(real.(sine_wave))
    lines!(ax, times, real.(sine_wave))
    reconstructed_data = reconstructed_data .+ real.(sine_wave)
end
reconstructed_data == data
ax = Axis(fig[5, 1])
lines!(ax, times, reconstructed_data)
xlims!(ax, 0, times[end])
ylims!(ax, -3, 3)



# Figure 11.6 
sample_rate = 200
nyquist_freq = sample_rate / 2
times = 0:(1/sample_rate):(1-(1/sample_rate))
data = randn(sample_rate)
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, times, data)
xlims!(ax, 0, times[end])
ylims!(ax, -3, 3)
fourier = fft(data)
frequencies = LinRange(0, sample_rate, length(times));
ax = Axis(fig[2, 1])
xlims!(ax, 0, times[end])
ax = Axis(fig[2, 1])
barplot!(ax, frequencies, abs.(fourier[1:length(frequencies)] .* 2))
xlims!(ax, 0 - 0.5, frequencies[end] + 0.5)
# reconstruct data
reconstructed_data = ifft(fourier)
ax = Axis(fig[3, 1])
lines!(ax, times, real.(reconstructed_data))
xlims!(ax, 0, times[end])
ylims!(ax, -3, 3)



# Figure 11.10
sample_rate = 256
times = -1:(1/sample_rate):1
signal = DataFrame(CSV.File("/home/ian/Desktop/EEGfun/data1.csv")).data[1:513]
s = 5 / (2 * pi * 30);
kernel = exp.((-times .^ 2) / (2 .* s^2)) ./ 30
# Determine the size for FFT (length of convolution result)
n = length(signal) + length(kernel) - 1
# Zero-pad signals to length n for FFT
padded_signal = [signal; zeros(n - length(signal))]
padded_kernel = [kernel; zeros(n - length(kernel))]
# Perform convolution in the time domain
time_domain_convolution = conv(signal, kernel)
# Compute the FFT of the padded signals
freq_domain_signal = fft(padded_signal)
freq_domain_kernel = fft(padded_kernel)
# Perform multiplication in the frequency domain
freq_domain_product = freq_domain_signal .* freq_domain_kernel
# Compute the inverse FFT to return to the time domain
freq_to_time_domain = real(ifft(freq_domain_product))
# Align FFT-based convolution result to match time-domain convolution
freq_to_time_domain = freq_to_time_domain[1:length(time_domain_convolution)]
lines(time_domain_convolution)
lines!(freq_to_time_domain .+ 0.1)  # added y offset for visual 
result1 = ifft(fft(signal) .* fft(kernel))
result2 = ifft(fft(padded_signal) .* fft(padded_kernel))
result3 = conv(signal, kernel)#[floor(Int, n/4):end-(floor(Int, n/4)-1)]
lines(real.(result2))
lines!(real.(result3) .+ 0.1)


# Chapter 12: Morlet Wavelets

# Figure 12.1
sample_rate = 1000
times = -1:(1/sample_rate):1
freq = 4
sine_wave = cos.(2 * pi .* freq .* times)
# sine_wave = @. cos(2*pi * freq * time)
# lines(sine_wave)
n = 4 # number of wavelets
s = n / (2 * pi * freq)
gaussian_win = exp.(-times .^ 2 ./ (2 * s^2))
# lines(gaussian_win)
lines(times, sine_wave .* gaussian_win)

# Figure 12.4
sample_rate = 250
times = -1:(1/sample_rate):1
fig = Figure()
ax = Axis(fig[1, 1])
num_wavelets = 80
lowest_freq = 2
highest_freq = 100
freqs = LinRange(lowest_freq, highest_freq, num_wavelets)
for freq in freqs
    # sine_wave = cos.(2 * pi .* freq .* time)
    sine_wave = exp.(2 * im * pi .* freq .* times)
    println("Current frequency: ", freq)
    n = 6 # number of wavelet
    s = n / (2 * pi * freq)
    gaussian_win = exp.(-times .^ 2 ./ (2 * s^2))
    # lines(gaussian_win)
    lines!(ax, times, real.(sine_wave) .* gaussian_win)
end


# Figure 12.5
signal = DataFrame(CSV.File("/home/ian/Desktop/EEGfun/data1.csv")).data[1:513]
fig = Figure()
ax = Axis(fig[1, 1])
# create wavelet
sample_rate = 256
times = -1:(1/sample_rate):1
freq = 6
lines!(ax, signal)
sine_wave = exp.(2 * im * pi .* freq .* times)
s = 4.5 / (2 * pi * freq)
gaussian_win = exp.(-times .^ 2 ./ (2 * s^2))
wavelet = sine_wave .* gaussian_win
halfwaveletsize = ceil(Int, length(wavelet) / 2);
n_conv = length(wavelet) + length(signal) - 1
wavelet_padded = [wavelet; zeros(floor(Int, (n_conv - length(wavelet))))]
signal_padded = [signal; zeros(floor(Int, (n_conv - length(signal))))]
# half of the wavelet size, useful for chopping off edges after convolution.
fft_wavelet = fft(wavelet_padded)
fft_signal = fft(signal_padded)
ift = ifft(fft_wavelet .* fft_signal) .* sqrt(s) / 10;
wavelet_conv_data = real(ift[halfwaveletsize:(end-halfwaveletsize)]);
lines!(ax, wavelet_conv_data)
filter_high = digitalfilter(Highpass(4), Butterworth(2); fs = 256)
filter_low = digitalfilter(Lowpass(8), Butterworth(6); fs = 256)
signal_filtered = filtfilt(filter_high, signal)
signal_filtered = filtfilt(filter_low, signal_filtered)
lines!(ax, signal_filtered)




find_idx_start_end(time, limits) = findmin(abs.(time .- limits[1]))[2], findmin(abs.(time .- limits[end]))[2]



# Figure 13.11
function tf_morlet(signal, times, sample_rate, freqs, cycles; tois = nothing)

    # data dimensions
    n_samples, n_trials = size(signal)
    n_freq = length(freqs)

    tois_idx = 1:n_samples
    if !isnothing(tois)
        tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in tois]
    end

    # reshape signal
    signal = reshape(signal, length(signal), 1)
    freqs = exp.(range(log(freqs[1]), log(freqs[end]), length = n_freq))

    wavelet_time = -2:(1/sample_rate):2
    if length(cycles) == 1
        cycles = cycles[1] ./ (2 * pi .* freqs) # variable cycles/frequency
    else
        cycles = exp.(range(log(cycles[1]), log(cycles[end]), length = n_freq)) ./ (2 * pi .* freqs) # variable cycles/frequency
    end

    n_convolution = length(wavelet_time) + length(signal) - 1
    n_conv_pow2 = nextpow(2, n_convolution)
    half_of_wavelet_size = floor(Int, (length(wavelet_time) - 1) / 2)
    signal_fft = fft([signal; zeros(n_conv_pow2 - length(signal))])

    tf_trials = zeros(n_freq, length(tois_idx), n_trials)
    @inbounds for idx_freq in eachindex(freqs)
        wavelet =
            sqrt(1 ./ (cycles[idx_freq] .* sqrt(pi))) .* exp.(2 * im * pi * freqs[idx_freq] .* wavelet_time) .*
            exp.(-wavelet_time .^ 2 ./ (2 * (cycles[idx_freq] .^ 2)))

        wavelet_fft = fft([wavelet; zeros(n_conv_pow2 - length(wavelet))])

        # convolution
        eegconv = ifft(wavelet_fft .* signal_fft)[1:n_convolution][half_of_wavelet_size+1:end-half_of_wavelet_size]

        # reshape to original signal size
        eegconv = reshape(eegconv, n_samples, n_trials)

        @views tf_trials[idx_freq, :, :] = abs2.(eegconv)[tois_idx, :]
    end

    return tf_trials, times[tois_idx], freqs

end



function apply_tf_baseline_db(tf_data, times, baseline_window)
    base_idx = findall(x -> baseline_window[1] ≤ x ≤ baseline_window[2], times)
    baseline_power = mean(tf_data[:, base_idx], dims = 2)
    tf_corrected = 10 .* log10.(tf_data ./ baseline_power)
    return tf_corrected
end

# Figure 13.11 Panel A
signal = DataFrame(CSV.File("/home/ian/Desktop/EEGfun/data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1000:(1000/sample_rate):1500-(1000/sample_rate)) ./ 1000
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80, [3, 10])
tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
fig = Figure()
ax = Axis(fig[1, 1], yscale = log10)
contourf!(ax, times_out, freqs_out, transpose(tf), levels = -4:0.2:4, colormap = :jet)
xlims!(ax, -0.2, 1.0)
ylims!(ax, 2, 80)
ax.yticks = round.(freqs_out)
display(fig)


# Figure 13.11 Panel B
signal = DataFrame(CSV.File("/home/ian/Desktop/EEGfun/data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1000:(1000/sample_rate):1500-(1000/sample_rate)) ./ 1000
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80, [3, 10])
tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
fig = Figure()
ax = Axis(fig[1, 1])
contourf!(ax, times_out, freqs_out, transpose(tf), levels = -4:0.2:4, colormap = :jet)
xlims!(ax, -0.2, 1.0)
ylims!(ax, 2, 80)
ax.yticks = round.(freqs_out)
display(fig)


function generate_signal(time, n_trials, freqs, amps, times, noise)
    # Parameters
    n_samples = length(time)  # 900 samples per trial
    # Initialize signal (900 samples × 76 trials)
    signal_trials = zeros(n_samples, n_trials) + (rand(n_samples, n_trials) * noise)
    # Add frequency components to all trials
    for trial = 1:n_trials
        for signal in zip(freqs, amps, times)
            signal_trials[(time.>=signal[3][1]).&(time.<signal[3][2]), trial] .+=
                signal[2] .* sin.(2 * pi * signal[1] * time[(time.>=signal[3][1]).&(time.<signal[3][2])])
        end
    end
    return signal_trials
end

# simulated signal
sampling_rate = 256
duration = [-2 3]
times = duration[1]:((1/sampling_rate)):(duration[2]-((1/sampling_rate)))
signal =
    generate_signal(times, 1, [10, 15, 25, 40], [1, 1, 1, 1], [[0.0, 0.25], [0.25, 0.5], [0.5, 1.0], [1.0, 1.5]], 0)
lines(times, vec(mean(signal, dims = 2)))
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80)
tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
# tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
fig = Figure()
ax = Axis(fig[1, 1])
# contourf!(ax, times_out, freqs_out, transpose(tf), levels = -4:0.2:4, colormap = :jet)
contourf!(ax, times_out, freqs_out, transpose(tf), colormap = :jet)
xlims!(ax, duration[1], duration[2])
ylims!(ax, 0, 80)
ax.yticks = round.(freqs_out)
display(fig)




# Figure 13.11 Panel B
signal = DataFrame(CSV.File("/home/ian/Desktop/EEGfun/data2.csv")).data
signal = reshape(signal, 640, 99)
sample_rate = 256
times = (-1000:(1000/sample_rate):1500-(1000/sample_rate)) ./ 1000
tf_trials, times_out, freqs_out = tf_morlet(signal, times, 256, 1:1:80, [3])
tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
tf = apply_tf_baseline_db(tf, times_out, (-0.5, -0.2))
fig = Figure()
ax = Axis(fig[1, 1], yscale = log10)
contourf!(ax, times_out, freqs_out, transpose(tf), levels = -4:0.2:4, colormap = :jet)
xlims!(ax, -0.2, 1.0)
ylims!(ax, 2, 80)
ax.yticks = round.(freqs_out)
display(fig)
















signal = DataFrame(CSV.File("dataFlank.csv")).data
# Figure 13.11
function dataFlank(signal)
    min_freq = 1
    max_freq = 40
    num_freq = 40
    sample_rate = 256
    time = -2:(1/sample_rate):2
    eeg_time = -1000:(1000/sample_rate):2000
    #freqs = exp.(range(log(min_freq), log(max_freq), length = num_freq))
    freqs = 1:40#exp.(range(log(min_freq), log(max_freq), length = num_freq))
    #s = exp.(range(log(3), log(10), length = num_freq)) ./ (2 * pi .* freqs) # variable cycles/frequency
    s = 7 ./ (2 * pi .* freqs)
    n_wavelet = length(time)
    n_data = length(signal)
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)
    half_of_wavelet_size = floor(Int, (n_wavelet - 1) / 2)
    signal_padded = [signal; zeros(n_conv_pow2 - length(signal))]
    eegfft = fft(signal_padded)
    eegpower = zeros(num_freq, div(n_data, 160))  # seems 99 trials in exported mat file
    base_idx = find_idx_start_end(eeg_time, [-500 -100])
    @inbounds for fi = 1:num_freq
        wavelet =
            sqrt(1 ./ (s[fi] .* sqrt(pi))) .* exp.(2 * im * pi * freqs[fi] .* time) .*
            exp.(-time .^ 2 ./ (2 * (s[fi] .^ 2)))
        #lines(real.(wavelet))
        wavelet_padded = [wavelet; zeros(n_conv_pow2 - length(wavelet))]
        wavelet_fft = fft(wavelet_padded)
        # convolution
        eegconv = ifft(wavelet_fft .* eegfft)[1:n_convolution][half_of_wavelet_size+1:end-half_of_wavelet_size]
        # eegconv = eegconv[1:n_convolution][half_of_wavelet_size:end-halfwaveletsize];
        eegconv = reshape(eegconv, 769, 160)
        temppower = mean((abs.(eegconv)) .^ 2, dims = 2)
        # eegpower[fi, :] = temppower
        eegpower[fi, :] = 10 .* log10.(temppower ./ mean(temppower[base_idx[1]:base_idx[2]])) # db baseline
    end
    fig = Figure()
    ax = Axis(fig[1, 1])
    contourf!(ax, eeg_time, freqs, transpose(eegpower), levels = -5:0.2:5, colormap = :jet)
    # contourf!(ax, eeg_time, freqs, transpose(eegpower), colormap = :jet)
    ax.yticks = round.(freqs)
    ylims!(ax, 0, 40)
    display(fig)
    return fig, ax
end
fig1, ax = dataFlank(signal)


# TODO: loop over trials vs. concat data
# TODO:: polyremoval
signal = DataFrame(CSV.File("dataFIC.csv")).data
# Figure 13.11
function dataFIC(signal)
    min_freq = 1
    max_freq = 60
    num_freq = 119
    sample_rate = 300
    time = -100:((1000/sample_rate)/1000):100
    eeg_time = -1000:(1000/sample_rate):(2000-(1000/sample_rate))
    freqs = 1:0.5:60
    s = 5 ./ (2 * pi .* freqs)
    n_wavelet = length(time)
    n_data = length(signal)
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)
    half_of_wavelet_size = floor(Int, (n_wavelet - 1) / 2)
    signal_padded = [signal; zeros(n_conv_pow2 - length(signal))]
    eegfft = fft(signal_padded)
    eegpower = zeros(num_freq, div(n_data, 76))
    base_idx = find_idx_start_end(eeg_time, [-500 -100])
    @inbounds for fi = 1:num_freq
        wavelet =
            sqrt(1 ./ (s[fi] .* sqrt(pi))) .* exp.(2 * im * pi * freqs[fi] .* time) .*
            exp.(-time .^ 2 ./ (2 * (s[fi] .^ 2)))
        wavelet_padded = [wavelet; zeros(n_conv_pow2 - length(wavelet))]
        wavelet_fft = fft(wavelet_padded)
        eegconv = ifft(wavelet_fft .* eegfft)[1:n_convolution][half_of_wavelet_size+1:end-half_of_wavelet_size]
        eegconv = reshape(eegconv, 900, 76)
        temppower = mean((abs.(eegconv)) .^ 2, dims = 2)
        eegpower[fi, :] = 10 .* log10.(temppower ./ mean(temppower[base_idx[1]:base_idx[2]]))
    end
    fig = Figure()
    ax = Axis(fig[1, 1]) #, yscale = log10)
    contourf!(ax, eeg_time, freqs, transpose(eegpower), levels = -3:0.1:3, colormap = :jet)
    xlims!(ax, -500, 1500)
    #ax.yticks = round.(freqs)
    #ylims!(ax, 2, 80)
    display(fig)
    return fig, ax
end
fig1, ax = dataFIC(signal)








# Figure 13.12
frequency_good = 6
frequency_bad = 2
sample_rate = 500;
cycles = 4
time = -0.5:(1/sample_rate):0.5
wavelet_good =
    exp.(2 * im * pi * frequency_good .* time) .* exp.(-time .^ 2 ./ (2 .* (cycles / (2 * pi * frequency_good)) .^ 2))
wavelet_bad =
    exp.(2 * im * pi * frequency_bad .* time) .* exp.(-time .^ 2 ./ (2 .* (cycles / (2 * pi * frequency_bad)) .^ 2))
lines(real.(wavelet_good))
lines!(real.(wavelet_bad))


# Figure 13.13
frequency = 10;
sample_rate = 500;
t = -0.5:1/sample_rate:0.5;
numcycles = [3 7];
fig = Figure()
# make wavelet
wavelet1 = exp.(2 * im * pi * frequency .* t) .* exp.(-t .^ 2 ./ (2 * (numcycles[1] / (2 * pi * frequency)) .^ 2));
ax = Axis(fig[1, 1])
lines!(ax, t, real(wavelet1))
ax.xlabel = "Time"
ax.ylabel = "Amplitude"
ax.title = "Wavelet $(frequency) Hz at $(numcycles[1]) cycles"
wavelet2 = exp.(2 * im * pi * frequency .* t) .* exp.(-t .^ 2 ./ (2 * (numcycles[2] / (2 * pi * frequency)) .^ 2));
ax = Axis(fig[1, 2])
lines!(ax, t, real(wavelet2))
ax.xlabel = "Time"
ax.ylabel = "Amplitude"
ax.title = "Wavelet $(frequency) Hz at $(numcycles[1]) cycles"
ax = Axis(fig[2, 1:2])
hz_wav = LinRange(0, sample_rate / 2, round(Int, length(wavelet1) / 2) + 1)
fft_wavelet1 = 2 * abs.(fft(wavelet1))
fft_wavelet2 = 2 * abs.(fft(wavelet2))
lines!(ax, hz_wav, fft_wavelet1[1:length(hz_wav)])
lines!(ax, hz_wav, fft_wavelet2[1:length(hz_wav)])
xlims!(ax, 0, 50)


# Figure 15.1
function detrend(x, y)
    X = hcat(ones(length(x)), x)  # Design matrix (with intercept)
    β = X \ y  # Solve for coefficients (m, b)
    return y - (X * β)
end

# signal = DataFrame(CSV.File("data3.csv")).data
# sample_rate = 256
# t = -1:1/sample_rate:(1.5-(1/sample_rate));
# dsignal = detrend(t, signal)
# fig = Figure()
# ax = Axis(fig[1, 1])
# lines!(ax, t, signal)
# lines!(ax, t, dsignal)
# timewin = 1000
# timewinidx = round(Int, timewin / (1000 / sample_rate))
# hann_win = 0.5 * (1 .- cos.(2 * pi * (0:timewinidx-1) / (timewinidx - 1)))
# # lines(hann_win)
# ax = Axis(fig[2, 1])
# stime = find_idx_start_end(t, [0.50])[1]
# lines!(ax, t[stime:stime+(timewinidx-1)], dsignal[stime:stime+(timewinidx-1)])
# lines!(ax, t[stime:stime+(timewinidx-1)], dsignal[stime:stime+(timewinidx-1)] .* hann_win)
# dfft = fft(dsignal[stime:stime+timewinidx-1] .* hann_win)
# f = LinRange(0, div(sample_rate, 2), floor(Int, length(hann_win) / 2 + 1))
# ax = Axis(fig[3, 1])
# lines!(ax, f[2:end], abs.(dfft[2:floor(Int, length(hann_win) / 2)+1]) .^ 2);
# # create TF matrix and input column of data at selected time point
# tf = zeros(floor(Int, length(hann_win) / 2), length(signal));
# tf[:, stime+div(timewinidx, 2):stime+div(timewinidx, 2)+20] =
#     repeat(abs.(dfft[2:floor(Int, length(hann_win) / 2)+1] .* 2) * 2, 1, 21)
# ax = Axis(fig[4, 1])
# heatmap!(ax, t, f, transpose(log10.(tf .+ 1)))
# display(fig)


function tf_hanning(signal, time, sample_rate, frequencies, time_window, time_steps)
    # Convert time_steps to indices
    times2saveidx = [findfirst(≈(t), time) for t in time_steps]
    # Calculate window length in samples
    timewinidx = round(Int, time_window * sample_rate)
    # Ensure window length is odd to center the taper
    if iseven(timewinidx)
        timewinidx += 1
    end
    # Create Hanning taper
    hann_win = 0.5 * (1 .- cos.(2 * pi * (0:timewinidx-1) / (timewinidx - 1)))
    # Preallocate output (n_frequencies × n_time_points × n_trials)
    n_frex = length(frequencies)
    n_timepoints = length(time_steps)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)
    # Precompute frequency bins
    frex_idx = round.(Int, frequencies .* timewinidx ./ sample_rate) .+ 1
    # Perform time-frequency analysis trial-by-trial
    for trial = 1:n_trials
        for (timepointi, center_idx) in enumerate(times2saveidx)
            # Calculate start and end indices for the window (centered)
            start_idx = center_idx - fld(timewinidx, 2)
            end_idx = center_idx + fld(timewinidx, 2)
            # Handle edge cases (window extends beyond signal boundaries)
            if start_idx < 1 || end_idx > size(signal, 1)
                # Zero-pad the data if the window extends beyond the signal
                pad_left = max(1 - start_idx, 0)
                pad_right = max(end_idx - size(signal, 1), 0)
                tmpdat = zeros(timewinidx)
                valid_start = max(start_idx, 1)
                valid_end = min(end_idx, size(signal, 1))
                tmpdat[(pad_left+1):(end_idx-start_idx+1-pad_right)] = signal[valid_start:valid_end, trial]
            else
                # Extract data segment for this trial
                tmpdat = signal[start_idx:end_idx, trial]
            end
            # Apply Hanning taper
            taperdat = tmpdat .* hann_win
            # Compute FFT
            fdat = fft(taperdat) ./ timewinidx  # Normalize by window length
            # Compute power for each frequency
            tf_trials[:, timepointi, trial] .= abs.(fdat[frex_idx]) .^ 2
        end
    end
    return tf_trials
end
signal = DataFrame(CSV.File("dataFIC.csv")).data
# # reshap into 99 trials
signal = reshape(signal, 900, 76)
function generate_signal(sampling_rate, duration, n_trials, freqs, amps, times, noise)
    # Parameters
    t = duration[1]:((1000/sampling_rate)/1000):(duration[2]-((1000/sampling_rate)/1000))
    n_samples = length(t)  # 900 samples per trial
    # Initialize signal (900 samples × 76 trials)
    signal_trials = zeros(n_samples, n_trials) + (rand(n_samples, n_trials) * noise)
    # Add frequency components to all trials
    for trial = 1:n_trials
        println(trial)
        for signal in zip(freqs, amps, times)
            signal_trials[(t.>=signal[3][1]).&(t.<signal[3][2]), trial] .+=
                signal[2] .* sin.(2 * pi * signal[1] * t[(t.>=signal[3][1]).&(t.<signal[3][2])])
        end
    end
    return signal_trials
end
sample_rate = 300
signal = generate_signal(
    sample_rate,
    [-1, 2.0],
    76,
    [10, 15, 25, 5],
    [1, 1, 1, 1],
    [[0.0, 0.25], [0.25, 0.5], [0.5, 1.0], [1.0, 1.5]],
    0,
);
# time = -1:1/sample_rate:(2-(1/sample_rate));
times = -1:1/sample_rate:(2-(1/sample_rate));
frequencies = 1:1:50
time_window = 0.25
time_steps = -0.5:0.05:1.5
tf_trials = tf_hanning(signal, times, sample_rate, frequencies, time_window, time_steps)
# Average power across trials
tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
# Debugging: Print TF matrix min/max
println("TF matrix - Min: ", minimum(tf), " Max: ", maximum(tf))
# Plot time-frequency results
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Frequency (Hz)")
heatmap!(ax, time_steps, frequencies, transpose(tf), colormap = :jet, interpolate = false)  # Disable interpolation
# heatmap!(ax, times2save, frex, transpose(log10.(tf)), colormap = :jet, interpolate = false)  # Disable interpolation
# Colorbar(fig[1, 2], limits = (minimum(log10.(tf)), maximum(log10.(tf))), label = "Power (dB)")
#xlims!(ax, -0.5, 1.5)  # Match Fieldtrip's time range
#ylims!(ax, 1, 31)      # Match Fieldtrip's frequency range
ax.yticks = frequencies       # Explicitly set frequency ticks to match data points
display(fig)



#






# Frequency-dependent window length (7 cycles per window)
function tf_hanning_cycles(signal, time, sample_rate, frequencies, time_steps, num_cycles)
    # Frequency-dependent window length (7 cycles per window)
    times2saveidx = [findfirst(≈(t), time) for t in time_steps]
    tf_trials = zeros(length(frequencies), length(time_steps), size(signal, 2))  # Power for each trial
    # Perform time-frequency analysis trial-by-trial
    for trial = 1:size(signal, 2)
        for fi = 1:length(frequencies)
            # Calculate frequency-dependent window length
            timewin = num_cycles / frequencies[fi]  # Window length in seconds
            timewinidx = round(Int, timewin * sample_rate)  # Window length in samples
            # Ensure window length is odd to center the taper
            if iseven(timewinidx)
                timewinidx += 1
            end
            hann_win = 0.5 * (1 .- cos.(2 * pi * (0:timewinidx-1) / (timewinidx - 1)))  # Hanning taper
            for timepointi = 1:length(time_steps)
                # Calculate start and end indices for the window
                start_idx = Int(times2saveidx[timepointi] - fld(timewinidx, 2))
                end_idx = Int(times2saveidx[timepointi] + fld(timewinidx, 2))
                #start_idx = Int(times2saveidx[timepointi])
                #end_idx = Int(times2saveidx[timepointi] + timewinidx -1)
                # Handle edge cases (window extends beyond signal boundaries)
                if start_idx < 1 || end_idx > size(signal_trials, 1)
                    # Zero-pad the data if the window extends beyond the signal
                    pad_left = max(1 - start_idx, 0)
                    pad_right = max(end_idx - size(signal_trials, 1), 0)
                    tmpdat = zeros(timewinidx)
                    valid_start = max(start_idx, 1)
                    valid_end = min(end_idx, size(signal_trials, 1))
                    tmpdat[(pad_left+1):(end_idx-start_idx+1-pad_right)] = signal_trials[valid_start:valid_end, trial]
                else
                    # Extract data segment for this trial
                    tmpdat = signal_trials[start_idx:end_idx, trial]
                end
                # Apply Hanning taper
                taperdat = tmpdat .* hann_win
                # Compute FFT
                fdat = fft(taperdat) ./ timewinidx  # Normalize by window length
                #fdat = fftshift(fdat, 2)
                # Map the desired frequency to the correct FFT bin
                bin = round(Int, frequencies[fi] * timewinidx / sample_rate) + 1
                # Compute power for this trial
                tf_trials[fi, timepointi, trial] = abs.(fdat[bin]) .^ 2  # Use th
            end
        end
    end
    # Average power across trials
    tf = mean(tf_trials, dims = 3)  # Average over trials
    tf = dropdims(tf, dims = 3)  # Remove singleton dimension
    # Debugging: Print TF matrix min/max
    println("TF matrix - Min: ", minimum(tf), " Max: ", maximum(tf))
    # Plot time-frequency results
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
        title = "Time-Frequency Analysis (7 Cycles per Window)",
    )
    # heatmap!(ax, times2save, frex, transpose(log10.(tf)), colormap = :jet, interpolate = false)  # Disable interpolation
    heatmap!(ax, time_steps, frequencies, transpose(tf), colormap = :jet, interpolate = false)  # Disable interpolation
    #Colorbar(fig[1, 2], limits = (minimum(log10.(tf)), maximum(log10.(tf))), label = "Power (dB)")
    #xlims!(ax, -0.5, 1.5)  # Match Fieldtrip's time range
    #ylims!(ax, 2, 30)      # Match Fieldtrip's frequency range
    ax.yticks = frequencies       # Explicitly set frequency ticks to match data points
    display(fig)
end
signal = DataFrame(CSV.File("dataFIC.csv")).data
# # reshap into 99 trials
signal = reshape(signal, 900, 76)
function generate_signal(sampling_rate, duration, n_trials, freqs, amps, times, noise)
    # Parameters
    t = duration[1]:((1000/sampling_rate)/1000):(duration[2]-((1000/sampling_rate)/1000))
    n_samples = length(t)  # 900 samples per trial
    # Initialize signal (900 samples × 76 trials)
    signal_trials = zeros(n_samples, n_trials) + (rand(n_samples, n_trials) * noise)
    # Add frequency components to all trials
    for trial = 1:n_trials
        println(trial)
        for signal in zip(freqs, amps, times)
            signal_trials[(t.>=signal[3][1]).&(t.<signal[3][2]), trial] .+=
                signal[2] .* sin.(2 * pi * signal[1] * t[(t.>=signal[3][1]).&(t.<signal[3][2])])
        end
    end
    return signal_trials
end
sample_rate = 300
# signal = generate_signal(sample_rate, [-1, 2.0], 76, [10, 15, 25, 5], [1, 1, 1, 1], [[0.0, 0.25], [0.25, 0.5], [0.5, 1.0], [1.0, 1.5]], 0);
# time = -1:1/sample_rate:(2-(1/sample_rate));
time = -1:1/sample_rate:(2-(1/sample_rate));
frequencies = 1:1:50
time_window = 0.25
time_steps = -0.5:0.05:1.5
num_cycles = 10
tf_hanning_cycles(signal, time, sample_rate, frequencies, time_steps, num_cycles)




function multitaper(
    signal,
    time,
    sample_rate,
    frequencies;
    foi = 1:2:30,
    toi = -0.5:0.05:1.5,
    tapsmofrq = 0.4,
    padding = 1,
    do_detrend = true,
)
    # Calculate frequency-specific parameters
    time_windows = 5 ./ foi  # cfg.t_ftimwin
    time_bandwidth = tapsmofrq * foi  # cfg.tapsmofrq
    # Calculate number of tapers based on time-bandwidth product
    num_tapers = round.(Int, 2 * time_bandwidth .- 1)
    # Convert time steps to indices
    times2saveidx = [findfirst(≈(t), time) for t in toi]
    # Preallocate output
    n_frex = length(foi)
    n_timepoints = length(toi)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)
    # Main analysis loop
    for (fi, freq) in enumerate(foi)
        # Frequency-specific parameters
        current_time_bandwidth = time_bandwidth[fi]
        current_num_tapers = num_tapers[fi]
        timewinidx = round(Int, time_windows[fi] * sample_rate)
        # Skip this frequency band if num_tapers < 1
        if current_num_tapers < 1
            @warn "Skipping frequency $freq Hz: num_tapers < 1 (NW = $(current_time_bandwidth / 2))"
            continue
        end
        # Ensure odd window size
        if iseven(timewinidx)
            timewinidx += 1
        end
        # Generate tapers for this frequency using DSP.jl
        NW = current_time_bandwidth / 2  # NW = (time_window * bandwidth) / 2
        tapers = DSP.dpss(timewinidx, NW, current_num_tapers)
        # Time-frequency analysis for this frequency
        for trial = 1:n_trials
            for (timepointi, center_idx) in enumerate(times2saveidx)
                # Extract data segment with bounds checking
                half_win = fld(timewinidx, 2)
                start_idx = max(1, center_idx - half_win)
                end_idx = min(size(signal, 1), center_idx + half_win)
                tmpdat = signal[start_idx:end_idx, trial]
                # Handle padding if needed
                if length(tmpdat) < timewinidx
                    pad_size = timewinidx - length(tmpdat)
                    if start_idx == 1
                        tmpdat = vcat(zeros(pad_size), tmpdat)
                    else
                        tmpdat = vcat(tmpdat, zeros(pad_size))
                    end
                elseif length(tmpdat) > timewinidx
                    tmpdat = tmpdat[1:timewinidx]
                end
                if do_detrend
                    tmpdat = detrend(1:length(tmpdat), tmpdat)
                end
                # Apply tapers and compute FFT
                taper_power = 0.0
                for taper = 1:current_num_tapers
                    tapered_data = tmpdat .* tapers[:, taper]
                    # Handle padding
                    if padding > 1
                        padded_size = round(Int, padding * length(tapered_data))
                        tapered_data = vcat(tapered_data, zeros(padded_size - length(tapered_data)))
                    end
                    fdat = fft(tapered_data)
                    frex = (0:(length(fdat)-1)) .* (sample_rate / length(fdat))
                    idx = argmin(abs.(frex .- freq))
                    println("Frequency: ", freq, " Hz → FFT bin: ", frex[idx], " Hz")
                    taper_power += abs2(fdat[idx])
                end
                tf_trials[fi, timepointi, trial] = taper_power / current_num_tapers
            end
        end
    end
    return tf_trials
end

sample_rate = 300
signal = generate_signal(sample_rate, [-1, 2.0], 76, [10, 10, 20], [1, 1, 1], [[0.0, 0.25], [0, 0.25], [0, 0.25]], 0);
# time = -1:1/sample_rate:(2-(1/sample_rate));
time = -1:1/sample_rate:(2-(1/sample_rate));
time_steps = -0.5:0.05:1.5
time = -1:1/sample_rate:(2-(1/sample_rate))
frequencies = 1:2:30
time_steps = -0.5:0.05:1.5
@btime tf_trials = multitaper(signal, time, sample_rate, frequencies, toi = time_steps, tapsmofrq = 0.5)
# Average power across trials
@btime tf = mean(tf_trials, dims = 3)  # Average over trials
tf = dropdims(tf, dims = 3)  # Remove singleton dimension
# Debugging: Print TF matrix min/max
println("TF matrix - Min: ", minimum(tf), " Max: ", maximum(tf))
# Plot time-frequency results
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Frequency (Hz)")
heatmap!(ax, time_steps, frequencies, transpose(tf), colormap = :jet, interpolate = true)  # Disable interpolation
# heatmap!(ax, times2save, frex, transpose(log10.(tf)), colormap = :jet, interpolate = false)  # Disable interpolation
# Colorbar(fig[1, 2], limits = (minimum(log10.(tf)), maximum(log10.(tf))), label = "Power (dB)")
#xlims!(ax, -0.5, 1.5)  # Match Fieldtrip's time range
ylims!(ax, 0, 30)      # Match Fieldtrip's frequency range
ax.yticks = frequencies       # Explicitly set frequency ticks to match data points
display(fig)







function multitaper(
    signal,
    time,
    sample_rate,
    frequencies;
    foi = 1:2:30,
    toi = -0.5:0.05:1.5,
    tapsmofrq = 0.4,
    padding = 2,  # Increase padding for better frequency resolution
    do_detrend = true,
)
    # Calculate frequency-specific parameters
    time_windows = 5 ./ foi  # cfg.t_ftimwin
    time_bandwidth = tapsmofrq * foi  # cfg.tapsmofrq

    # Calculate number of tapers based on time-bandwidth product
    num_tapers = round.(Int, 2 * time_bandwidth .- 1)

    # Convert time steps to indices
    times2saveidx = [findfirst(≈(t), time) for t in toi]

    # Preallocate output
    n_frex = length(foi)
    n_timepoints = length(toi)
    n_trials = size(signal, 2)
    tf_trials = zeros(n_frex, n_timepoints, n_trials)

    # Main analysis loop (parallelized over frequencies)
    for fi = 1:n_frex
        # Frequency-specific parameters
        freq = foi[fi]
        current_time_bandwidth = time_bandwidth[fi]
        current_num_tapers = num_tapers[fi]
        timewinidx = round(Int, time_windows[fi] * sample_rate)

        # Skip this frequency band if num_tapers < 1
        if current_num_tapers < 1
            @warn "Skipping frequency $freq Hz: num_tapers < 1 (NW = $(current_time_bandwidth / 2))"
            continue
        end

        # Ensure odd window size
        if iseven(timewinidx)
            timewinidx += 1
        end

        # Generate tapers for this frequency using DSP.jl
        NW = current_time_bandwidth / 2  # NW = (time_window * bandwidth) / 2
        tapers = dpss(timewinidx, NW, current_num_tapers)

        # Precompute FFT frequencies
        padded_size = padding > 1 ? round(Int, padding * timewinidx) : timewinidx
        frex = (0:(padded_size-1)) .* (sample_rate / padded_size)

        # Time-frequency analysis for this frequency
        for trial = 1:n_trials
            for (timepointi, center_idx) in enumerate(times2saveidx)
                # Extract data segment with bounds checking
                half_win = fld(timewinidx, 2)
                start_idx = max(1, center_idx - half_win)
                end_idx = min(size(signal, 1), center_idx + half_win)
                tmpdat = signal[start_idx:end_idx, trial]

                # Ensure tmpdat has the same length as tapers
                if length(tmpdat) < timewinidx
                    pad_size = timewinidx - length(tmpdat)
                    if start_idx == 1
                        tmpdat = vcat(zeros(pad_size), tmpdat)
                    else
                        tmpdat = vcat(tmpdat, zeros(pad_size))
                    end
                elseif length(tmpdat) > timewinidx
                    tmpdat = tmpdat[1:timewinidx]
                end

                # Detrend if requested
                if do_detrend
                    tmpdat = detrend(1:length(tmpdat), tmpdat)
                end

                # Apply all tapers at once
                tapered_data = tmpdat .* tapers

                # Apply padding
                if padding > 1
                    tapered_data = vcat(tapered_data, zeros(padded_size - timewinidx, current_num_tapers))
                end

                # Compute FFT for all tapers
                fdat = fft(tapered_data, 1)
                idx = argmin(abs.(frex .- freq))
                taper_power = mean(abs2.(fdat[idx, :]))  # Average power across tapers

                tf_trials[fi, timepointi, trial] = taper_power
            end
        end
    end

    return tf_trials
end























# Figure 15.2
timewin = 400
times2save = -300:50:1000


# Baseline normalisation
# Decibel conversion

signal = DataFrame(CSV.File("data_baseline.csv")).data
# wavelet parameters
min_freq = 2;
max_freq = 128;
num_frex = 30;
sample_rate = 256;
# other wavelet parameters
freqs = exp.(range(log(min_freq), log(max_freq), length = num_freq))
time = -1:(1/sample_rate):1;
eeg_time = -1:(1/sample_rate):(1.5-(1/sample_rate))
half_of_wavelet_size = floor(Int, (length(time) - 1) / 2);
# FFT parameters (use next-power-of-2)
n_wavelet = length(time);
n_data = length(signal)
n_convolution = n_wavelet + n_data - 1;
n_conv_pow2 = nextpow(2, n_convolution)
wavelet_cycles = 4;
# FFT of data (note: this doesn't change on frequency iteration)
signal_padded = [signal; zeros(n_conv_pow2 - length(signal))]
fft_data = fft(signal_padded);
# initialize output time-frequency data
tf_data = zeros(length(freqs), n_data);
for fi = 1:length(freqs)
    # create wavelet and get its FFT
    wavelet =
        (pi * freqs[fi] * sqrt(pi))^-0.5 * exp.(2 * im * pi * freqs[fi] .* time) .*
        exp.(-time .^ 2 ./ (2 * (wavelet_cycles / (2 * pi * freqs[fi]))^2)) / freqs[fi]
    wavelet_padded = [wavelet; zeros(n_conv_pow2 - length(wavelet))]
    fft_wavelet = fft(wavelet_padded)
    # run convolution
    convolution_result_fft = ifft(fft_wavelet .* fft_data)
    convolution_result_fft = convolution_result_fft[1:n_convolution] # note: here we remove the extra points from the power-of-2 FFT
    convolution_result_fft = convolution_result_fft[half_of_wavelet_size+1:end-half_of_wavelet_size]
    # put power data into time-frequency matrix
    tf_data[fi, :] = abs.(convolution_result_fft) .^ 2
end
# plot results
# dB-correct
base_idx = find_idx_start_end(eeg_time, [-0.5 -0.2])
baseline_power = mean(tf_data[:, base_idx[1]:base_idx[2]], dims = 2);
tf_data_db_base = 10 .* log10.(tf_data ./ baseline_power)
fig = Figure()
ax = Axis(fig[1, 1], yscale = log10)
contourf!(ax, eeg_time, freqs, transpose(tf_data_db_base), levels = -25:1:25, colormap = :jet)
ax.yticks = round.(freqs)
xlims!(ax, -0.5, 1.5)
# percentage change
pctchange = 100 * (tf_data - repeat(baseline_power, 1, n_data)) ./ repeat(baseline_power, 1, n_data)
ax = Axis(fig[1, 2], yscale = log10)
contourf!(ax, eeg_time, freqs, transpose(pctchange), levels = -1000:100:1000, colormap = :jet)
ax.yticks = round.(freqs)
xlims!(ax, -0.5, 1.5)
# relative
baselinediv = tf_data ./ repeat(baseline_power, 1, n_data);
ax = Axis(fig[2, 1], yscale = log10)
contourf!(ax, eeg_time, freqs, transpose(baselinediv), levels = -10:1:10, colormap = :jet)
ax.yticks = round.(freqs)
xlims!(ax, -0.5, 1.5)
# Z-transform
baseline_power = tf_data[:, base_idx[1]:base_idx[2]]
baselineZ =
    (tf_data - repeat(mean(baseline_power, dims = 2), 1, n_data)) ./ repeat(std(baseline_power, dims = 2), 1, n_data);
ax = Axis(fig[2, 2], yscale = log10)
contourf!(ax, eeg_time, freqs, transpose(baselineZ), levels = -10:1:10, colormap = :jet)
ax.yticks = round.(freqs)
xlims!(ax, -0.5, 1.5)















# TODO: brush up on filters + plot filter response
# filters
fs = 1000
o = 2
f = 0.1
filter = digitalfilter(Highpass(f), Butterworth(o); fs = sample_rate)
H, w = freqresp(filter::FilterCoefficients)
lines(abs.(H))
# filters
fs = 1000
o = 20
f = 30
filter = digitalfilter(Lowpass(f), Butterworth(o); fs = sample_rate)
H, w = freqresp(filter::FilterCoefficients)
lines!(abs.(H))



function ft_morlet(
    signal::AbstractMatrix{<:Real};  # Trials in columns (time × trials)
    sample_rate::Real = 256,
    foi::AbstractVector{<:Real} = 1:2:30,
    toi::AbstractVector{<:Real} = -0.5:0.01:1.5,
    n_cycles::Real = 7,
)
    # Create time vector for wavelet
    time = -2:(1/sample_rate):2

    # Convert toi to sample indices
    toi_samples = round.(Int, toi .* sample_rate) .+ floor(Int, size(signal, 1) / 2)

    # Calculate wavelet parameters
    s = n_cycles ./ (2 * pi .* foi)
    n_wavelet = length(time)
    n_data = size(signal, 1)
    n_trials = size(signal, 2)
    n_convolution = n_wavelet + n_data - 1
    n_conv_pow2 = nextpow(2, n_convolution)
    half_of_wavelet_size = floor(Int, (n_wavelet - 1) / 2)

    # Initialize time-frequency matrix (frequencies × time × trials)
    tf_data = zeros(length(foi), length(toi), n_trials)

    # Main time-frequency decomposition
    @inbounds for trial = 1:n_trials
        # Prepare signal for current trial
        signal_padded = [signal[:, trial]; zeros(n_conv_pow2 - n_data)]
        eegfft = fft(signal_padded)

        for fi = 1:length(foi)
            # Create wavelet
            wavelet =
                sqrt(1 ./ (s[fi] .* sqrt(pi))) .* exp.(2 * im * pi * foi[fi] .* time) .*
                exp.(-time .^ 2 ./ (2 * (s[fi] .^ 2)))

            # FFT of wavelet
            wavelet_padded = [wavelet; zeros(n_conv_pow2 - length(wavelet))]
            wavelet_fft = fft(wavelet_padded)

            # Convolution
            eegconv = ifft(wavelet_fft .* eegfft)
            eegconv = eegconv[1:n_convolution][half_of_wavelet_size+1:end-half_of_wavelet_size]

            # Store power for this frequency and trial
            tf_data[fi, :, trial] = abs2.(eegconv[toi_samples])
        end
    end

    return (;
        tf_data = tf_data,  # (frequencies × time × trials)
        frequencies = foi,
        times = toi,
    )
end


