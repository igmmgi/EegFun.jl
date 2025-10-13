###########################################################
# mapcols(x -> filtfilt(filter, x), dat.data[:, 3:end])
# mapcols!(x -> filtfilt(filter, x), dat.data[:, 3:end])
#
#
# TODO: Would it be useful to keep basefile name somewhere?
# Base.basename(splitext(dat.filename)[1]) 

# # save / load
# save_object("$(subject)_continuous.jld2", dat)
# dat1 = load_object("3_continuous.jld2")




# older implementation
# function tf_morlet(signal, times, sample_rate, freqs, cycles; tois = nothing)
# 
#     # data dimensions
#     n_samples, n_trials = size(signal)
#     n_freq = length(freqs)
# 
#     tois_idx = 1:n_samples
#     if !isnothing(tois)
#         tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in tois]
#     end
# 
#     # reshape signal
#     signal = reshape(signal, length(signal), 1)
#     freqs = exp.(range(log(freqs[1]), log(freqs[end]), length = n_freq))
# 
#     wavelet_time = -2:(1/sample_rate):2
#     if length(cycles) == 1
#         cycles = cycles[1] ./ (2 * pi .* freqs) # variable cycles/frequency
#     else
#         cycles = exp.(range(log(cycles[1]), log(cycles[end]), length = n_freq)) ./ (2 * pi .* freqs) # variable cycles/frequency
#     end
# 
#     n_convolution = length(wavelet_time) + length(signal) - 1
#     n_conv_pow2 = nextpow(2, n_convolution)
#     half_of_wavelet_size = floor(Int, (length(wavelet_time) - 1) / 2)
#     signal_fft = fft([signal; zeros(n_conv_pow2 - length(signal))])
# 
#     tf_trials = zeros(n_freq, length(tois_idx), n_trials)
#     @inbounds for idx_freq in eachindex(freqs)
#         wavelet =
#             sqrt(1 ./ (cycles[idx_freq] .* sqrt(pi))) .* exp.(2 * im * pi * freqs[idx_freq] .* wavelet_time) .*
#             exp.(-wavelet_time .^ 2 ./ (2 * (cycles[idx_freq] .^ 2)))
# 
#         wavelet_fft = fft([wavelet; zeros(n_conv_pow2 - length(wavelet))])
# 
#         # convolution
#         eegconv = ifft(wavelet_fft .* signal_fft)[1:n_convolution][half_of_wavelet_size+1:end-half_of_wavelet_size]
# 
#         # reshape to original signal size
#         eegconv = reshape(eegconv, n_samples, n_trials)
# 
#         @views tf_trials[idx_freq, :, :] = abs2.(eegconv)[tois_idx, :]
#     end
# 
#     return tf_trials, times[tois_idx], freqs
# 
# end

# older implementation
# function tf_hanning(signal, times, sample_rate, frequencies, time_steps; window_length = nothing, cycles = nothing)
#     if isnothing(window_length) && isnothing(cycles)
#         throw(ArgumentError("Must specify either window_length or cycles"))
#     end
# 
#     # Convert time_steps to indices
#     tois_idx = [findfirst(≈(t, atol = (1000 / sample_rate) / 1000), times) for t in time_steps]
# 
#     # Pre-allocate output and temporary arrays
#     n_frex = length(frequencies)
#     n_timepoints = length(time_steps)
#     n_trials = size(signal, 2)
#     tf_trials = zeros(n_frex, n_timepoints, n_trials)
# 
#     # Pre-compute maximum window length for pre-allocation
#     max_timewinidx = if !isnothing(cycles)
#         round(Int, cycles / minimum(frequencies) * sample_rate)
#     else
#         round(Int, window_length * sample_rate)
#     end
#     max_timewinidx += iseven(max_timewinidx)
# 
#     # Pre-allocate reusable arrays outside all loops
#     tmpdat = zeros(max_timewinidx, n_trials)
#     fdat = zeros(ComplexF64, max_timewinidx, n_trials)
# 
#     for (fi, freq) in enumerate(frequencies)
#         # Calculate frequency-dependent or fixed window length
#         timewinidx = if !isnothing(cycles)
#             round(Int, cycles / freq * sample_rate)
#         else
#             round(Int, window_length * sample_rate)
#         end
#         timewinidx += iseven(timewinidx)
# 
#         # Pre-compute window and frequency index for this frequency
#         hann_win = 0.5 * (1 .- cos.(2π * (0:timewinidx-1) / (timewinidx - 1)))
#         frex_idx = round(Int, freq * timewinidx / sample_rate) + 1
# 
#         for (timepointi, center_idx) in enumerate(tois_idx)
#             # Calculate window indices
#             half_win = fld(timewinidx, 2)
#             start_idx = center_idx - half_win
#             end_idx = center_idx + half_win
# 
#             # Handle edge cases and data extraction in one step
#             fill!(view(tmpdat, 1:timewinidx, :), 0)
#             if start_idx < 1 || end_idx > size(signal, 1)
#                 pad_left = max(1 - start_idx, 0)
#                 valid_start = max(start_idx, 1)
#                 valid_end = min(end_idx, size(signal, 1))
#                 valid_length = valid_end - valid_start + 1
#                 tmpdat[(pad_left+1):(pad_left+valid_length), :] = view(signal, valid_start:valid_end, :)
#             else
#                 tmpdat[1:timewinidx, :] = view(signal, start_idx:end_idx, :)
#             end
# 
#             # Apply Hanning taper and compute FFT in-place
#             tmpdat[1:timewinidx, :] .*= hann_win
#             fdat[1:timewinidx, :] = tmpdat[1:timewinidx, :]
#             fft!(view(fdat, 1:timewinidx, :), 1)  # FFT along time dimension
# 
#             # Store power for all trials at once
#             tf_trials[fi, timepointi, :] = abs2.(fdat[frex_idx, :]) ./ timewinidx^2
#         end
#     end
# 
#     return tf_trials, times[tois_idx], frequencies
# end

# set_aog_theme!()
# fig = Figure()
# all_data = eegfun.all_data(epochs[1])
# mydata = stack(all_data, [:Fp1, :Fp2], variable_name = :channel, value_name = :value)
# plt =
#     data(mydata) *
#     mapping(:time => "Time [ms]", :value => "Amplitude [μV]", color = :channel => nonnumeric) *
#     visual(Lines) *
#     mapping(layout = :epoch => nonnumeric) 
# # plt = paginate(plt, layout = 4)
# # draw(plt, 2)
# draw(plt)


# Sample data
x = 0:0.1:10
y = sin.(x)
z = cos.(x)
fig1 = Figure()
ax = Axis(fig1[1, 1])
lines!(ax, x, y)
fig2 = Figure()
ax = Axis(fig2[1, 1])
lines!(ax, x, z)
display(GLMakie.Screen(), fig1)
display(GLMakie.Screen(), fig2)
