###########################################################
# mapcols(x -> filtfilt(filter, x), dat.data[:, 3:end])
# mapcols!(x -> filtfilt(filter, x), dat.data[:, 3:end])
#
#
# TODO: Would it be useful to keep basefile name somewhere?
# Base.basename(splitext(dat.filename)[1]) 


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
