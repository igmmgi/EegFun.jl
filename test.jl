# Analyzin Neural TimeSeries Data
# Chapter 10
using LinearAlgebra
using DSP 
using BenchmarkTools
using GLMakie

a = [1, 2, 3]
b = [4, 5, 6]

@btime sum(a .* b)
@btime dot(a, b) # faster

# Figure 10.2
# TODO: length of my_conv result vs. DSP conv function?
function my_conv(signal, kernel)
  # step 1: pad signal with zeros at start/end size of kernel - 1
  kernel_length = length(kernel)
  signal = [zeros(kernel_length-1); signal; zeros(kernel_length-1)]
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
time = LinRange(0, 2, 2000)
freq = 20
amp =2
signal = amp*sin.(2*pi.*freq*time)
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, time, signal)
ax.xlabel = "Time (S)"
ax.ylabel = "Amplitude"

# Figure 11.2/11.3
time = LinRange(0, 2, 2000)
freqs = [5, 10, 15, 30]
amps = [5, 1, 2, 8]
fig = Figure()
signal_sum = zeros(length(time))
for (idx, plt) in enumerate(freqs)
  signal = amps[idx]*sin.(2*pi.*freqs[idx]*time)
  signal_sum .+= signal
  ax = Axis(fig[idx,1])
  ylims!(ax,-maximum(amps)*1.1, maximum(amps)*1.1 )
  lines!(ax, time, signal)
end
ax = Axis(fig[length(freqs)+1, 1])
lines!(ax, time, signal_sum)
lines!(ax, time, signal_sum .+ (rand(length(signal_sum)) .- 0.5 ).* 50) # add some noise


# Figure 11.4
N = 10 
data = rand(N)
fourier = zeros(N)
time = collect((1:N-1)/N)

fi = 5 

lines(real.(exp.(-1*im .*2 .* pi*(fi.-1) .* time)))


for fi = 1:N
  sine_wave = exp(-im.*2*pi*(fi.-1).*time)



