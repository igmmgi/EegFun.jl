
using EegFun
using CUDA
using BenchmarkTools
using Test

#######################################################################
@info EegFun.section("TEST 1: Synthetic Signal with Known Frequencies")
#######################################################################

# Generate synthetic signal
sample_rate = 1000.0
times, signal = EegFun.generate_signal(
    400,                                    # n_trials
    [-1.0, 3.0],                            # time_interval
    sample_rate,                            # sample_rate
    [5.0, 25, 35.0],                        # frequencies
    [5.0, 5.0, 5.0],                        # amplitudes
    [[0.1, 0.5], [0.6, 1.0], [1.1, 1.5]],   # time intervals for each freq 
    0.0,                                    # noise amplitude
);
epochs_synthetic = EegFun.signal_to_data(times, signal, :Channel1, sample_rate)

println("\n--- Correctness Tests (CPU vs GPU) ---")

function test_equiv(a, b)
    if size(a) != size(b)
        return false
    end
    max_diff = 0.0
    for i in eachindex(a)
        if isnan(a[i]) && isnan(b[i])
            continue
        elseif isnan(a[i]) || isnan(b[i])
            return false
        else
            diff = abs(a[i] - b[i]) / max(abs(a[i]), abs(b[i]), 1e-10)
            max_diff = max(max_diff, diff)
            if diff > 1e-4
                return false
            end
        end
    end
    return true
end

@testset "CPU vs GPU Equivalence" begin
    # tf_morlet
    tf_cpu = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, use_gpu = false)
    tf_gpu = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, use_gpu = true)
    @test test_equiv(tf_cpu.data_power[!, :Channel1], tf_gpu.data_power[!, :Channel1])
    
    tf_cpu_pad = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, pad = :both, use_gpu = false)
    tf_gpu_pad = EegFun.tf_morlet(epochs_synthetic, frequencies = 1:1:40, cycles = 3, pad = :both, use_gpu = true)
    @test test_equiv(tf_cpu_pad.data_power[!, :Channel1], tf_gpu_pad.data_power[!, :Channel1])

    # tf_multitaper
    tf_mt_cpu = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5, use_gpu = false)
    tf_mt_gpu = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5, use_gpu = true)
    @test test_equiv(tf_mt_cpu.data_power[!, :Channel1], tf_mt_gpu.data_power[!, :Channel1])

    tf_mt_cpu_pad = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5, pad = :both, use_gpu = false)
    tf_mt_gpu_pad = EegFun.tf_multitaper(epochs_synthetic, frequencies = 1:1:40, cycles = 5, pad = :both, use_gpu = true)
    @test test_equiv(tf_mt_cpu_pad.data_power[!, :Channel1], tf_mt_gpu_pad.data_power[!, :Channel1])

    # tf_stft
    tf_stft_cpu = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, use_gpu = false)
    tf_stft_gpu = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, use_gpu = true)
    @test test_equiv(tf_stft_cpu.data_power[!, :Channel1], tf_stft_gpu.data_power[!, :Channel1])

    tf_stft_cpu_pad = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, pad = :both, use_gpu = false)
    tf_stft_gpu_pad = EegFun.tf_stft(epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, pad = :both, use_gpu = true)
    @test test_equiv(tf_stft_cpu_pad.data_power[!, :Channel1], tf_stft_gpu_pad.data_power[!, :Channel1])
end

println("\n--- Benchmarks ---")

println("tf_morlet (CPU):")
@btime EegFun.tf_morlet($epochs_synthetic, frequencies = 1:1:40, cycles = 3, use_gpu = false)
println("tf_morlet (GPU):")
@btime EegFun.tf_morlet($epochs_synthetic, frequencies = 1:1:40, cycles = 3, use_gpu = true)

println("tf_multitaper (CPU):")
@btime EegFun.tf_multitaper($epochs_synthetic, frequencies = 1:1:40, cycles = 5, use_gpu = false)
println("tf_multitaper (GPU):")
@btime EegFun.tf_multitaper($epochs_synthetic, frequencies = 1:1:40, cycles = 5, use_gpu = true)

println("tf_stft (CPU):")
@btime EegFun.tf_stft($epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, use_gpu = false)
println("tf_stft (GPU):")
@btime EegFun.tf_stft($epochs_synthetic, frequencies = 1:1:40, window_length = 0.5, use_gpu = true)


