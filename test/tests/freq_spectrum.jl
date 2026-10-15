using Test
using EegFun
using DSP
using DataFrames

@testset "Frequency Spectrum (Welch PSD)" begin
    fs = 256
    duration = 4.0
    n_samples = Int(duration * fs)
    t = range(0, stop = duration, length = n_samples)

    # 1. ContinuousData with synthetic frequencies
    c_dat = EegFun.create_test_continuous_data(n = n_samples, fs = fs, n_channels = 3, seed = 1)
    ch_labels = EegFun.channel_labels(c_dat)
    ch1 = ch_labels[1]
    ch2 = ch_labels[2]

    # Inject 15 Hz on ch1 and 40 Hz on ch2
    c_dat.data[!, ch1] .= 10.0 .* sin.(2π .* 15.0 .* t)
    c_dat.data[!, ch2] .= 10.0 .* sin.(2π .* 40.0 .* t)

    # Compute spectrum up to 60 Hz
    spec = EegFun.freq_spectrum(c_dat; max_freq = 60.0)

    @test spec isa EegFun.SpectrumData
    @test spec.sample_rate == fs
    @test spec.method == :welch
    @test maximum(spec.data.freq) <= 60.0
    @test :freq in propertynames(spec.data)
    @test ch1 in propertynames(spec.data)
    @test ch2 in propertynames(spec.data)

    # Check peak frequency detection
    peak_idx1 = argmax(spec.data[!, ch1])
    peak_freq1 = spec.data.freq[peak_idx1]
    @test abs(peak_freq1 - 15.0) <= 1.0

    peak_idx2 = argmax(spec.data[!, ch2])
    peak_freq2 = spec.data.freq[peak_idx2]
    @test abs(peak_freq2 - 40.0) <= 1.0

    # 2. Channel selection
    spec_ch1 = EegFun.freq_spectrum(c_dat; channel_selection = EegFun.channels(ch1))
    @test ch1 in propertynames(spec_ch1.data)
    @test !(ch2 in propertynames(spec_ch1.data))

    # 3. EpochData
    ep_dat = EegFun.create_test_epoch_data(n = 256, fs = fs, n_epochs = 5, n_channels = 3, seed = 1)
    spec_ep = EegFun.freq_spectrum(ep_dat)
    @test spec_ep isa EegFun.SpectrumData
    @test :freq in propertynames(spec_ep.data)

    # 4. ErpData
    erp_dat = EegFun.create_test_erp_data(n_channels = 3, seed = 1)
    spec_erp = EegFun.freq_spectrum(erp_dat)
    @test spec_erp isa EegFun.SpectrumData
    @test :freq in propertynames(spec_erp.data)

    # 5. Parameter validation
    @test_throws EegFun.EegFunError EegFun.freq_spectrum(c_dat; overlap = -0.1)
    @test_throws EegFun.EegFunError EegFun.freq_spectrum(c_dat; overlap = 1.0)
    @test_throws EegFun.EegFunError EegFun.freq_spectrum(c_dat; overlap = 1.5)
    @test_throws EegFun.EegFunError EegFun.freq_spectrum(c_dat; channel_selection = EegFun.channels(Symbol[]))
end
