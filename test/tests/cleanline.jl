using Test
using EegFun
using Statistics

@testset "CleanLine Line Noise Removal" begin
    # Create 6 seconds of continuous data at 256 Hz
    fs = 256
    n_samples = 6 * fs
    t = range(0, stop = 6, length = n_samples)
    dat = EegFun.create_test_continuous_data(n = n_samples, fs = fs, n_channels = 3, seed = 42)

    ch_labels = EegFun.channel_labels(dat)
    ch1 = ch_labels[1]
    ch2 = ch_labels[2]
    ch3 = ch_labels[3]

    # Inject strong 50 Hz line noise into channel 1 only
    noise_50hz = 10.0 .* sin.(2π .* 50.0 .* t)
    dat.data[!, ch1] .+= noise_50hz

    # 1. CleanLineParams struct
    params = EegFun.CleanLineParams(line_frequencies = [50.0], bandwidth = 2.0)
    @test params isa EegFun.CleanLineParams
    @test params.line_frequencies == [50.0]
    @test params.bandwidth == 2.0

    # 2. Out-of-place cleanline
    orig_dat = copy(dat)
    cl_dat = EegFun.cleanline(
        dat;
        line_frequencies = [50.0],
        sliding_win_length = 2.0,
        sliding_win_step = 1.0,
        bandwidth = 2.0,
        channel_selection = EegFun.channels(ch1),
    )

    # Check input was not mutated
    @test dat.data[!, ch1] == orig_dat.data[!, ch1]
    # Check unselected channels are completely identical
    @test cl_dat.data[!, ch2] == dat.data[!, ch2]
    @test cl_dat.data[!, ch3] == dat.data[!, ch3]

    # Check 50Hz noise was reconstructed and removed
    diff_signal = orig_dat.data[!, ch1] .- cl_dat.data[!, ch1]
    @test cor(diff_signal, noise_50hz) > 0.95
    @test std(cl_dat.data[!, ch1]) < std(orig_dat.data[!, ch1])

    # 3. In-place cleanline!
    mut_dat = copy(dat)
    res = EegFun.cleanline!(
        mut_dat;
        line_frequencies = [50.0],
        sliding_win_length = 2.0,
        sliding_win_step = 1.0,
        bandwidth = 2.0,
        channel_selection = EegFun.channels(ch1),
    )
    @test res === mut_dat
    @test isapprox(mut_dat.data[!, ch1], cl_dat.data[!, ch1], atol = 1e-10)

    # 4. Empty channel selection returns dat and does nothing
    no_ch_dat = copy(dat)
    @test EegFun.cleanline!(no_ch_dat, channel_selection = EegFun.channels(Symbol[])) === no_ch_dat
    @test no_ch_dat.data == dat.data
end
