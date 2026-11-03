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

    # 1. Dead CleanLineParams removed; CleanLineConfig is the active config struct
    @test !isdefined(EegFun, :CleanLineParams)
    cfg = EegFun.CleanLineConfig(line_frequencies = [50.0], bandwidth = 2.0)
    @test cfg isa EegFun.CleanLineConfig
    @test cfg.line_frequencies == [50.0]
    @test cfg.bandwidth == 2.0

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

    # 5. Ergonomic channel selection
    ergo_dat = EegFun.cleanline(
        dat;
        line_frequencies = [50.0],
        sliding_win_length = 2.0,
        sliding_win_step = 1.0,
        bandwidth = 2.0,
        channel_selection = ch1,
    )
    @test isapprox(ergo_dat.data[!, ch1], cl_dat.data[!, ch1], atol = 1e-10)

    # 6. Vector dispatch
    vec_cl = EegFun.cleanline(
        [dat, dat];
        line_frequencies = [50.0],
        sliding_win_length = 2.0,
        sliding_win_step = 1.0,
        bandwidth = 2.0,
        channel_selection = [ch1],
    )
    @test length(vec_cl) == 2
    @test isapprox(vec_cl[1].data[!, ch1], cl_dat.data[!, ch1], atol = 1e-10)

    vec_mut = [copy(dat), copy(dat)]
    EegFun.cleanline!(
        vec_mut;
        line_frequencies = [50.0],
        sliding_win_length = 2.0,
        sliding_win_step = 1.0,
        bandwidth = 2.0,
        channel_selection = [ch1],
    )
    @test isapprox(vec_mut[1].data[!, ch1], cl_dat.data[!, ch1], atol = 1e-10)
    @test isapprox(vec_mut[2].data[!, ch1], cl_dat.data[!, ch1], atol = 1e-10)

    # 7. Phase 2 regressions: tail samples cleaned, Nyquist harmonics skipped, k_tapers <= 1 guard, Float32/view inputs
    @testset "Tail window cleaning, Nyquist skip, k_tapers guard, and Real views" begin
        # Non-multiple length: 5.5s with win_length=2.0, win_step=1.0 leaves 0.5s tail after last regular step
        n_tail = Int(round(5.5 * fs))
        t_tail = (0:(n_tail-1)) ./ fs
        tail_dat = EegFun.create_test_continuous_data(n = n_tail, fs = fs, n_channels = 1, seed = 7)
        tail_ch = EegFun.channel_labels(tail_dat)[1]
        tail_dat.data[!, tail_ch] .= 15.0 .* sin.(2π .* 50.0 .* t_tail)

        cleaned_tail = EegFun.cleanline(
            tail_dat;
            line_frequencies = [50.0, 150.0], # 150 Hz >= fs/2 (128 Hz) must be skipped cleanly!
            sliding_win_length = 2.0,
            sliding_win_step = 1.0,
            p_value = 0.2,
        )
        # Verify tail samples (last 100 samples) had 50 Hz substantially attenuated
        @test std(cleaned_tail.data[!, tail_ch][(end-100):end]) < 0.5 * std(tail_dat.data[!, tail_ch][(end-100):end])

        # k_tapers <= 1 guard throws EegFunError instead of ZeroDivisionError in F-test DOF (2K - 2)
        @test_throws EegFun.EegFunError EegFun.cleanline(tail_dat; k_tapers = 1)
        @test_throws EegFun.EegFunError EegFun.cleanline(tail_dat; k_tapers = 0)

        # Relaxed _cleanline_channel! accepts Float32 SubArray view, Int sample_rate, and Vector{Int} line_frequencies
        buf32 = Float32.(15.0 .* sin.(2π .* 50.0 .* t_tail) .+ 0.2 .* randn(n_tail))
        v32 = @view buf32[1:end]
        df32 = DataFrame(time = collect(t_tail), sample = collect(1:n_tail), copycols = false)
        df32[!, tail_ch] = v32
        lay32 = EegFun.Layout(DataFrame(label = [tail_ch], inc = [0.0], azi = [0.0]), nothing, nothing, nothing)
        dat32 = EegFun.ContinuousData("f32_view", df32, lay32, Int(fs), EegFun.AnalysisInfo())
        EegFun.cleanline!(dat32; line_frequencies = [50], sliding_win_length = 2.0, sliding_win_step = 1.0, p_value = 0.2)
        @test std(v32) < Float32(5.0)
    end
end
