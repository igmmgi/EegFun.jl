# Tests for IIR and FIR high-pass, low-pass, band-pass, and notch filtering.

using Test
using DataFrames
using CairoMakie
using Statistics
using JLD2
using EegFun

@testset "filter" begin

    @testset "create_filter" begin
        fs = 1000.0
        # IIR low-pass
        fi_iir = EegFun.create_lowpass_filter(40.0, fs; order = 4, transition_width = 0.1)
        @test fi_iir isa EegFun.FilterInfo
        @test fi_iir.filter_type == :lp
        @test fi_iir.filter_method == :iir
        @test fi_iir.cutoff_freq == 40.0
        @test fi_iir.sample_rate == fs
        @test fi_iir.order == 4
        @test isnothing(fi_iir.n_taps)

        # FIR high-pass
        fi_fir = EegFun.create_highpass_filter(1.0, fs; filter_method = :fir, transition_width = 0.25)
        @test fi_fir isa EegFun.FilterInfo
        @test fi_fir.filter_type == :hp
        @test fi_fir.filter_method == :fir
        @test !isnothing(fi_fir.n_taps)
        @test fi_fir.n_taps % 2 == 1  # odd taps
        @test fi_fir.n_taps >= 101

        # FIR low-pass and tap sizing monotonicity
        fi_lp_wide = EegFun.create_lowpass_filter(40.0, fs; filter_method = :fir, transition_width = 0.2)
        fi_lp_narrow = EegFun.create_lowpass_filter(40.0, fs; filter_method = :fir, transition_width = 0.05)
        @test !isnothing(fi_lp_wide.n_taps) && !isnothing(fi_lp_narrow.n_taps)
        @test fi_lp_narrow.n_taps > fi_lp_wide.n_taps

    end

    @testset "filter_data! application and metadata" begin
        dat = EegFun.create_test_continuous_data()
        dat_orig = copy(dat)

        # High-pass to remove DC; check mean is reduced towards ~0 for channel Ch1
        EegFun.highpass_filter!(dat, 1.0; order = 1, filter_method = :iir, channel_selection = EegFun.channels([:Ch1]))
        @test abs(mean(dat.data.Ch1)) < abs(mean(dat_orig.data.Ch1))
        # Only selected channel modified
        @test !all(dat.data.Ch1 .== dat_orig.data.Ch1)
        @test all(dat.data.Ch2 .== dat_orig.data.Ch2)
        # Analysis info updated
        @test dat.analysis_info.hp_filter == 1.0

        # Low-pass; update lp field and modify both channels when selecting both
        EegFun.lowpass_filter!(dat, 30.0; order = 3, filter_method = :iir, channel_selection = EegFun.channels([:Ch1, :Ch2]))
        @test dat.analysis_info.lp_filter == 30.0
        @test !all(dat.data.Ch2 .== dat_orig.data.Ch2)
    end

    @testset "non-mutating lowpass_filter" begin
        dat = EegFun.create_test_continuous_data()
        dat_orig = copy(dat)
        dat2 = EegFun.lowpass_filter(dat, 30.0; order = 3)
        # Original unchanged
        @test all(dat.data.Ch1 .== dat_orig.data.Ch1)
        # Copy modified
        @test !all(dat2.data.Ch1 .== dat.data.Ch1)
    end

    @testset "no channels selected returns early" begin
        dat = EegFun.create_test_continuous_data()
        dat_orig = copy(dat)
        # channel_selection picks none
        result = EegFun.highpass_filter!(dat, 1.0; channel_selection = EegFun.channels(Symbol[]))
        @test result === dat
        # Data and analysis_info unchanged
        @test all(dat.data.Ch1 .== dat_orig.data.Ch1)
        @test dat.analysis_info.hp_filter == dat_orig.analysis_info.hp_filter
        @test dat.analysis_info.lp_filter == dat_orig.analysis_info.lp_filter
    end

    @testset "single-pass vs zero-phase" begin
        dat1 = EegFun.create_test_continuous_data()
        dat2 = copy(dat1)
        # Single-pass introduces phase; zero-phase differs from single-pass
        EegFun.lowpass_filter!(dat1, 20.0; filter_func = :filt)
        EegFun.lowpass_filter!(dat2, 20.0; filter_func = :filtfilt)
        @test !all(dat1.data.Ch1 .== dat2.data.Ch1)
    end

    @testset "EpochData filtering" begin
        # Build two epochs from the same base
        base = EegFun.create_test_continuous_data(; n = 1000, fs = 500)
        df1 = copy(base.data, copycols = true)
        df2 = copy(base.data, copycols = true)
        # Keep originals for comparison (distinct objects)
        df1o = copy(df1, copycols = true)
        df2o = copy(df2, copycols = true)
        # Add epoch identifiers to the versions that will be filtered
        df1.epoch = fill(1, nrow(df1))
        df2.epoch = fill(2, nrow(df2))
        ep = EegFun.EpochData(base.file, 1, "condition_1", [df1, df2], base.layout, base.sample_rate, EegFun.AnalysisInfo())
        EegFun.highpass_filter!(ep, 0.5)
        @test ep.analysis_info.hp_filter == 0.5
        @test !all(ep.data[1].Ch1 .== df1o.Ch1)
        @test !all(ep.data[2].Ch1 .== df2o.Ch1)
    end

    @testset "ErpData filtering" begin
        # Build ERP from base
        base = EegFun.create_test_continuous_data(; n = 2000, fs = 1000)
        erp_df = select(base.data, [:time, :Ch1, :Ch2])
        erp = EegFun.ErpData(
            base.file,
            1,
            "condition_1",
            copy(erp_df, copycols = true),
            base.layout,
            base.sample_rate,
            EegFun.AnalysisInfo(),
            25,
        )
        erp_orig = copy(erp)
        EegFun.lowpass_filter!(erp, 30.0; order = 3)
        @test erp.analysis_info.lp_filter == 30.0
        @test !all(erp.data.Ch1 .== erp_orig.data.Ch1)
        # Non-mutating path
        erp2 = EegFun.highpass_filter(erp_orig, 0.5)
        @test erp_orig.analysis_info.hp_filter == 0.0  # unchanged
        @test erp2.analysis_info.hp_filter == 0.5
    end

    @testset "Vector of EegData filtering" begin
        base1 = EegFun.create_test_continuous_data(; n = 1000, fs = 500)
        base2 = EegFun.create_test_continuous_data(; n = 1000, fs = 500)
        vec_dats = [base1, base2]

        filtered_vec = EegFun.lowpass_filter(vec_dats, 30.0)
        @test length(filtered_vec) == 2
        @test filtered_vec[1].analysis_info.lp_filter == 30.0
        @test filtered_vec[2].analysis_info.lp_filter == 30.0
        @test vec_dats[1].analysis_info.lp_filter == 0.0

        EegFun.highpass_filter!(vec_dats, 1.0)
        @test vec_dats[1].analysis_info.hp_filter == 1.0
        @test vec_dats[2].analysis_info.hp_filter == 1.0
    end

    @testset "filter characteristics" begin
        fs = 1000.0
        fi = EegFun.create_lowpass_filter(40.0, fs; order = 4, transition_width = 0.1)
        chars = EegFun.get_filter_characteristics(fi; npoints = 256)
        @test chars.filter_type == :lp
        @test isapprox(chars.transition_band, 4.0; atol = 1e-6)  # 40 Hz * 0.1 = 4.0 Hz
        @test any(abs.(chars.cutoff_freq_3db .- 40.0) .< 5.0)  # near cutoff
        @test chars.stopband_atten < -10  # should be attenuated
        # print helper should not error
        @test isnothing(EegFun.print_filter_characteristics(fi; npoints = 128))
    end

    @testset "plot_filter_response (no display)" begin
        fs = 1000.0
        fi = EegFun.create_highpass_filter(1.0, fs)
        res1 = EegFun.plot_filter_response(fi; xscale = :linear, display_plot = false)
        fig = res1 isa NamedTuple ? res1.fig : res1[1]
        axes = res1 isa NamedTuple ? res1.axes : res1[2]
        @test fig isa Figure
        @test length(axes) == 2
        # log scale path
        res2 = EegFun.plot_filter_response(fi; xscale = :log, display_plot = false)
        fig2 = res2 isa NamedTuple ? res2.fig : res2[1]
        axes2 = res2 isa NamedTuple ? res2.axes : res2[2]
        @test fig2 isa Figure
        @test length(axes2) == 2
    end

    @testset "FIR transition width in Hz, odd tap count, and zero-phase filtfilt dB" begin
        # 1. FIR highpass at 0.5 Hz (fs = 250 Hz): default transition width is 0.5 Hz
        fs = 250
        fi_hp = EegFun.create_highpass_filter(0.5, fs; filter_method = :fir)
        @test fi_hp.transition_band == 0.5
        @test fi_hp.n_taps == 2001 # odd tap count ceil(4 * 250 / 0.5) + 1 = 2001

        # 6100 samples: (6100 - 1) ÷ 3 = 2033 >= 2001 fits without widening transition band
        dat_fir_fit = EegFun.create_test_continuous_data(; n = 6100, fs = fs)
        filtered_fir = EegFun.highpass_filter(dat_fir_fit, 0.5; filter_method = :fir)
        @test filtered_fir.analysis_info.hp_filter == 0.5

        # Short data (2200 samples -> max_odd_taps = 733) widens transition_band cleanly in Hz without crashing filtfilt
        dat_fir_short = EegFun.create_test_continuous_data(; n = 2200, fs = fs)
        filtered_short = @test_logs (:warn, r"Widening") match_mode = :any EegFun.highpass_filter(dat_fir_short, 0.5; filter_method = :fir)
        @test filtered_short.analysis_info.hp_filter == 0.5

        # 2. Zero-phase filtfilt doubles dB attenuation relative to single-pass filt in get_filter_characteristics
        fi_lp = EegFun.create_lowpass_filter(30.0, 1000.0; order = 2)
        ch_ff = EegFun.get_filter_characteristics(fi_lp; npoints = 512, filter_func = :filtfilt)
        ch_f = EegFun.get_filter_characteristics(fi_lp; npoints = 512, filter_func = :filt)
        # At the single-pass -3 dB frequency, zero-phase filtfilt has -6 dB (exact 2x dB scaling)
        @test isapprox(ch_ff.cutoff_freq_6db, ch_f.cutoff_freq_3db; rtol = 1e-6)
        # Stopband attenuation in dB is nearly doubled (much deeper than single-pass)
        @test ch_ff.stopband_atten < 1.8 * ch_f.stopband_atten

        # 3. Empty stopband_mask guard (cutoff very close to Nyquist so cutoff + transition >= Nyquist)
        fi_near_nyq = EegFun.create_lowpass_filter(480.0, 1000.0; order = 2, transition_width = 0.2)
        ch_near_nyq = EegFun.get_filter_characteristics(fi_near_nyq; npoints = 128)
        @test isnan(ch_near_nyq.stopband_atten)

        # 4. One-pass zero-phase FIR (:filt_zerophase, MNE/EEGLAB style) preserves exact peak latency & linear ramps
        fs_zp = 256
        n_zp = 385
        dat_zp = EegFun.create_test_continuous_data(; n = n_zp, fs = fs_zp)
        t_zp = (0:(n_zp-1)) ./ fs_zp
        # Synthetic Gaussian ERP peak centered at sample 180 + high-frequency 60 Hz ripple
        peak_idx = 180
        clean_peak = 10.0 .* exp.(-0.5 .* ((t_zp .- t_zp[peak_idx]) ./ 0.04) .^ 2)
        dat_zp.data.Ch1 .= clean_peak .+ 2.0 .* sin.(2π .* 60.0 .* t_zp)
        # Linear ramp on Ch2 to verify C1 odd-reflection boundary preservation
        dat_zp.data.Ch2 .= 3.0 .+ 5.0 .* t_zp

        dat_fir_zp = EegFun.lowpass_filter(dat_zp, 30.0; filter_method = :fir, filter_func = :filt_zerophase)
        dat_fir_causal = EegFun.lowpass_filter(dat_zp, 30.0; filter_method = :fir, filter_func = :filt)

        # Zero-phase 1-pass FIR keeps the peak at exact sample 180 (0 samples shift), whereas causal filt shifts right by ~68 samples
        @test argmax(dat_fir_zp.data.Ch1) == peak_idx
        @test argmax(dat_fir_causal.data.Ch1) > peak_idx + 50
        # Odd-reflection preserves exact linear ramp across the entire signal including endpoints
        @test isapprox(dat_fir_zp.data.Ch2, dat_zp.data.Ch2; atol = 1e-10)

        # Short-ERP highpass FIR tap clamping works for both :filt_zerophase and :filt (no 10,241-tap wipeout)
        dat_hp_zp = @test_logs (:warn, r"Widening") match_mode = :any EegFun.highpass_filter(
            dat_zp,
            0.1;
            filter_method = :fir,
            filter_func = :filt_zerophase,
        )
        dat_hp_causal =
            @test_logs (:warn, r"Widening") match_mode = :any EegFun.highpass_filter(dat_zp, 0.1; filter_method = :fir, filter_func = :filt)
        @test std(dat_hp_zp.data.Ch1) > 1.0
        @test std(dat_hp_causal.data.Ch1) > 0.5

        # IIR + :filt_zerophase logs warning and falls back to zero-phase :filtfilt
        dat_iir_zp = @test_logs (:warn, r"requires filter_method = :fir") match_mode = :any EegFun.lowpass_filter(
            dat_zp,
            30.0;
            filter_method = :iir,
            order = 2,
            filter_func = :filt_zerophase,
        )
        dat_iir_ff = EegFun.lowpass_filter(dat_zp, 30.0; filter_method = :iir, order = 2, filter_func = :filtfilt)
        @test dat_iir_zp.data.Ch1 == dat_iir_ff.data.Ch1

        # Edge cases in _filt_fir_zerophase! (1-tap trivial kernel & reflect_limited when d >= n) and invalid filter_func
        dst_1tap = zeros(5)
        EegFun._filt_fir_zerophase!(dst_1tap, [2.0], [1.0, 2.0, 3.0, 4.0, 5.0])
        @test dst_1tap == [2.0, 4.0, 6.0, 8.0, 10.0]
        dst_short = zeros(3)
        EegFun._filt_fir_zerophase!(dst_short, [0.25, 0.5, 0.25, 0.0, 0.0][1:5], [4.0, 4.0, 4.0])
        @test length(dst_short) == 3
        fi_fir_char = EegFun.create_lowpass_filter(30.0, 256.0; filter_method = :fir)
        EegFun.print_filter_characteristics(fi_fir_char; filter_func = :filt_zerophase)
        @test_throws Exception EegFun.lowpass_filter(dat_zp, 30.0; filter_func = :invalid_func)
    end

end

@testset "Batch Filter" begin
    isdefined(@__MODULE__, :run_batch_jld2_contract_tests!) ||
        Base.include(@__MODULE__, joinpath(pkgdir(EegFun), "test", "tests", "test_helpers.jl"))

    run_batch_jld2_contract_tests!(
        (pat; kwargs...) -> EegFun.lowpass_filter(pat, 30.0; kwargs...);
        make_data = _ -> EegFun.create_test_batch_erp_data(n_conditions = 2),
        pattern = "erps",
        log_file = "filter_lp.log",
        log_snippets = ["filtering started", "lowpass", "cutoff: 30.0 Hz"],
    )

    mktempdir() do test_dir
        orig_erps = EegFun.create_test_batch_erp_data(n_conditions = 2)
        jldsave(joinpath(test_dir, "1_erps.jld2"); data = orig_erps)

        # High-pass batch filter & invalid cutoff error
        hp_out = joinpath(test_dir, "filtered_hp")
        res_hp = EegFun.highpass_filter("erps", 1.0; input_dir = test_dir, output_dir = hp_out)
        @test res_hp.success == 1 && res_hp.errors == 0
        @test load(joinpath(hp_out, "1_erps.jld2"), "data")[1] isa EegFun.ErpData
        @test_throws Exception EegFun.lowpass_filter("erps", -10.0; input_dir = test_dir)

        # Frequency attenuation & time preservation
        lp_out = joinpath(test_dir, "filtered_lp")
        EegFun.lowpass_filter("erps", 30.0; input_dir = test_dir, output_dir = lp_out)
        filtered_erps = load(joinpath(lp_out, "1_erps.jld2"), "data")
        @test std(filtered_erps[1].data.Ch1) < std(orig_erps[1].data.Ch1)
        @test filtered_erps[1].data.time == orig_erps[1].data.time

        # EpochData batch support
        ep_dir = mkpath(joinpath(test_dir, "epochs_in"))
        jldsave(
            joinpath(ep_dir, "1_epochs.jld2");
            data = [EegFun.create_test_epoch_data(condition = 1), EegFun.create_test_epoch_data(condition = 2)],
        )
        ep_out = joinpath(test_dir, "epochs_out")
        res_ep = EegFun.lowpass_filter("epochs", 30.0; input_dir = ep_dir, output_dir = ep_out)
        @test res_ep.success == 1 && res_ep.errors == 0
        filtered_epochs = load(joinpath(ep_out, "1_epochs.jld2"), "data")
        @test length(filtered_epochs) == 2 && filtered_epochs[1] isa EegFun.EpochData
    end

    @testset "Unwrapped phase_delay in get_filter_characteristics" begin
        fi_fir = EegFun.create_lowpass_filter(40.0, 1000.0; filter_method = :fir, order = 201)
        chars = EegFun.get_filter_characteristics(fi_fir; filter_func = :filt)
        @test isapprox(chars.phase_delay, chars.group_delay; rtol = 0.05)
    end
end
