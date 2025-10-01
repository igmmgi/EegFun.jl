using Test
using DataFrames
using eegfun
using Makie
using Statistics

@testset "filter" begin

    # Helper to create simple ContinuousData for filtering tests
    function create_test_data(; n::Int = 2000, fs::Int = 1000)
        t = collect(0:(n-1)) ./ fs
        # Signal: DC offset + 5 Hz + 100 Hz components
        x = 0.5 .+ sin.(2π .* 5 .* t) .+ 0.2 .* sin.(2π .* 100 .* t)
        y = 2 .* x .+ 0.1 .* randn(length(x))
        df = DataFrame(time = t, triggers = zeros(Int, n), A = x, B = y)
        layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        dat = eegfun.ContinuousData(copy(df, copycols=true), layout, fs, eegfun.AnalysisInfo())
        return dat
    end

    @testset "create_filter" begin
        fs = 1000.0
        # IIR low-pass
        fi_iir = eegfun.create_filter("lp", "iir", 40.0, fs; order = 4, transition_width = 0.1)
        @test fi_iir isa eegfun.FilterInfo
        @test fi_iir.filter_type == "lp"
        @test fi_iir.filter_method == "iir"
        @test fi_iir.cutoff_freq == 40.0
        @test fi_iir.sample_rate == fs
        @test fi_iir.order == 4
        @test fi_iir.n_taps === nothing

        # FIR high-pass
        fi_fir = eegfun.create_filter("hp", "fir", 1.0, fs; transition_width = 0.25)
        @test fi_fir isa eegfun.FilterInfo
        @test fi_fir.filter_type == "hp"
        @test fi_fir.filter_method == "fir"
        @test fi_fir.n_taps !== nothing
        @test fi_fir.n_taps % 2 == 1  # odd taps
        @test fi_fir.n_taps >= 101

        # FIR low-pass and tap sizing monotonicity
        fi_lp_wide = eegfun.create_filter("lp", "fir", 40.0, fs; transition_width = 0.2)
        fi_lp_narrow = eegfun.create_filter("lp", "fir", 40.0, fs; transition_width = 0.05)
        @test fi_lp_wide.n_taps !== nothing && fi_lp_narrow.n_taps !== nothing
        @test fi_lp_narrow.n_taps > fi_lp_wide.n_taps

        # Invalid parameters return nothing via @minimal_error
        @test eegfun.create_filter("bp", "iir", 10.0, fs) === nothing
        @test eegfun.create_filter("lp", "foo", 10.0, fs) === nothing
        @test eegfun.create_filter("lp", "iir", -1.0, fs) === nothing
        @test eegfun.create_filter("lp", "iir", fs, fs) === nothing
        @test eegfun.create_filter("lp", "iir", 10.0, 0.0) === nothing
        @test eegfun.create_filter("lp", "iir", 10.0, fs; order = 0) === nothing
    end

    @testset "filter_data! application and metadata" begin
        dat = create_test_data()
        dat_orig = copy(dat)

        # High-pass to remove DC; check mean is reduced towards ~0 for channel A
        eegfun.filter_data!(dat, "hp", 1.0; order = 1, filter_method = "iir", channel_selection = eegfun.channels([:A]))
        @test abs(mean(dat.data.A)) < abs(mean(dat_orig.data.A))
        # Only selected channel modified
        @test !all(dat.data.A .== dat_orig.data.A)
        @test all(dat.data.B .== dat_orig.data.B)
        # Analysis info updated
        @test dat.analysis_info.hp_filter == 1.0

        # Low-pass; update lp field and modify both channels when selecting both
        eegfun.filter_data!(dat, "lp", 30.0; order = 3, filter_method = "iir", channel_selection = eegfun.channels([:A, :B]))
        @test dat.analysis_info.lp_filter == 30.0
        @test !all(dat.data.B .== dat_orig.data.B)
    end

    @testset "non-mutating filter_data" begin
        dat = create_test_data()
        dat2 = eegfun.filter_data(dat, "hp", 1.0; order = 1)
        # Original unchanged
        @test all(dat.data.A .== create_test_data().data.A)
        # Copy modified
        @test !all(dat2.data.A .== dat.data.A)
    end

    @testset "no channels selected returns early" begin
        dat = create_test_data()
        dat_orig = copy(dat)
        # channel_selection picks none
        result = eegfun.filter_data!(dat, "hp", 1.0; channel_selection = eegfun.channels(Symbol[]))
        @test result === nothing
        # Data and analysis_info unchanged
        @test all(dat.data.A .== dat_orig.data.A)
        @test dat.analysis_info.hp_filter == dat_orig.analysis_info.hp_filter
        @test dat.analysis_info.lp_filter == dat_orig.analysis_info.lp_filter
    end

    @testset "single-pass vs zero-phase" begin
        dat1 = create_test_data()
        dat2 = copy(dat1)
        # Single-pass introduces phase; zero-phase differs from single-pass
        eegfun.filter_data!(dat1, "lp", 20.0; filter_func = "filt")
        eegfun.filter_data!(dat2, "lp", 20.0; filter_func = "filtfilt")
        @test !all(dat1.data.A .== dat2.data.A)
    end

    @testset "EpochData filtering" begin
        # Build two epochs from the same base
        base = create_test_data(; n = 1000, fs = 500)
        df1 = copy(base.data, copycols=true)
        df2 = copy(base.data, copycols=true)
        # Keep originals for comparison (distinct objects)
        df1o = copy(df1, copycols=true)
        df2o = copy(df2, copycols=true)
        # Add epoch identifiers to the versions that will be filtered
        df1.epoch = fill(1, nrow(df1))
        df2.epoch = fill(2, nrow(df2))
        ep = eegfun.EpochData([df1, df2], base.layout, base.sample_rate, eegfun.AnalysisInfo())
        eegfun.filter_data!(ep, "hp", 0.5)
        @test ep.analysis_info.hp_filter == 0.5
        @test !all(ep.data[1].A .== df1o.A)
        @test !all(ep.data[2].B .== df2o.B)
    end

    @testset "ErpData filtering" begin
        # Build ERP from base
        base = create_test_data(; n = 2000, fs = 1000)
        erp_df = select(base.data, [:time, :A, :B])
        erp = eegfun.ErpData(copy(erp_df, copycols=true), base.layout, base.sample_rate, eegfun.AnalysisInfo(), 25)
        erp_orig = copy(erp)
        eegfun.filter_data!(erp, "lp", 30.0; order = 3)
        @test erp.analysis_info.lp_filter == 30.0
        @test !all(erp.data.A .== erp_orig.data.A)
        # Non-mutating path
        erp2 = eegfun.filter_data(erp_orig, "hp", 0.5)
        @test erp_orig.analysis_info.hp_filter == 0.0  # unchanged
        @test erp2.analysis_info.hp_filter == 0.5
    end

    @testset "filter characteristics" begin
        fs = 1000.0
        fi = eegfun.create_filter("lp", "iir", 40.0, fs; order = 4, transition_width = 0.1)
        chars = eegfun.get_filter_characteristics(fi; npoints = 256)
        @test chars.filter_type == "lp"
        @test isapprox(chars.transition_width, 0.1; atol = 1e-6)
        @test any(abs.(chars.cutoff_freq_3db .- 40.0) .< 5.0)  # near cutoff
        @test chars.stopband_atten < -10  # should be attenuated
        # print helper should not error
        @test eegfun.print_filter_characteristics(fi; npoints = 128) === nothing
    end

    @testset "plot_filter_response (no display)" begin
        fs = 1000.0
        fi = eegfun.create_filter("hp", "iir", 1.0, fs)
        fig, axes = eegfun.plot_filter_response(fi; xscale = :linear, display_plot = false)
        @test fig isa Figure
        @test length(axes) == 3
        # log scale path
        fig2, axes2 = eegfun.plot_filter_response(fi; xscale = :log, display_plot = false)
        @test fig2 isa Figure
        @test length(axes2) == 3
    end

end


