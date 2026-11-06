using Test
using Makie
using EegFun

@testset "plot_filter" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "plot_filter_response" begin
        # Creates a filter and plots response (IIR & FIR across filtfilt, filt_zerophase, filt)
        filter_info = EegFun.create_highpass_filter(0.1, 1000.0)
        fig = EegFun.plot_filter_response(filter_info; display_plot = false)
        @test fig.fig isa Figure

        fi_fir = EegFun.create_lowpass_filter(30.0, 256.0; filter_method = :fir)
        @test EegFun.plot_filter_response(fi_fir; filter_func = :filtfilt, display_plot = false).fig isa Figure
        @test EegFun.plot_filter_response(fi_fir; filter_func = :filt_zerophase, display_plot = false).fig isa Figure
        @test EegFun.plot_filter_response(fi_fir; filter_func = :filt, display_plot = false).fig isa Figure
    end

    @testset "plot_erp_filter_gui" begin
        # Zero-arg call (loads bundled example data)
        res0 = EegFun.plot_erp_filter_gui(display_plot = false)
        @test res0.fig isa Figure
        @test haskey(res0, :screen)

        # Single ErpData (including toggling LP and HP filters on across IIR/FIR and filter functions)
        erp_data = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
        res1 = EegFun.plot_erp_filter_gui(erp_data[1]; channel = :Cz, display_plot = false)
        @test res1.fig isa Figure
        toggles = filter(x -> x isa Makie.Toggle, res1.fig.content)
        menus = filter(x -> x isa Makie.Menu, res1.fig.content)
        @test length(toggles) >= 3
        @test length(menus) >= 3
        toggles[2].active[] = true
        toggles[3].active[] = true
        @test length(filter(p -> p isa Makie.Lines, res1.axes[1].scene.plots)) == 2

        # Exercise FIR mode + filt_zerophase / filt / filtfilt + mirror toggle
        menus[2].i_selected[] = 2 # FIR (Hamming)
        menus[3].i_selected[] = 2 # filt_zerophase (1-pass FIR)
        @test length(filter(p -> p isa Makie.Lines, res1.axes[1].scene.plots)) == 2
        menus[3].i_selected[] = 3 # filt (1-pass causal)
        @test length(filter(p -> p isa Makie.Lines, res1.axes[1].scene.plots)) == 2
        toggles[1].active[] = true # mirror enabled
        menus[2].i_selected[] = 1  # back to Butterworth (IIR) while filt/filt_zerophase tested
        menus[3].i_selected[] = 2
        @test length(filter(p -> p isa Makie.Lines, res1.axes[1].scene.plots)) == 2

        # Vector of ErpData
        res2 = EegFun.plot_erp_filter_gui(erp_data; channel = :Cz, display_plot = false)
        @test res2.fig isa Figure

        # Filepath string
        res3 = EegFun.plot_erp_filter_gui(EegFun.example_path("data/julia/erps/example1_erps.jld2"); channel = :Cz, display_plot = false)
        @test res3.fig isa Figure
    end
end
