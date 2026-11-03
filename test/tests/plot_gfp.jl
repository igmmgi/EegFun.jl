using Test
using Makie
using EegFun

@testset "plot_gfp" begin
    dat = EegFun.create_test_erp_data(n_channels = 3)
    dat2 = EegFun.create_test_erp_data(n_channels = 3, condition = 2)

    @testset "single ErpData" begin
        res = EegFun.plot_gfp(dat; color = nothing, linestyle = nothing, display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
        @test length(res.axes) == 1
        @test res.axes[1].ylabel[] == "GFP (%)"

        res_raw = EegFun.plot_gfp(dat; normalize = false, display_plot = false)
        @test res_raw.axes[1].ylabel[] == "GFP (μV)"
    end

    @testset "vector of ErpData" begin
        res = EegFun.plot_gfp([dat, dat2]; display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
    end

    @testset "channel selection" begin
        res = EegFun.plot_gfp(dat; channel_selection = EegFun.channels([1, 2]), display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
    end

    @testset "show_erp_traces" begin
        res = EegFun.plot_gfp(dat; show_erp_traces = true, display_plot = false)
        @test res.fig isa Figure
        @test length(res.axes) == 2
    end

    @testset "show_dissimilarity" begin
        res = EegFun.plot_gfp(dat; show_dissimilarity = true, display_plot = false)
        @test res.fig isa Figure
        @test length(res.axes) == 2
    end

    @testset "all panels" begin
        res = EegFun.plot_gfp(dat; show_erp_traces = true, show_dissimilarity = true, display_plot = false)
        @test res.fig isa Figure
        @test length(res.axes) == 3
    end

    @testset "GFP (μV) vs GFP (%) label heuristic" begin
        gfp_df = EegFun.gfp(dat; normalize = false)
        # Ensure values lie in 0..100 μV — must STILL label "GFP (μV)" when normalize == false!
        gfp_df.gfp .= clamp.(abs.(gfp_df.gfp) .+ 2.5, 1.0, 20.0)
        res_uv = EegFun.plot_gfp(gfp_df; normalize = false, display_plot = false)
        @test res_uv.axes[1].ylabel[] == "GFP (μV)"

        res_pct = EegFun.plot_gfp(gfp_df; normalize = true, display_plot = false)
        @test res_pct.axes[1].ylabel[] == "GFP (%)"
    end

end
