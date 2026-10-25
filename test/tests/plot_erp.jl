using Test
using Makie
using EegFun

@testset "plot_erp" begin
    dat = EegFun.create_test_erp_data(n_channels = 3)
    dat2 = EegFun.create_test_erp_data(n_channels = 3, condition = 2)

    @testset "single ErpData" begin
        result = EegFun.plot_erp(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "vector of ErpData" begin
        result = EegFun.plot_erp([dat, dat2]; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp(dat; channel_selection = EegFun.channels([1, 2]), display_plot = false)
        @test result.fig isa Figure
    end

    @testset "baseline interval" begin
        result = EegFun.plot_erp(dat; baseline_interval = (-0.1, 0.0), display_plot = false)
        @test result.fig isa Figure
    end

    @testset "legend kwargs and nbanks" begin
        leg_kw = Dict(:framevisible => false)
        result = EegFun.plot_erp([dat, dat2]; legend_kwargs = leg_kw, legend_nbanks = 2, layout = :grid, display_plot = false)
        @test result.fig isa Figure
        @test haskey(leg_kw, :framevisible)
    end

    @testset "empty data check" begin
        @test_throws EegFun.EegFunError EegFun.plot_erp(dat; condition_selection = EegFun.conditions("nonexistent"), display_plot = false)
    end

    @testset "channel_plot_order and AbstractVector views" begin
        # String channel_plot_order
        res = EegFun.plot_erp(dat; channel_plot_order = ["Ch2", "Ch1"], display_plot = false)
        @test res.fig isa Figure

        # AbstractVector view
        vec_dats = [dat, dat2]
        res_view = EegFun.plot_erp(@view(vec_dats[1:2]); channel_plot_order = [:Ch1, :Ch2], display_plot = false)
        @test res_view.fig isa Figure
    end
end
