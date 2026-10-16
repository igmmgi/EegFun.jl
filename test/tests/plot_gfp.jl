using Test
using Makie
using EegFun

@testset "plot_gfp" begin
    dat = EegFun.create_test_erp_data(n_channels = 3)
    dat2 = EegFun.create_test_erp_data(n_channels = 3, condition = 2)

    @testset "single ErpData" begin
        fig = EegFun.plot_gfp(dat; display_plot = false)
        @test fig isa Figure
    end

    @testset "vector of ErpData" begin
        fig = EegFun.plot_gfp([dat, dat2]; display_plot = false)
        @test fig isa Figure
    end

    @testset "channel selection" begin
        fig = EegFun.plot_gfp(dat; channel_selection = EegFun.channels([1, 2]), display_plot = false)
        @test fig isa Figure
    end

    @testset "show_erp_traces" begin
        fig = EegFun.plot_gfp(dat; show_erp_traces = true, display_plot = false)
        @test fig isa Figure
        @test length(fig.content) >= 2
    end

    @testset "show_dissimilarity" begin
        fig = EegFun.plot_gfp(dat; show_dissimilarity = true, display_plot = false)
        @test fig isa Figure
        @test length(fig.content) >= 2
    end

    @testset "all panels" begin
        fig = EegFun.plot_gfp(dat; show_erp_traces = true, show_dissimilarity = true, display_plot = false)
        @test fig isa Figure
        @test length(fig.content) >= 3
    end

end
