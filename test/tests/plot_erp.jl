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
end
