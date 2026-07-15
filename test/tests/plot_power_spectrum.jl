using Test
using Makie
using EegFun

@testset "plot_power_spectrum" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "plot_channel_spectrum" begin
        result = EegFun.plot_channel_spectrum(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "plot_channel_spectrum with kwargs" begin
        result = EegFun.plot_channel_spectrum(dat; max_freq = 50.0, display_plot = false)
        @test result.fig isa Figure
    end
end
