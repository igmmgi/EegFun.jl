using Test
using Makie
using EegFun

@testset "plot_filter" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "plot_filter_response" begin
        # Creates a filter and plots response
        filter_info = EegFun.create_highpass_filter(0.1, 1000.0)

        fig = EegFun.plot_filter_response(filter_info; display_plot = false)
        @test fig.fig isa Figure
    end
end
