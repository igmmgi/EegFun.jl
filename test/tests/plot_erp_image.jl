using Test
using Makie
using EegFun

@testset "plot_erp_image" begin
    dat = EegFun.create_test_epoch_data(n_epochs = 20, n_channels = 3)

    @testset "basic functionality" begin
        result = EegFun.plot_erp_image(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp_image(dat; channel_selection = EegFun.channels([1]), display_plot = false)
        @test result.fig isa Figure
    end

end
