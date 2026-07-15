using Test
using Makie
using EegFun

@testset "plot_epochs" begin
    dat = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)

    @testset "basic functionality" begin
        result = EegFun.plot_epochs(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel and epoch selection" begin
        result =
            EegFun.plot_epochs(dat; channel_selection = EegFun.channels([1]), epoch_selection = EegFun.epochs([1, 2]), display_plot = false)
        @test result.fig isa Figure
    end

end
