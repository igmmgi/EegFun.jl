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

    @testset "grid layout" begin
        result = EegFun.plot_epochs(dat; layout = :grid, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) == 3
    end

    @testset "topo layout" begin
        result = EegFun.plot_epochs(dat; layout = :topo, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) >= 3
    end

    @testset "channel average" begin
        result = EegFun.plot_epochs(dat; average_channels = true, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) == 1
    end

    @testset "multi-condition vector" begin
        dat2 = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)
        result = EegFun.plot_epochs([dat, dat2]; layout = :grid, display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) == 3
    end

end
