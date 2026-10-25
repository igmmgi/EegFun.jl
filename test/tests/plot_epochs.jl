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

    @testset "baseline_interval and legend_kwargs" begin
        res_bl = EegFun.plot_epochs(dat; baseline_interval = (0.0, 0.2), legend_kwargs = (; framevisible = false), layout = :grid, display_plot = false)
        @test res_bl.fig isa Figure
    end

    @testset "empty data check" begin
        @test_throws EegFun.EegFunError EegFun.plot_epochs(dat; condition_selection = EegFun.conditions("nonexistent"), display_plot = false)
    end

    @testset "Float32 data support" begin
        dat32 = EegFun.create_test_epoch_data(n_epochs = 4, n_channels = 2)
        for df in dat32.data
            for col in [:Ch1, :Ch2]
                df[!, col] = Float32.(df[!, col])
            end
        end
        res32 = EegFun.plot_epochs(dat32; display_plot = false)
        @test res32.fig isa Figure
    end
end
