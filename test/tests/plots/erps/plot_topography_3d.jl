# Tests for 3D head-mesh scalp topography plotting.

using Test
using DataFrames
using CairoMakie
using EegFun

@testset "Plot Topography 3D" begin
    # Create test data with a layout that includes spherical coordinates
    dat = EegFun.create_test_erp_data(n_channels = 6)
    # Give them varying coordinates so they don't snap to the same point
    dat.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    dat.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]

    @testset "plot_topography_3d basic functionality" begin
        # Test basic rendering (headless)
        result = EegFun.plot_topography_3d(dat, display_plot = false)

        @test result.fig isa Figure
        @test first(result.axes) isa LScene
        @test result.axes == [result.ax]

        # Test colorrange
        res_cr = EegFun.plot_topography_3d(dat; colorrange = (-3.0, 3.0), display_plot = false)
        @test res_cr.fig isa Figure

        # Test vector of SingleDataFrameEeg (multi-condition 3D topoplot)
        dat2 = copy(dat)
        dat2.condition = 2
        dat2.condition_name = "condition_2"
        res_vec = EegFun.plot_topography_3d([dat, dat2]; display_plot = false)
        @test res_vec.fig isa Figure
        @test length(res_vec.axes) == 2
        @test all(a -> a isa LScene, res_vec.axes)
    end

    @testset "plot_topography_3d with kwargs" begin
        # Test applying label and point styles
        result = EegFun.plot_topography_3d(
            dat,
            display_plot = false,
            label_plot = true,
            label_fontsize = 16,
            label_color = :red,
            point_plot = true,
            point_markersize = 15,
            point_color = :blue,
            colormap = :viridis,
        )

        @test result.fig isa Figure
        @test first(result.axes) isa LScene
    end

    @testset "plot_topography for AbstractVector{<:EpochData}" begin
        ep1 = EegFun.create_test_epoch_data(n_channels = 6, n_epochs = 5)
        ep2 = EegFun.create_test_epoch_data(n_channels = 6, n_epochs = 5, condition = 2)
        ep1.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        ep1.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        ep2.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        ep2.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        vec_eps = [ep1, ep2]
        res = EegFun.plot_topography(@view(vec_eps[1:2]); display_plot = false)
        @test res.fig isa Figure
    end
end
