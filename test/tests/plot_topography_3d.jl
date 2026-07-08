using Test
using DataFrames
using Makie
using EegFun

@testset "Plot Topography 3D" begin
    # Create test data with a layout that includes spherical coordinates
    dat = EegFun.create_test_erp_data(n_channels=6)
    # Give them varying coordinates so they don't snap to the same point
    dat.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    dat.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]

    @testset "plot_topography_3d basic functionality" begin
        # Test basic rendering (headless)
        result = EegFun.plot_topography_3d(dat, display_plot = false)
        
        @test result.fig isa Figure
        @test first(result.axes) isa LScene
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
            colormap = :viridis
        )
        
        @test result.fig isa Figure
        @test first(result.axes) isa LScene
    end
end
