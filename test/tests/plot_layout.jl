using DataFrames
using Makie

@testset "Plot Layout Tests" begin

    layout = EegFun.create_test_layout(; n_channels = 6, layout_type = :topo)

    @testset "plot_layout_2d! basic functionality" begin
        fig = Figure()
        ax = Axis(fig[1, 1])

        # Test basic functionality
        @test EegFun.plot_layout_2d!(fig, ax, layout, display_plot = false) === nothing

        # Test with neighbours
        @test EegFun.plot_layout_2d!(fig, ax, layout, neighbours = true, display_plot = false) === nothing
    end

    @testset "plot_layout_2d basic functionality" begin
        # Test basic functionality
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test with neighbours
        fig, ax = EegFun.plot_layout_2d(layout, neighbours = true, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "plot_layout_2d with prefixed kwargs" begin
        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            head_color = :red,
            head_linewidth = 3,
            point_plot = true,
            point_color = :blue,
            point_marker = :square,
            point_markersize = 15,
            label_plot = true,
            label_color = :green,
            label_fontsize = 18,
            label_xoffset = 2,
            label_yoffset = 2,
        )
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "plot_layout_3d! basic functionality" begin
        fig = Figure()
        ax = Axis3(fig[1, 1])

        # Test basic functionality
        result = EegFun.plot_layout_3d!(fig, ax, layout, display_plot = false)
        @test result isa Tuple && length(result) == 2

        # Test with neighbours
        result = EegFun.plot_layout_3d!(fig, ax, layout, neighbours = true, display_plot = false)
        @test result isa Tuple && length(result) == 2
    end

    @testset "plot_layout_3d basic functionality" begin
        # Test basic functionality
        fig, ax = EegFun.plot_layout_3d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis3

        # Test with neighbours
        fig, ax = EegFun.plot_layout_3d(layout, neighbours = true, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis3
    end

    @testset "plot_layout_3d with prefixed kwargs" begin
        fig, ax = EegFun.plot_layout_3d(
            layout,
            display_plot = false,
            head_color = :red,
            head_linewidth = 3,
            point_plot = true,
            point_color = :blue,
            point_marker = :square,
            point_markersize = 15,
            label_plot = true,
            label_color = :green,
            label_fontsize = 18,
            label_xoffset = 2,
            label_yoffset = 2,
            label_zoffset = 2,
        )
        @test fig isa Figure
        @test ax isa Axis3
    end

    @testset "add_topo_rois! basic functionality" begin
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)

        # Test basic ROI
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2]], topo_border_size = 10) === nothing

        # Test multiple ROIs
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2], [:O1, :O2]], roi_border_size = 5) === nothing
    end

    @testset "add_topo_rois! with prefixed kwargs" begin
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)

        # Test with custom ROI styling
        @test EegFun.add_topo_rois!(
            ax,
            layout,
            [[:Fp1, :Fp2]],
            roi_border_size = 10,
            roi_linecolor = :red,
            roi_linewidth = 3,
            roi_fill = true,
            roi_fillcolor = :blue,
            roi_fillalpha = 0.3,
        ) === nothing
    end

    @testset "add_topo_rois! with multiple ROIs and different styles" begin
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)

        # Test multiple ROIs with different styles
        @test EegFun.add_topo_rois!(
            ax,
            layout,
            [[:Fp1, :Fp2], [:O1, :O2]],
            roi_border_size = 8,
            roi_color = [:red, :blue],
            roi_linewidth = [2, 4],
            roi_fill = [true, false],
            roi_fillcolor = [:green, :yellow],
            roi_fillalpha = [0.2, 0.5],
        ) === nothing
    end

    @testset "add_topo_rois! error handling" begin
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)

        # Test with non-existent electrodes (should warn but not fail)
        @test EegFun.add_topo_rois!(ax, layout, [[:NonExistent]], roi_border_size = 10) === nothing

        # Test with mismatched array lengths (should throw error)
        @test_throws ArgumentError EegFun.add_topo_rois!(
            ax,
            layout,
            [[:Fp1], [:Fp2]],
            roi_border_size = 10,
            roi_color = [:red],  # Only one color for two ROIs
        )
    end

    @testset "kwargs validation" begin
        # Test that invalid parameters are handled gracefully
        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            invalid_param = :test,  # Should be ignored
            head_color = :red,       # Should work
        )
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "display_plot parameter" begin
        # Test display_plot = false (should not throw)
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test display_plot = true (might throw in headless environment)
        try
            fig, ax = EegFun.plot_layout_2d(layout, display_plot = true)
            @test fig isa Figure
            @test ax isa Axis
        catch e
            # Expected in headless test environment
            @test e isa MethodError
        end
    end

    @testset "edge cases and boundary conditions" begin
        # Test with empty layout
        empty_df = DataFrame(label = Symbol[], x2 = Float64[], y2 = Float64[], x3 = Float64[], y3 = Float64[], z3 = Float64[])
        empty_layout = EegFun.Layout(empty_df, nothing, nothing)

        fig, ax = EegFun.plot_layout_2d(empty_layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test with single electrode
        single_df = DataFrame(label = [:Cz], x2 = [0.0], y2 = [0.0], x3 = [0.0], y3 = [0.0], z3 = [0.0])
        single_layout = EegFun.Layout(single_df, nothing, nothing)

        fig, ax = EegFun.plot_layout_2d(single_layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test with very large border_size
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2]], roi_border_size = 1000) === nothing

        # Test with zero border_size
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2]], roi_border_size = 0) === nothing

        # Test with negative border_size (should handle gracefully)
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2]], roi_border_size = -5) === nothing
    end

    @testset "component control parameters" begin
        # Test with all components disabled
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false, point_plot = false, label_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test with only points, no labels
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false, point_plot = true, label_plot = false)
        @test fig isa Figure
        @test ax isa Axis

        # Test with only labels, no points
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false, point_plot = false, label_plot = true)
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "extreme parameter values" begin
        # Test with extreme font sizes
        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            label_fontsize = 1,  # Very small
            label_plot = true,
        )
        @test fig isa Figure
        @test ax isa Axis

        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            label_fontsize = 100,  # Very large
            label_plot = true,
        )
        @test fig isa Figure
        @test ax isa Axis

        # Test with extreme marker sizes
        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            point_markersize = 1,  # Very small
            point_plot = true,
        )
        @test fig isa Figure
        @test ax isa Axis

        fig, ax = EegFun.plot_layout_2d(
            layout,
            display_plot = false,
            point_markersize = 100,  # Very large
            point_plot = true,
        )
        @test fig isa Figure
        @test ax isa Axis

        # Test with extreme offsets
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false, label_xoffset = 1000, label_yoffset = -1000, label_plot = true)
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "ROI edge cases" begin
        fig, ax = EegFun.plot_layout_2d(layout, display_plot = false)

        # Test with empty ROI list
        @test EegFun.add_topo_rois!(ax, layout, Vector{Symbol}[], roi_border_size = 10) === nothing

        # Test with single electrode ROI
        @test EegFun.add_topo_rois!(ax, layout, [[:Cz]], roi_border_size = 10) === nothing

        # Test with all electrodes in one ROI
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]], roi_border_size = 10) === nothing

        # Test with overlapping ROIs
        @test EegFun.add_topo_rois!(ax, layout, [[:Fp1, :Cz], [:Cz, :Pz]], roi_border_size = 10) === nothing

        # Test with extreme alpha values
        @test EegFun.add_topo_rois!(
            ax,
            layout,
            [[:Fp1, :Fp2]],
            roi_border_size = 10,
            roi_fill = true,
            roi_fillalpha = 0.0,  # Completely transparent
        ) === nothing

        @test EegFun.add_topo_rois!(
            ax,
            layout,
            [[:Fp1, :Fp2]],
            roi_border_size = 10,
            roi_fill = true,
            roi_fillalpha = 1.0,  # Completely opaque
        ) === nothing
    end

    @testset "constants and defaults" begin
        # Test that constants are properly defined
        @test haskey(EegFun.PLOT_LAYOUT_HEAD_KWARGS, :head_color)
        @test haskey(EegFun.PLOT_LAYOUT_HEAD_KWARGS, :head_linewidth)
        @test haskey(EegFun.PLOT_LAYOUT_POINT_KWARGS, :point_plot)
        @test haskey(EegFun.PLOT_LAYOUT_POINT_KWARGS, :point_color)
        @test haskey(EegFun.PLOT_LAYOUT_LABEL_KWARGS, :label_plot)
        @test haskey(EegFun.PLOT_LAYOUT_LABEL_KWARGS, :label_fontsize)
        @test haskey(EegFun.PLOT_LAYOUT_ROI_KWARGS, :roi_linecolor)
        @test haskey(EegFun.PLOT_LAYOUT_ROI_KWARGS, :roi_fill)

        # Test default values
        @test EegFun.PLOT_LAYOUT_HEAD_KWARGS[:head_color][1] == :black
        @test EegFun.PLOT_LAYOUT_HEAD_KWARGS[:head_linewidth][1] == 2
        @test EegFun.PLOT_LAYOUT_POINT_KWARGS[:point_plot][1] == true
        @test EegFun.PLOT_LAYOUT_POINT_KWARGS[:point_color][1] == :black
        @test EegFun.PLOT_LAYOUT_LABEL_KWARGS[:label_plot][1] == true
        @test EegFun.PLOT_LAYOUT_LABEL_KWARGS[:label_fontsize][1] == 20
    end
end
