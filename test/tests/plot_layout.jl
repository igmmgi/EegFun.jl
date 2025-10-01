@testset "Plot Layout Tests" begin
    using DataFrames
    using Makie: Figure, Axis, Axis3
    
    # Create test layout data
    function create_test_layout()
        df = DataFrame(
            label = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2],
            x2 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y2 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            x3 = [0.0, 0.0, 0.0, 0.0, -0.5, 0.5],
            y3 = [1.0, 1.0, 0.0, -1.0, -1.0, -1.0],
            z3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        return eegfun.Layout(df, nothing, nothing)
    end
    
    layout = create_test_layout()
    
    @testset "plot_layout_2d! basic functionality" begin
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        # Test basic functionality
        @test eegfun.plot_layout_2d!(fig, ax, layout, display_plot = false) === nothing
        
        # Test with neighbours
        @test eegfun.plot_layout_2d!(fig, ax, layout, neighbours = true, display_plot = false) === nothing
    end
    
    @testset "plot_layout_2d basic functionality" begin
        # Test basic functionality
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with neighbours
        fig, ax = eegfun.plot_layout_2d(layout, neighbours = true, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
    end
    
    @testset "plot_layout_2d with prefixed kwargs" begin
        fig, ax = eegfun.plot_layout_2d(layout, 
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
            label_yoffset = 2
        )
        @test fig isa Figure
        @test ax isa Axis
    end
    
    @testset "plot_layout_3d! basic functionality" begin
        fig = Figure()
        ax = Axis3(fig[1, 1])
        
        # Test basic functionality
        result = eegfun.plot_layout_3d!(fig, ax, layout, display_plot = false)
        @test result isa Tuple && length(result) == 2
        
        # Test with neighbours
        result = eegfun.plot_layout_3d!(fig, ax, layout, neighbours = true, display_plot = false)
        @test result isa Tuple && length(result) == 2
    end
    
    @testset "plot_layout_3d basic functionality" begin
        # Test basic functionality
        fig, ax = eegfun.plot_layout_3d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis3
        
        # Test with neighbours
        fig, ax = eegfun.plot_layout_3d(layout, neighbours = true, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis3
    end
    
    @testset "plot_layout_3d with prefixed kwargs" begin
        fig, ax = eegfun.plot_layout_3d(layout, 
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
            label_zoffset = 2
        )
        @test fig isa Figure
        @test ax isa Axis3
    end
    
    @testset "add_topo_rois! basic functionality" begin
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        
        # Test basic ROI
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], border_size = 10) === nothing
        
        # Test multiple ROIs
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2], [:O1, :O2]], border_size = 5) === nothing
    end
    
    @testset "add_topo_rois! with prefixed kwargs" begin
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        
        # Test with custom ROI styling
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], 
            border_size = 10,
            roi_color = :red,
            roi_linewidth = 3,
            roi_fill = true,
            roi_fillcolor = :blue,
            roi_fillalpha = 0.3
        ) === nothing
    end
    
    @testset "add_topo_rois! with multiple ROIs and different styles" begin
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        
        # Test multiple ROIs with different styles
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2], [:O1, :O2]], 
            border_size = 8,
            roi_color = [:red, :blue],
            roi_linewidth = [2, 4],
            roi_fill = [true, false],
            roi_fillcolor = [:green, :yellow],
            roi_fillalpha = [0.2, 0.5]
        ) === nothing
    end
    
    @testset "add_topo_rois! error handling" begin
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        
        # Test with non-existent electrodes (should warn but not fail)
        @test eegfun.add_topo_rois!(ax, layout.data, [[:NonExistent]], border_size = 10) === nothing
        
        # Test with mismatched array lengths (should throw error)
        @test_throws ArgumentError eegfun.add_topo_rois!(ax, layout.data, [[:Fp1], [:Fp2]], 
            border_size = 10,
            roi_color = [:red]  # Only one color for two ROIs
        )
    end
    
    @testset "kwargs validation" begin
        # Test that invalid parameters are handled gracefully
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            invalid_param = :test,  # Should be ignored
            head_color = :red       # Should work
        )
        @test fig isa Figure
        @test ax isa Axis
    end
    
    @testset "display_plot parameter" begin
        # Test display_plot = false (should not throw)
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        
        # Test display_plot = true (might throw in headless environment)
        try
            fig, ax = eegfun.plot_layout_2d(layout, display_plot = true)
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
        empty_layout = eegfun.Layout(empty_df, nothing, nothing)
        
        fig, ax = eegfun.plot_layout_2d(empty_layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with single electrode
        single_df = DataFrame(label = [:Cz], x2 = [0.0], y2 = [0.0], x3 = [0.0], y3 = [0.0], z3 = [0.0])
        single_layout = eegfun.Layout(single_df, nothing, nothing)
        
        fig, ax = eegfun.plot_layout_2d(single_layout, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with very large border_size
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], border_size = 1000) === nothing
        
        # Test with zero border_size
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], border_size = 0) === nothing
        
        # Test with negative border_size (should handle gracefully)
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], border_size = -5) === nothing
    end
    
    @testset "component control parameters" begin
        # Test with all components disabled
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            point_plot = false,
            label_plot = false
        )
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with only points, no labels
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            point_plot = true,
            label_plot = false
        )
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with only labels, no points
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            point_plot = false,
            label_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
    end
    
    @testset "extreme parameter values" begin
        # Test with extreme font sizes
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            label_fontsize = 1,  # Very small
            label_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
        
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            label_fontsize = 100,  # Very large
            label_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with extreme marker sizes
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            point_markersize = 1,  # Very small
            point_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
        
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            point_markersize = 100,  # Very large
            point_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with extreme offsets
        fig, ax = eegfun.plot_layout_2d(layout, 
            display_plot = false,
            label_xoffset = 1000,
            label_yoffset = -1000,
            label_plot = true
        )
        @test fig isa Figure
        @test ax isa Axis
    end
    
    @testset "ROI edge cases" begin
        fig, ax = eegfun.plot_layout_2d(layout, display_plot = false)
        
        # Test with empty ROI list
        @test eegfun.add_topo_rois!(ax, layout.data, Vector{Symbol}[], border_size = 10) === nothing
        
        # Test with single electrode ROI
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Cz]], border_size = 10) === nothing
        
        # Test with all electrodes in one ROI
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]], border_size = 10) === nothing
        
        # Test with overlapping ROIs
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Cz], [:Cz, :Pz]], border_size = 10) === nothing
        
        # Test with extreme alpha values
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], 
            border_size = 10,
            roi_fill = true,
            roi_fillalpha = 0.0  # Completely transparent
        ) === nothing
        
        @test eegfun.add_topo_rois!(ax, layout.data, [[:Fp1, :Fp2]], 
            border_size = 10,
            roi_fill = true,
            roi_fillalpha = 1.0  # Completely opaque
        ) === nothing
    end
    
    @testset "constants and defaults" begin
        # Test that constants are properly defined
        @test haskey(eegfun.PLOT_LAYOUT_HEAD_KWARGS, :head_color)
        @test haskey(eegfun.PLOT_LAYOUT_HEAD_KWARGS, :head_linewidth)
        @test haskey(eegfun.PLOT_LAYOUT_POINT_KWARGS, :point_plot)
        @test haskey(eegfun.PLOT_LAYOUT_POINT_KWARGS, :point_color)
        @test haskey(eegfun.PLOT_LAYOUT_LABEL_KWARGS, :label_plot)
        @test haskey(eegfun.PLOT_LAYOUT_LABEL_KWARGS, :label_fontsize)
        @test haskey(eegfun.PLOT_LAYOUT_ROI_KWARGS, :roi_color)
        @test haskey(eegfun.PLOT_LAYOUT_ROI_KWARGS, :roi_fill)
        
        # Test default values
        @test eegfun.PLOT_LAYOUT_HEAD_KWARGS[:head_color][1] == :black
        @test eegfun.PLOT_LAYOUT_HEAD_KWARGS[:head_linewidth][1] == 2
        @test eegfun.PLOT_LAYOUT_POINT_KWARGS[:point_plot][1] == true
        @test eegfun.PLOT_LAYOUT_POINT_KWARGS[:point_color][1] == :black
        @test eegfun.PLOT_LAYOUT_LABEL_KWARGS[:label_plot][1] == true
        @test eegfun.PLOT_LAYOUT_LABEL_KWARGS[:label_fontsize][1] == 20
    end
end
