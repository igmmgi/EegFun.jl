using Test
using eegfun
using GLMakie
using DataFrames

# Test data setup
# Use generic create_test_layout from test_utils.jl
# create_test_layout(; n_channels::Int = 4, layout_type::Symbol = :grid)

# Use generic create_test_erp_data from test_utils.jl
# create_test_erp_data(participant, condition, n_timepoints, n_channels)

@testset "Layout System Tests" begin

    @testset "PlotLayout Struct" begin
        # Test struct creation
        layout = eegfun.PlotLayout(:grid, 2, 2, [(0.0, 0.0)], [:Fp1], Dict{Symbol,Any}())
        @test layout.type == :grid
        @test layout.rows == 2
        @test layout.cols == 2
        @test layout.channels == [:Fp1]
    end

    @testset "create_single_layout" begin
        channels = [:Ch1, :Ch2, :Ch3]
        layout = eegfun.create_single_layout(channels)

        @test layout.type == :single
        @test layout.rows == 1
        @test layout.cols == 1
        @test layout.channels == channels
        @test layout.positions == [(0.0, 0.0)]
    end

    @testset "create_grid_layout" begin
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]

        # Test auto-calculated grid
        layout = eegfun.create_grid_layout(channels)
        @test layout.type == :grid
        @test layout.rows * layout.cols >= length(channels)

        # Test specified rows
        layout = eegfun.create_grid_layout(channels, rows = 2)
        @test layout.rows == 2
        @test layout.cols == 2  # Should calculate cols automatically

        # Test specified cols
        layout = eegfun.create_grid_layout(channels, cols = 2)
        @test layout.cols == 2
        @test layout.rows == 2  # Should calculate rows automatically

        # Test both specified
        layout = eegfun.create_grid_layout(channels, rows = 2, cols = 2)
        @test layout.rows == 2
        @test layout.cols == 2

        # Test error for too small grid
        @test_throws ArgumentError eegfun.create_grid_layout(channels, rows = 1, cols = 1)
    end

    @testset "create_topo_layout" begin
        test_layout = create_test_layout()
        channels = [:Ch1, :Ch2]

        layout = eegfun.create_topo_layout(test_layout, channels)
        @test layout.type == :topo
        @test layout.channels == channels
        @test length(layout.positions) == length(channels)
        @test haskey(layout.metadata, :plot_width)
        @test haskey(layout.metadata, :plot_height)
        @test haskey(layout.metadata, :margin)
    end

    @testset "create_layout function" begin
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        test_layout = create_test_layout()

        # Test symbol layouts
        @test eegfun.create_layout(:single, channels, test_layout).type == :single
        @test eegfun.create_layout(:grid, channels, test_layout).type == :grid
        @test eegfun.create_layout(:topo, channels, test_layout).type == :topo

        # Test custom grid dimensions
        layout = eegfun.create_layout([2, 2], channels, test_layout)
        @test layout.type == :grid
        @test layout.rows == 2
        @test layout.cols == 2

        # Test PlotLayout passthrough
        existing_layout = eegfun.PlotLayout(:grid, 2, 2, [], channels, Dict())
        @test eegfun.create_layout(existing_layout, channels, test_layout) == existing_layout
    end

    @testset "apply_layout! function" begin
        fig = Figure()
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]

        # Test single layout
        layout = eegfun.create_single_layout(channels)
        axes, assignments =
            eegfun.apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)
        @test length(axes) == 1
        @test length(assignments) == 1
        @test assignments[1] == channels[1]  # First channel

        # Test grid layout
        layout = eegfun.create_grid_layout(channels, rows = 2, cols = 2)
        axes, assignments =
            eegfun.apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)
        @test length(axes) == 4
        @test length(assignments) == 4
        @test all(ax -> ax isa Axis, axes)

        # Test that each axis has correct grid properties
        for ax in axes
            @test hasfield(typeof(ax), :xgridvisible)
            @test hasfield(typeof(ax), :ygridvisible)
        end
    end

    @testset "Axis Property Functions" begin
        fig = Figure()
        # Create axis with the grid properties we want to test
        ax = Axis(
            fig[1, 1],
            xgridvisible = true,
            ygridvisible = false,
            xminorgridvisible = true,
            yminorgridvisible = false,
        )

        # Test _apply_axis_properties!
        kwargs = Dict(
            :title => "Test Title",
            :xlabel => "X Label",
            :ylabel => "Y Label",
            :xlim => (0, 10),
            :ylim => (-5, 5),
            :yreversed => true,
            :xgrid => true,
            :ygrid => false,
            :xminorgrid => true,
            :yminorgrid => false,
        )

        eegfun._apply_axis_properties!(ax; kwargs...)

        @test ax.title[] == "Test Title"
        @test ax.xlabel[] == "X Label"
        @test ax.ylabel[] == "Y Label"
        @test ax.yreversed[] == true
        @test ax.xgridvisible[] == true
        @test ax.ygridvisible[] == false
        @test ax.xminorgridvisible[] == true
        @test ax.yminorgridvisible[] == false
    end

    @testset "Grid Axis Properties" begin
        fig = Figure()
        layout = eegfun.PlotLayout(:grid, 2, 2, [], [:Fp1], Dict())

        # Test _set_grid_axis_properties! for leftmost column (should keep label)
        ax1 = Axis(fig[1, 1])
        ax1.ylabel = "Test Label"  # Set a default label to test against
        eegfun._set_grid_axis_properties!(ax1, layout, :Ch1, 1, 1, 2, 2)
        @test ax1.title[] == "Ch1"  # Should set channel name as title
        @test ax1.ylabel[] != ""  # Should keep label (leftmost column)

        # Test _set_grid_axis_properties! for non-leftmost column (should hide label)
        ax2 = Axis(fig[1, 2])
        eegfun._set_grid_axis_properties!(ax2, layout, :Ch1, 1, 2, 2, 2)
        @test ax2.ylabel[] == ""  # Should hide label (not leftmost)
        @test ax2.yticklabelsvisible[] == false
    end

    @testset "Layout Axis Properties" begin
        fig = Figure()
        channels = [:Ch1, :Ch2]
        layout = eegfun.create_grid_layout(channels, rows = 1, cols = 2)
        axes, _ = eegfun.apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)

        # Test _apply_layout_axis_properties!
        kwargs = Dict(:title => "")
        eegfun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should set channel names as titles when user title is empty
        @test axes[1].title[] == "Ch1"
        @test axes[2].title[] == "Ch2"

        # Test with user title
        kwargs = Dict(:title => "User Title")
        eegfun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should respect user title, not override with channel names
        @test axes[1].title[] == "User Title"
        @test axes[2].title[] == "User Title"
    end

    @testset "Topographic Layout Properties" begin
        fig = Figure()
        test_layout = create_test_layout()
        channels = [:Ch1, :Ch2]
        layout = eegfun.create_topo_layout(test_layout, channels)
        axes, _ = eegfun.apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)

        # Test _apply_layout_axis_properties! for topo
        kwargs = Dict(:title => "")
        eegfun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should set channel names as titles when user title is empty
        @test axes[1].title[] == "Ch1"
        @test axes[2].title[] == "Ch2"

    end

    @testset "Edge Cases" begin
        # Test empty channel list
        @test_throws MethodError eegfun.create_grid_layout([], rows = 1, cols = 1)

        # Test single channel
        layout = eegfun.create_grid_layout([:Ch1], rows = 1, cols = 1)
        @test layout.rows == 1
        @test layout.cols == 1

        # Test prime number of channels
        layout = eegfun.create_grid_layout([:Ch1, :Ch2, :Ch3, :Ch4, :Ch5])
        @test layout.rows * layout.cols >= 5
    end

    @testset "Integration with plot_erp" begin
        # Test that the layout system works with actual plotting
        erp_data = create_test_erp_data(1, 1, 501, 4)

        # Test single layout
        fig, axes = eegfun.plot_erp(erp_data, layout = :single)
        @test length(axes) == 1

        # Test grid layout
        fig, axes = eegfun.plot_erp(erp_data, layout = :grid)
        @test length(axes) == 4  # 4 channels

        # Test topo layout
        fig, axes = eegfun.plot_erp(erp_data, layout = :topo)
        @test length(axes) == 4  # 4 channels

        # Test custom grid dimensions
        fig, axes = eegfun.plot_erp(erp_data, layout = [2, 2])
        @test length(axes) == 4
    end
end

println("Layout system tests completed!")
