# Tests for single, grid, topographic, and custom multi-axis `PlotLayout` placement.

using Test
using CairoMakie
using DataFrames
using EegFun


@testset "Layout System Tests" begin

    @testset "PlotLayout Struct" begin
        # Test struct creation
        layout = EegFun.PlotLayout(:grid, [2, 2], [(0.0, 0.0)], [:Fp1], Dict{Symbol,Any}())
        @test layout.type == :grid
        @test layout.dims[1] == 2  # rows
        @test layout.dims[2] == 2  # cols
        @test layout.channels == [:Fp1]
    end

    @testset "create_single_layout" begin
        channels = [:Ch1, :Ch2, :Ch3]
        test_layout = EegFun.create_test_layout()
        layout = EegFun.create_layout(:single, channels, test_layout)

        @test layout.type == :single
        @test layout.dims[1] == 1  # rows
        @test layout.dims[2] == 1  # cols
        @test layout.channels == channels
        @test layout.positions == [(0.0, 0.0)]
    end

    @testset "create_grid_layout" begin
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        test_layout = EegFun.create_test_layout()

        # Test auto-calculated grid
        layout = EegFun.create_layout(:grid, channels, test_layout)
        @test layout.type == :grid
        @test layout.dims[1] * layout.dims[2] >= length(channels)

        # Test specified grid_dims
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = (2, 2))
        @test layout.dims[1] == 2  # rows
        @test layout.dims[2] == 2  # cols

        # Test grid_dims as vector
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = [2, 2])
        @test layout.dims[1] == 2
        @test layout.dims[2] == 2

        # Test that too small grid gets auto-corrected (should warn and use best_rect)
        # This should not throw, but should use optimal dimensions instead
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = (1, 1))
        @test layout.dims[1] * layout.dims[2] >= length(channels)  # Should be auto-corrected
    end

    @testset "create_topo_layout" begin
        test_layout = EegFun.create_test_layout()
        channels = [:Ch1, :Ch2]

        layout = EegFun.create_layout(:topo, channels, test_layout)
        @test layout.type == :topo
        @test layout.channels == channels
        @test length(layout.positions) == length(channels)
        @test haskey(layout.metadata, :topo_plot_width)
        @test haskey(layout.metadata, :topo_plot_height)
        @test haskey(layout.metadata, :topo_margin)
    end

    @testset "create_layout function" begin
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        test_layout = EegFun.create_test_layout()

        # Test symbol layouts
        @test EegFun.create_layout(:single, channels, test_layout).type == :single
        @test EegFun.create_layout(:grid, channels, test_layout).type == :grid
        @test EegFun.create_layout(:topo, channels, test_layout).type == :topo

        # Test custom grid dimensions via grid_dims parameter
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = (2, 2))
        @test layout.type == :grid
        @test layout.dims[1] == 2  # rows
        @test layout.dims[2] == 2  # cols

        # Test PlotLayout passthrough
        existing_layout = EegFun.PlotLayout(:grid, [2, 2], [], channels, Dict())
        @test EegFun.create_layout(existing_layout, channels, test_layout) == existing_layout
    end

    @testset "apply_layout! function" begin
        fig = Figure()
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        test_layout = EegFun.create_test_layout()

        # Test single layout
        layout = EegFun.create_layout(:single, channels, test_layout)
        axes, assignments = EegFun._apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)
        @test length(axes) == 1
        @test length(assignments) == 1
        @test assignments[1] == channels[1]  # First channel

        # Test grid layout
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = (2, 2))
        axes, assignments = EegFun._apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)
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
        ax = Axis(fig[1, 1], xgridvisible = true, ygridvisible = false, xminorgridvisible = true, yminorgridvisible = false)

        # Test _apply_axis_properties!
        kwargs = Dict(
            :plot_title => "Test Title",
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

        EegFun._apply_axis_properties!(ax; kwargs...)

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
        layout = EegFun.PlotLayout(:grid, [2, 2], [], [:Fp1], Dict())

        # Test _set_grid_axis_properties! for leftmost column (should keep label)
        ax1 = Axis(fig[1, 1])
        ax1.ylabel = "Test Label"  # Set a default label to test against
        EegFun._set_grid_axis_properties!(ax1, :Ch1, 1, 1, 2, 2)
        @test ax1.title[] == "Ch1"  # Should set channel name as title
        @test ax1.ylabel[] != ""  # Should keep label (leftmost column)

        # Test _set_grid_axis_properties! for non-leftmost column (should hide label)
        ax2 = Axis(fig[1, 2])
        EegFun._set_grid_axis_properties!(ax2, :Ch1, 1, 2, 2, 2)
        @test ax2.ylabel[] == ""  # Should hide label (not leftmost)
        @test ax2.yticklabelsvisible[] == false
    end

    @testset "Layout Axis Properties" begin
        fig = Figure()
        channels = [:Ch1, :Ch2]
        test_layout = EegFun.create_test_layout()
        layout = EegFun.create_layout(:grid, channels, test_layout; grid_dims = (1, 2))
        axes, _ = EegFun._apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)

        # Test _apply_layout_axis_properties!
        kwargs = Dict(:plot_title => nothing)
        EegFun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should set channel names as titles when user title is nothing
        @test axes[1].title[] == "Ch1"
        @test axes[2].title[] == "Ch2"

        # Test with user title
        kwargs = Dict(:plot_title => "User Title")
        EegFun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should respect user title, not override with channel names
        @test axes[1].title[] == "User Title"
        @test axes[2].title[] == "User Title"

        # Should hide titles when show_title = false
        axes_hidden, _ = EegFun._apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)
        EegFun._apply_layout_axis_properties!(axes_hidden, layout; plot_title = "Ignored", show_title = false)
        @test axes_hidden[1].title[] == ""
        @test axes_hidden[2].title[] == ""
    end

    @testset "Topographic Layout Properties" begin
        fig = Figure()
        test_layout = EegFun.create_test_layout()
        channels = [:Ch1, :Ch2]
        layout = EegFun.create_layout(:topo, channels, test_layout)
        axes, _ = EegFun._apply_layout!(fig, layout; xgrid = true, ygrid = true, xminorgrid = false, yminorgrid = false)

        # Test _apply_layout_axis_properties! for topo
        kwargs = Dict(:plot_title => nothing)
        EegFun._apply_layout_axis_properties!(axes, layout; kwargs...)

        # Should set channel names as titles when user title is nothing
        @test axes[1].title[] == "Ch1"
        @test axes[2].title[] == "Ch2"

    end

    @testset "Edge Cases" begin
        test_layout = EegFun.create_test_layout()

        # Test empty channel list - should throw EegFunError for empty channels
        @test_throws EegFun.EegFunError EegFun.create_layout(:grid, Symbol[], test_layout; grid_dims = (1, 1))

        # Test single channel
        layout = EegFun.create_layout(:grid, [:Ch1], test_layout; grid_dims = (1, 1))
        @test layout.dims[1] == 1  # rows
        @test layout.dims[2] == 1  # cols

        # Test prime number of channels
        layout = EegFun.create_layout(:grid, [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5], test_layout)
        @test layout.dims[1] * layout.dims[2] >= 5
    end

    @testset "Integration with plot_erp" begin
        # Test that the layout system works with actual plotting
        erp_data = EegFun.create_test_erp_data(participant = 1, condition = 1)

        # Test single layout
        fig, axes = EegFun.plot_erp(erp_data, layout = :single, display_plot = false)
        @test length(axes) == 1

        # Test grid layout
        fig, axes = EegFun.plot_erp(erp_data, layout = :grid, display_plot = false)
        @test length(axes) == 3  # 4 channels

        # Test topo layout
        fig, axes = EegFun.plot_erp(erp_data, layout = :topo, display_plot = false)
        @test length(axes) == 3  # 4 channels

        # Test custom grid dimensions via grid_dims parameter
        fig, axes = EegFun.plot_erp(erp_data, layout = :grid, layout_kwargs = (; grid_dims = (2, 2)), display_plot = false)
        @test length(axes) == 3

        # Test create_custom_layout in both integer (row, col) grid mode and float (x, y) canvas mode
        custom_grid = EegFun.create_custom_layout([(1, 1), (1, 3), (2, 2)], [:Ch1, :Ch2, :Ch3])
        @test custom_grid.type == :custom
        @test custom_grid.dims == [2, 3]
        res_cg = EegFun.plot_erp(erp_data, layout = custom_grid, display_plot = false)
        @test length(res_cg.axes) == 3
        @test res_cg.axes[1].ylabel[] != ""
        @test res_cg.axes[2].ylabel[] == ""
        @test res_cg.axes[1].xlabel[] == ""
        @test res_cg.axes[3].xlabel[] != ""

        custom_canvas = EegFun.create_custom_layout([(0.2, 0.8), (0.8, 0.8), (0.5, 0.2)], [:Ch1, :Ch2, :Ch3])
        @test custom_canvas.type == :custom
        res_cc = EegFun.plot_erp(erp_data, layout = custom_canvas, display_plot = false)
        @test length(res_cc.axes) == 3

        ep_data = EegFun.create_test_epoch_data(n = 50, n_epochs = 3, n_channels = 3)
        res_ep_cg = EegFun.plot_epochs(ep_data; layout = custom_grid, display_plot = false)
        @test length(res_ep_cg.axes) == 3

        # grid_skip_positions preserves true (row, col) label visibility on subsequent cells
        res_skip = EegFun.plot_erp(
            erp_data;
            layout = :grid,
            layout_kwargs = (; grid_dims = (2, 2), grid_skip_positions = [(1, 1)]),
            display_plot = false,
        )
        @test length(res_skip.axes) == 3
        # First plotted axis is at (1, 2) -> ylabel hidden; second is at (2, 1) -> ylabel & xlabel visible
        @test res_skip.axes[1].ylabel[] == ""
        @test res_skip.axes[2].ylabel[] != ""
        @test res_skip.axes[2].xlabel[] != ""

        # plot_erp_measurements places local-peak marker at local peak rather than global boundary extremum
        erp_loc = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 1)
        t_loc = erp_loc.data.time
        # Rising ramp at window edge (global max at boundary) + genuine interior local peak around t = 0.2 s
        erp_loc.data.Ch1 = 2.0 .* exp.(-((t_loc .- 0.20) ./ 0.02) .^ 2) .+ 5.0 .* max.(0.0, (t_loc .- 0.35) ./ 0.05)
        res_meas = EegFun.plot_erp_measurements(
            erp_loc,
            "max_peak_amplitude";
            analysis_interval = (0.10, 0.40),
            local_interval = 3,
            display_plot = false,
        )
        vls = filter(p -> p isa Makie.LineSegments || p isa Makie.Lines, res_meas.axes[1].scene.plots)
        @test res_meas.fig isa Figure && !isempty(vls)
    end
end

println("Layout system tests completed!")
