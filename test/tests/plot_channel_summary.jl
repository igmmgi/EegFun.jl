using Test
using DataFrames
using Statistics
using Makie

@testset "plot_channel_summary" begin

    @testset "plot_channel_summary! basic functionality" begin
        df = EegFun.create_test_summary_data()

        # Test basic plotting
        fig = Figure()
        ax = Axis(fig[1, 1])

        # Should not throw an error
        @test isnothing(EegFun.plot_channel_summary!(fig, ax, df, :std))

        # Test that the axis has been modified
        @test ax.title[] == ""  # Default title
        @test ax.xlabel[] == "Electrode"  # Default xlabel
        @test ax.ylabel[] == "std"  # Should be set to the column name
    end

    @testset "plot_channel_summary! with custom kwargs" begin
        df = EegFun.create_test_summary_data()

        fig = Figure()
        ax = Axis(fig[1, 1])

        # Test with custom parameters
        @test isnothing(
            EegFun.plot_channel_summary!(
                fig,
                ax,
                df,
                :range,
                plot_title = "Custom Title",
                xlabel = "Custom X Label",
                bar_color = :red,
                sort_values = true,
            ),
        )

        # Test that custom parameters were applied
        @test ax.title[] == "Custom Title"
        @test ax.xlabel[] == "Custom X Label"
    end

    @testset "plot_channel_summary! with averaging" begin
        df = EegFun.create_test_summary_data_with_epochs()

        fig = Figure()
        ax = Axis(fig[1, 1])

        # Test with averaging over epochs
        @test isnothing(EegFun.plot_channel_summary!(fig, ax, df, :std, average_over = :epoch, error_color = :blue, error_linewidth = 3))

        # Should have error bars when averaging
        @test ax.ylabel[] == "std (± 95% CI n=3)"
    end

    @testset "plot_channel_summary! input validation" begin
        df = EegFun.create_test_summary_data()

        fig = Figure()
        ax = Axis(fig[1, 1])

        # Test missing channel column - should throw now
        df_no_channel = select(df, Not(:channel))
        @test_throws Exception EegFun.plot_channel_summary!(fig, ax, df_no_channel, :std)

        # Test missing data column - should throw now
        @test_throws Exception EegFun.plot_channel_summary!(fig, ax, df, :nonexistent)

        # Test invalid averaging column - should throw now
        @test_throws Exception EegFun.plot_channel_summary!(fig, ax, df, :std, average_over = :nonexistent)
    end

    @testset "plot_channel_summary basic functionality" begin
        df = EegFun.create_test_summary_data()

        # Test basic plotting (non-mutating version)
        result = EegFun.plot_channel_summary(df, :std, display_plot = false)

        @test result.fig isa Figure
        @test first(result.axes) isa Axis
        @test first(result.axes).ylabel[] == "std"
    end

    @testset "plot_channel_summary with custom kwargs" begin
        df = EegFun.create_test_summary_data()

        # Test with custom parameters
        result = EegFun.plot_channel_summary(
            df,
            :range,
            plot_title = "Test Title",
            bar_color = :green,
            sort_values = true,
            display_plot = false,  # Don't display during testing
        )

        @test result.fig isa Figure
        @test first(result.axes) isa Axis
        @test first(result.axes).title[] == "Test Title"
        @test first(result.axes).ylabel[] == "range"
    end

    @testset "plot_channel_summary with averaging" begin
        df = EegFun.create_test_summary_data_with_epochs()

        # Test with averaging
        result = EegFun.plot_channel_summary(df, :var, average_over = :epoch, error_color = :red, display_plot = false)

        @test result.fig isa Figure
        @test first(result.axes) isa Axis
        @test first(result.axes).ylabel[] == "var (± 95% CI n=3)"
    end

    @testset "plot_channel_summary sorting functionality" begin
        df = EegFun.create_test_summary_data()

        # Test sorting by values
        result = EegFun.plot_channel_summary(df, :std, sort_values = true, display_plot = false)

        @test result.fig isa Figure
        @test first(result.axes) isa Axis

        # Test without sorting
        result2 = EegFun.plot_channel_summary(df, :std, sort_values = false, display_plot = false)

        @test result2.fig isa Figure
        @test first(result2.axes) isa Axis
    end

    @testset "plot_channel_summary display control" begin
        df = EegFun.create_test_summary_data()

        # Test with display_plot = false
        result = EegFun.plot_channel_summary(df, :std, display_plot = false)
        @test result.fig isa Figure
        @test first(result.axes) isa Axis

        # Test with display_plot = true (should work since GLMakie is available)
        result = EegFun.plot_channel_summary(df, :std, display_plot = true)
        @test result.fig isa Figure
        @test first(result.axes) isa Axis
    end

    @testset "plot_channel_summary different columns" begin
        df = EegFun.create_test_summary_data()

        # Test plotting different statistical measures
        columns_to_test = [:min, :max, :std, :var, :range, :zvar]

        for col in columns_to_test
            result = EegFun.plot_channel_summary(df, col, display_plot = false)
            @test result.fig isa Figure
            @test first(result.axes) isa Axis
            @test first(result.axes).ylabel[] == string(col)
        end
    end

    @testset "plot_channel_summary edge cases" begin
        # Test with single channel
        df_single = DataFrame(channel = [:Fp1], min = [-1.0], max = [1.0], std = [0.5], var = [0.25], range = [2.0], zvar = [0.0])

        result = EegFun.plot_channel_summary(df_single, :std, display_plot = false)
        @test result.fig isa Figure
        @test first(result.axes) isa Axis
        @test first(result.axes).ylabel[] == "std"

        # Test with two channels
        df_two = DataFrame(
            channel = [:Fp1, :Fp2],
            min = [-1.0, -0.8],
            max = [1.0, 0.9],
            std = [0.5, 0.4],
            var = [0.25, 0.16],
            range = [2.0, 1.7],
            zvar = [0.0, -0.5],
        )

        result = EegFun.plot_channel_summary(df_two, :var, display_plot = false)
        @test result.fig isa Figure
        @test first(result.axes) isa Axis
    end

    @testset "plot_channel_summary consistency between versions" begin
        df = EegFun.create_test_summary_data()

        # Test that both versions produce the same visual result
        result1 = EegFun.plot_channel_summary(df, :std, display_plot = false)

        fig2 = Figure()
        ax2 = Axis(fig2[1, 1])
        EegFun.plot_channel_summary!(fig2, ax2, df, :std)

        # Both should have the same ylabel
        @test first(result1.axes).ylabel[] == ax2.ylabel[]
        @test first(result1.axes).title[] == ax2.title[]
        @test first(result1.axes).xlabel[] == ax2.xlabel[]
    end

    @testset "DEFAULT_CHANNEL_SUMMARY_KWARGS" begin
        # Test that the default kwargs are properly defined
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :sort_values)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :average_over)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :display_plot)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_color)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :plot_title)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :xlabel)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_width)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_alpha)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :error_color)
        @test haskey(EegFun.PLOT_CHANNEL_SUMMARY_KWARGS, :error_linewidth)

        # Test that defaults are reasonable
        @test EegFun.PLOT_CHANNEL_SUMMARY_KWARGS[:sort_values][1] == false
        @test EegFun.PLOT_CHANNEL_SUMMARY_KWARGS[:display_plot][1] == true
        @test EegFun.PLOT_CHANNEL_SUMMARY_KWARGS[:bar_color][1] == :steelblue
        @test EegFun.PLOT_CHANNEL_SUMMARY_KWARGS[:xlabel][1] == "Electrode"
    end

end # plot_channel_summary testset
