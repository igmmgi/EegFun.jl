using Test
using DataFrames
using Statistics
using Makie
using eegfun

@testset "plot_channel_summary" begin

    # Helper function to create test data for plotting
    function create_test_summary_data()
        # Create a DataFrame with channel summary statistics
        df = DataFrame(
            channel = [:Fp1, :Fp2, :Cz, :Pz, :Oz],
            min = [-2.1, -1.8, -0.5, -0.3, -0.2],
            max = [2.3, 1.9, 0.8, 0.4, 0.3],
            std = [1.2, 1.1, 0.4, 0.2, 0.15],
            var = [1.44, 1.21, 0.16, 0.04, 0.0225],
            range = [4.4, 3.7, 1.3, 0.7, 0.5],
            zvar = [1.2, 0.8, -0.5, -1.1, -1.4]
        )
        return df
    end

    function create_test_summary_data_with_epochs()
        # Create a DataFrame with epoch information for testing averaging
        df = DataFrame(
            channel = repeat([:Fp1, :Fp2, :Cz], 3),
            epoch = repeat([1, 2, 3], 3),
            min = [-2.1, -1.8, -0.5, -2.0, -1.7, -0.4, -2.2, -1.9, -0.6],
            max = [2.3, 1.9, 0.8, 2.2, 1.8, 0.7, 2.4, 2.0, 0.9],
            std = [1.2, 1.1, 0.4, 1.1, 1.0, 0.3, 1.3, 1.2, 0.5],
            var = [1.44, 1.21, 0.16, 1.21, 1.0, 0.09, 1.69, 1.44, 0.25],
            range = [4.4, 3.7, 1.3, 4.2, 3.5, 1.1, 4.6, 3.9, 1.5],
            zvar = [1.2, 0.8, -0.5, 0.9, 0.6, -0.8, 1.5, 1.1, -0.2]
        )
        return df
    end

    @testset "plot_channel_summary! basic functionality" begin
        df = create_test_summary_data()
        
        # Test basic plotting
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        # Should not throw an error
        @test eegfun.plot_channel_summary!(fig, ax, df, :std) === nothing
        
        # Test that the axis has been modified
        @test ax.title[] == ""  # Default title
        @test ax.xlabel[] == "Electrode"  # Default xlabel
        @test ax.ylabel[] == "std"  # Should be set to the column name
    end

    @testset "plot_channel_summary! with custom kwargs" begin
        df = create_test_summary_data()
        
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        # Test with custom parameters
        @test eegfun.plot_channel_summary!(fig, ax, df, :range,
            title = "Custom Title",
            xlabel = "Custom X Label",
            bar_color = :red,
            sort_values = true
        ) === nothing
        
        # Test that custom parameters were applied
        @test ax.title[] == "Custom Title"
        @test ax.xlabel[] == "Custom X Label"
    end

    @testset "plot_channel_summary! with averaging" begin
        df = create_test_summary_data_with_epochs()
        
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        # Test with averaging over epochs
        @test eegfun.plot_channel_summary!(fig, ax, df, :std,
            average_over = :epoch,
            error_color = :blue,
            error_linewidth = 3
        ) === nothing
        
        # Should have error bars when averaging
        @test ax.ylabel[] == "std (± 95% CI n=3)"
    end

    @testset "plot_channel_summary! input validation" begin
        df = create_test_summary_data()
        
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        # Test missing channel column - should log error but not throw
        df_no_channel = select(df, Not(:channel))
        @test eegfun.plot_channel_summary!(fig, ax, df_no_channel, :std) === nothing
        
        # Test missing data column - should log error but not throw
        @test eegfun.plot_channel_summary!(fig, ax, df, :nonexistent) === nothing
        
        # Test invalid averaging column - should log error but not throw
        @test eegfun.plot_channel_summary!(fig, ax, df, :std, average_over = :nonexistent) === nothing
    end

    @testset "plot_channel_summary basic functionality" begin
        df = create_test_summary_data()
        
        # Test basic plotting (non-mutating version)
        fig, ax = eegfun.plot_channel_summary(df, :std, display_plot = false)
        
        @test fig isa Figure
        @test ax isa Axis
        @test ax.ylabel[] == "std"
    end

    @testset "plot_channel_summary with custom kwargs" begin
        df = create_test_summary_data()
        
        # Test with custom parameters
        fig, ax = eegfun.plot_channel_summary(df, :range,
            title = "Test Title",
            bar_color = :green,
            sort_values = true,
            display_plot = false  # Don't display during testing
        )
        
        @test fig isa Figure
        @test ax isa Axis
        @test ax.title[] == "Test Title"
        @test ax.ylabel[] == "range"
    end

    @testset "plot_channel_summary with averaging" begin
        df = create_test_summary_data_with_epochs()
        
        # Test with averaging
        fig, ax = eegfun.plot_channel_summary(df, :var,
            average_over = :epoch,
            error_color = :red,
            display_plot = false
        )
        
        @test fig isa Figure
        @test ax isa Axis
        @test ax.ylabel[] == "var (± 95% CI n=3)"
    end

    @testset "plot_channel_summary input validation" begin
        df = create_test_summary_data()
        
        # Test missing channel column - should log error but not throw
        df_no_channel = select(df, Not(:channel))
        @test_throws MethodError eegfun.plot_channel_summary(df_no_channel, :std)
        
        # Test missing data column - should log error but not throw  
        @test_throws MethodError eegfun.plot_channel_summary(df, :nonexistent)
    end

    @testset "plot_channel_summary sorting functionality" begin
        df = create_test_summary_data()
        
        # Test sorting by values
        fig, ax = eegfun.plot_channel_summary(df, :std, sort_values = true, display_plot = false)
        
        @test fig isa Figure
        @test ax isa Axis
        
        # Test without sorting
        fig2, ax2 = eegfun.plot_channel_summary(df, :std, sort_values = false, display_plot = false)
        
        @test fig2 isa Figure
        @test ax2 isa Axis
    end

    @testset "plot_channel_summary display control" begin
        df = create_test_summary_data()
        
        # Test with display_plot = false
        fig, ax = eegfun.plot_channel_summary(df, :std, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        
        # Test with display_plot = true (should throw due to no Makie backend in tests)
        @test_throws MethodError eegfun.plot_channel_summary(df, :std, display_plot = true)
    end

    @testset "plot_channel_summary different columns" begin
        df = create_test_summary_data()
        
        # Test plotting different statistical measures
        columns_to_test = [:min, :max, :std, :var, :range, :zvar]
        
        for col in columns_to_test
            fig, ax = eegfun.plot_channel_summary(df, col, display_plot = false)
            @test fig isa Figure
            @test ax isa Axis
            @test ax.ylabel[] == string(col)
        end
    end

    @testset "plot_channel_summary edge cases" begin
        # Test with single channel
        df_single = DataFrame(
            channel = [:Fp1],
            min = [-1.0],
            max = [1.0],
            std = [0.5],
            var = [0.25],
            range = [2.0],
            zvar = [0.0]
        )
        
        fig, ax = eegfun.plot_channel_summary(df_single, :std, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
        @test ax.ylabel[] == "std"
        
        # Test with two channels
        df_two = DataFrame(
            channel = [:Fp1, :Fp2],
            min = [-1.0, -0.8],
            max = [1.0, 0.9],
            std = [0.5, 0.4],
            var = [0.25, 0.16],
            range = [2.0, 1.7],
            zvar = [0.0, -0.5]
        )
        
        fig, ax = eegfun.plot_channel_summary(df_two, :var, display_plot = false)
        @test fig isa Figure
        @test ax isa Axis
    end

    @testset "plot_channel_summary consistency between versions" begin
        df = create_test_summary_data()
        
        # Test that both versions produce the same visual result
        fig1, ax1 = eegfun.plot_channel_summary(df, :std, display_plot = false)
        
        fig2 = Figure()
        ax2 = Axis(fig2[1, 1])
        eegfun.plot_channel_summary!(fig2, ax2, df, :std)
        
        # Both should have the same ylabel
        @test ax1.ylabel[] == ax2.ylabel[]
        @test ax1.title[] == ax2.title[]
        @test ax1.xlabel[] == ax2.xlabel[]
    end

    @testset "DEFAULT_CHANNEL_SUMMARY_KWARGS" begin
        # Test that the default kwargs are properly defined
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :sort_values)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :average_over)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :display_plot)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_color)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :title)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :xlabel)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_width)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :bar_alpha)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :error_color)
        @test haskey(eegfun.PLOT_CHANNEL_SUMMARY_KWARGS, :error_linewidth)
        
        # Test that defaults are reasonable
        @test eegfun.PLOT_CHANNEL_SUMMARY_KWARGS[:sort_values][1] == false
        @test eegfun.PLOT_CHANNEL_SUMMARY_KWARGS[:display_plot][1] == true
        @test eegfun.PLOT_CHANNEL_SUMMARY_KWARGS[:bar_color][1] == :steelblue
        @test eegfun.PLOT_CHANNEL_SUMMARY_KWARGS[:xlabel][1] == "Electrode"
    end

end # plot_channel_summary testset
