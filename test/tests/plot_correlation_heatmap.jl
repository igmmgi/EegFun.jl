using Test
using Makie
using EegFun

@testset "plot_correlation_heatmap" begin
    using DataFrames
    using Statistics

    labels = ["Ch1", "Ch2", "Ch3"]
    matrix = cor(rand(100, 3))
    corr_df = DataFrame(matrix, labels)
    insertcols!(corr_df, 1, :row => labels)

    @testset "basic functionality" begin
        result = EegFun.plot_correlation_heatmap(corr_df; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        # Actually plot_correlation_heatmap doesn't accept channel_selection as per the kwargs list!
        result = EegFun.plot_correlation_heatmap(corr_df; mask_range = (-0.1, 0.1), display_plot = false)
        @test result.fig isa Figure
    end
end
