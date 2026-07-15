using Test
using Makie
using EegFun

@testset "plot_decoding" begin
    dat = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 10, n_channels = 3)

    res = EegFun.decode_libsvm(dat, n_iterations = 2)

    @testset "basic functionality" begin
        result = EegFun.plot_decoding(res; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    stats_res = EegFun.test_against_chance([res, res])

    @testset "with stats" begin
        result = EegFun.plot_decoding(res, stats_res; display_plot = false)
        @test result isa Figure
    end
end
