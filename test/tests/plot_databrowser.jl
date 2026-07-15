using Test
using Makie
using EegFun

@testset "plot_databrowser" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "basic functionality" begin
        fig = EegFun.plot_databrowser(dat; display_plot = false)
        @test fig.fig isa Figure
    end

    @testset "with ICA" begin
        ica_res = EegFun.run_ica(dat)
        fig = EegFun.plot_databrowser(dat, ica_res; display_plot = false)
        @test fig.fig isa Figure
    end
end
