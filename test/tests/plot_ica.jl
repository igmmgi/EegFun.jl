using Test
using Makie
using EegFun

@testset "plot_ica" begin
    # Needs a real setup with ICA components
    dat = EegFun.create_test_continuous_data(n = 1000, n_channels = 4)
    dat.layout = EegFun.create_test_layout(n_channels = 4)

    # Run mock ICA (since real ICA takes a while, we test what we can)
    ica_res = EegFun.run_ica(dat)

    @testset "plot_ica_component_activation" begin
        result = EegFun.plot_ica_component_activation(dat, ica_res; display_plot = false)
        @test result.fig isa Figure
    end

    @testset "plot_ica_component_spectrum" begin
        result = EegFun.plot_ica_component_spectrum(dat, ica_res; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "plot_topography ICA" begin
        result = EegFun.plot_topography(ica_res; display_plot = false)
        @test result.fig isa Figure
    end
end
