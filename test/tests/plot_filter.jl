using Test
using Makie
using EegFun

@testset "plot_filter" begin
    dat = EegFun.create_test_continuous_data(n_channels = 3)

    @testset "plot_filter_response" begin
        # Creates a filter and plots response
        filter_info = EegFun.create_highpass_filter(0.1, 1000.0)

        fig = EegFun.plot_filter_response(filter_info; display_plot = false)
        @test fig.fig isa Figure
    end

    @testset "plot_erp_filter_gui" begin
        # Zero-arg call (loads bundled example data)
        res0 = EegFun.plot_erp_filter_gui(display_plot = false)
        @test res0.fig isa Figure
        @test haskey(res0, :screen)

        # Single ErpData
        erp_data = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
        res1 = EegFun.plot_erp_filter_gui(erp_data[1]; channel = :Cz, display_plot = false)
        @test res1.fig isa Figure

        # Vector of ErpData
        res2 = EegFun.plot_erp_filter_gui(erp_data; channel = :Cz, display_plot = false)
        @test res2.fig isa Figure

        # Filepath string
        res3 = EegFun.plot_erp_filter_gui(EegFun.example_path("data/julia/erps/example1_erps.jld2"); channel = :Cz, display_plot = false)
        @test res3.fig isa Figure
    end
end
