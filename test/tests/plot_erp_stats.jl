using Test
using Makie
using EegFun

@testset "plot_erp_stats" begin
    dat1 = EegFun.create_test_erp_data(participant = 1, n_channels = 3)
    dat2 = EegFun.create_test_erp_data(participant = 2, n_channels = 3)
    dat3 = EegFun.create_test_erp_data(participant = 1, condition = 2, n_channels = 3)
    dat4 = EegFun.create_test_erp_data(participant = 2, condition = 2, n_channels = 3)

    all_dat = [dat1, dat2, dat3, dat4]

    stat_data = EegFun.prepare_stats(all_dat; design = :paired)
    res = EegFun.analytic_test(stat_data)

    @testset "basic functionality" begin
        result = EegFun.plot_erp_stats(res; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp_stats(res; channel_selection = EegFun.channels([1]), display_plot = false)
        @test result.fig isa Figure

        # Symbol channel selection and String channel_plot_order
        res_sym = EegFun.plot_erp_stats(res; channel_selection = :Ch1, channel_plot_order = ["Ch1"], display_plot = false)
        @test res_sym.fig isa Figure

        # String channel selection and String channel_plot_order
        res_str = EegFun.plot_erp_stats(res; channel_selection = ["Ch1", "Ch2"], channel_plot_order = ["Ch2", "Ch1"], display_plot = false)
        @test res_str.fig isa Figure
    end
end
