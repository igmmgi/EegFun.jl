using Test
using Makie
using EegFun

@testset "plot_rsa" begin
    dat = EegFun.create_test_epoch_data_vector(conditions = 1:3, n_epochs = 5, n_channels = 3)

    model = EegFun.create_rdm([1, 1, 2])

    res = EegFun.rsa(dat)
    res = EegFun.compare_models(res, [model])
    # Populate significant p_values so `:star5` peak scatter path is exercised
    res.p_values = fill(0.01, length(res.times), 1)

    @testset "basic functionality and return shape" begin
        result = EegFun.plot_rsa(res; display_plot = false)
        @test result isa NamedTuple
        @test result.fig isa Figure
        @test result.axes isa Vector{Axis}
        @test length(result.axes) == 1
        @test result.axes[1].title[] == "Model Correlations"
    end

    @testset "tuple xlim/ylim/colorrange and show_title" begin
        res_cor = EegFun.plot_rsa(
            res;
            plot_type = :correlations,
            xlim = (-0.1, 0.4),
            ylim = (-0.5, 1.0),
            show_title = false,
            display_plot = false,
        )
        @test res_cor.fig isa Figure
        @test res_cor.axes[1].title[] == ""

        res_rdm = EegFun.plot_rsa(res; plot_type = :rdm, colorrange = (0.0, 1.5), show_title = true, display_plot = false)
        @test res_rdm.fig isa Figure
        @test startswith(res_rdm.axes[1].title[], "RDM at ")

        res_tc = EegFun.plot_rdm_timecourse(res; xlim = (-0.1, 0.4), display_plot = false)
        @test res_tc.fig isa Figure
        @test res_tc.axes[1].title[] == "RDM Timecourse"
    end
end
