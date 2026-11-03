using Test
using DataFrames
using Makie
using EegFun

@testset "plot_topography" begin
    erp1 = EegFun.create_test_erp_data(n_channels = 6, condition = 1)
    erp2 = EegFun.create_test_erp_data(n_channels = 6, condition = 2)
    for erp in (erp1, erp2)
        erp.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        erp.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        EegFun.polar_to_cartesian_xy!(erp.layout)
    end

    @testset "ErpData and Vector{ErpData} return shape & kwargs" begin
        res1 = EegFun.plot_topography(erp1; colorrange = (-5.0, 5.0), display_plot = false)
        @test res1 isa NamedTuple
        @test res1.fig isa Figure
        @test res1.axes isa Vector{Axis}
        @test length(res1.axes) == 1

        # Backwards-compatible ylim alias mapped to colorrange, and n_topos alias mapped to n_topo
        res_alias = EegFun.plot_topography(erp1; ylim = (-4.0, 4.0), n_topos = 3, display_plot = false)
        @test res_alias.fig isa Figure
        @test length(res_alias.axes) == 3

        res_vec = EegFun.plot_topography([erp1, erp2]; n_topo = 2, colorbar_plot = true, display_plot = false)
        @test res_vec.fig isa Figure
        @test length(res_vec.axes) == 4
    end

    @testset "EpochData and Vector{EpochData} forward kwargs" begin
        ep1 = EegFun.create_test_epoch_data(n_channels = 6, n_epochs = 4, condition = 1)
        ep2 = EegFun.create_test_epoch_data(n_channels = 6, n_epochs = 4, condition = 2)
        for ep in (ep1, ep2)
            ep.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
            ep.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
            EegFun.polar_to_cartesian_xy!(ep.layout)
        end

        res_ep = EegFun.plot_topography(ep1, 2; colorrange = (-3.0, 3.0), display_plot = false)
        @test res_ep.fig isa Figure
        @test length(res_ep.axes) == 1

        res_vec_ep = EegFun.plot_topography(
            [ep1, ep2],
            2;
            colormap = :viridis,
            colorrange = (-2.5, 2.5),
            colorbar_plot = false,
            display_plot = false,
        )
        @test length(res_vec_ep) == 2
        @test res_vec_ep[1].fig isa Figure
        @test length(res_vec_ep[1].axes) == 1
        cf = only(filter(p -> p isa Makie.Contourf, res_vec_ep[1].axes[1].scene.plots))
        @test cf.colorrange[] == (-2.5, 2.5)
    end

    @testset "_scale_topo_levels! updates Contourf levels and colorrange in sync" begin
        res = EegFun.plot_topography(erp1; colorrange = (-4.0, 4.0), num_levels = 10, point_plot = true, display_plot = false)
        EegFun._scale_topo_levels!(res.axes[1], 0.5)
        cf = only(filter(p -> p isa Makie.Contourf, res.axes[1].scene.plots))
        @test cf.colorrange[] == (-2.0, 2.0)
        @test first(cf.levels[]) ≈ -2.0
        @test last(cf.levels[]) ≈ 2.0
    end

    @testset "dark theme background mask resolution" begin
        Makie.with_theme(Makie.Theme(backgroundcolor = Makie.RGBf(0.1, 0.1, 0.12))) do
            res_dark = EegFun.plot_topography(erp1; display_plot = false)
            @test res_dark.fig isa Figure
            polys = filter(p -> p isa Makie.Poly, res_dark.axes[1].scene.plots)
            @test !isempty(polys)
            @test polys[end].color[] == Makie.to_color(Makie.RGBf(0.1, 0.1, 0.12))
        end
    end
end
