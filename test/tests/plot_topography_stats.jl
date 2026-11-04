using Test
using DataFrames
using Makie
using EegFun

@testset "plot_topography_stats" begin
    dats = EegFun.ErpData[]
    for p = 1:3, c = 1:2
        erp = EegFun.create_test_erp_data(participant = p, condition = c, n_channels = 6)
        erp.layout.data.inc = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        erp.layout.data.azi = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]
        EegFun.polar_to_cartesian_xy!(erp.layout)
        push!(dats, erp)
    end

    stat_data = EegFun.prepare_stats(dats; design = :paired)
    analytic_res = EegFun.analytic_test(stat_data)
    perm_res = EegFun.permutation_test(stat_data; n_permutations = 10)

    @testset "AnalyticResult and PermutationResult return shapes & modes" begin
        for res_stat in (analytic_res, perm_res)
            for mode in (:tvalues, :difference)
                out = EegFun.plot_topography_stats(res_stat; topo_data = mode, n_topos = 3, colorbar_plot = true, display_plot = false)
                @test out isa NamedTuple
                @test out.fig isa Figure
                @test out.axes isa Vector{Axis}
                @test length(out.axes) == 3
                @test haskey(out, :colorbar)
            end
        end

        # colorrange
        out_cr = EegFun.plot_topography_stats(analytic_res; colorrange = (-3.0, 3.0), n_topos = 2, display_plot = false)
        @test out_cr.fig isa Figure
        cf = only(filter(p -> p isa Makie.Contourf, out_cr.axes[1].scene.plots))
        @test cf.colorrange[] == (-3.0, 3.0)
    end

    @testset "Subsets layout to result.electrodes instead of zero-filling unanalyzed channels" begin
        # Add an unanalyzed channel :ExtraCh to a copy of the layout on result.data
        extra_layout = copy(analytic_res.data[1].layout)
        push!(extra_layout.data, (label = :ExtraCh, inc = 15.0, azi = 30.0, x2 = 0.15, y2 = 0.25); cols = :subset)
        d1_extra = copy(analytic_res.data[1])
        d2_extra = copy(analytic_res.data[2])
        d1_extra.layout = extra_layout
        d2_extra.layout = copy(extra_layout)
        res_subset = EegFun.AnalyticResult(
            analytic_res.test_info,
            [d1_extra, d2_extra],
            analytic_res.stat_matrix,
            analytic_res.masks,
            analytic_res.electrodes,
            analytic_res.time_points,
            analytic_res.critical_t,
            analytic_res.se,
            analytic_res.se_cond1,
            analytic_res.se_cond2,
            analytic_res.se_diff,
        )
        out_sub = EegFun.plot_topography_stats(res_subset; n_topos = 2, point_plot = true, display_plot = false)
        @test out_sub.fig isa Figure
        scatters = filter(p -> p isa Makie.Scatter, out_sub.axes[1].scene.plots)
        # Only the 6 analyzed electrodes are plotted, NOT 7!
        @test any(sc -> length(sc[1][]) == 6, scatters)
        @test !any(sc -> length(sc[1][]) == 7, scatters)
    end
end
