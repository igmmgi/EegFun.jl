using Test
using CairoMakie
using EegFun

@testset "TF Plots" begin

    # ────────────────────────────────────────────────────────────
    # Load real TF data once for all plotting tests
    # ────────────────────────────────────────────────────────────
    TF_TEST_DIR = joinpath(@__DIR__, "..", "..", "resources", "data", "julia", "tf")
    tf_participant1 = EegFun.read_data(joinpath(TF_TEST_DIR, "tf_data_1.jld2"))
    tf_participant2 = EegFun.read_data(joinpath(TF_TEST_DIR, "tf_data_2.jld2"))
    tf_all = vcat(tf_participant1, tf_participant2)

    # ────────────────────────────────────────────────────────────
    # 1. plot_topography (TF) — single TimeFreqData
    # ────────────────────────────────────────────────────────────
    @testset "plot_topography (single TF)" begin
        tf = tf_participant1[1]

        # Basic call — alpha band
        out = EegFun.plot_topography(tf, freq_range = (8.0, 13.0), display_plot = false)
        @test out.fig isa Figure

        # With interval selection
        out2 = EegFun.plot_topography(tf, freq_range = (4.0, 8.0), interval_selection = EegFun.times(0.0, 0.5), display_plot = false)
        @test out2.fig isa Figure
    end

    @testset "plot_topography (single TF) - baseline correction" begin
        tf = tf_participant1[1]

        # With baseline applied at plot time
        out = EegFun.plot_topography(
            tf,
            freq_range = (8.0, 13.0),
            interval_selection = EegFun.times(0.0, 0.5),
            baseline_interval = (-0.5, 0.0),
            baseline_method = :db,
            display_plot = false,
        )
        @test out.fig isa Figure

        # Pre-baselined data
        tf_bl = EegFun.tf_baseline(tf, (-0.5, 0.0); method = :db)
        out2 = EegFun.plot_topography(tf_bl, freq_range = (8.0, 13.0), display_plot = false)
        @test out2.fig isa Figure
    end

    # ────────────────────────────────────────────────────────────
    # 2. plot_topography (Vector{TimeFreqData})
    # ────────────────────────────────────────────────────────────
    @testset "plot_topography (multi TF)" begin
        # Two conditions
        out = EegFun.plot_topography(tf_participant1, freq_range = (8.0, 13.0), display_plot = false)
        @test out.fig isa Figure

        # Single element vector
        out2 = EegFun.plot_topography([tf_participant1[1]], freq_range = (8.0, 13.0), display_plot = false)
        @test out2.fig isa Figure
    end

    @testset "plot_topography (multi TF) - baseline" begin
        out = EegFun.plot_topography(
            tf_participant1,
            freq_range = (4.0, 13.0),
            baseline_interval = (-0.5, 0.0),
            baseline_method = :percent,
            display_plot = false,
        )
        @test out.fig isa Figure
    end

    @testset "plot_topography (TF) - empty/invalid" begin
        @test_throws Exception EegFun.plot_topography(EegFun.TimeFreqData[], freq_range = (4.0, 12.0))
    end

    # ────────────────────────────────────────────────────────────
    # 3. plot_topography_stats (TF)
    # ────────────────────────────────────────────────────────────
    # Prepare stats once for all topo_stats tests
    prepared_tf_stats = EegFun.prepare_stats(tf_all; design = :paired)
    tf_analytic_result = EegFun.analytic_test(prepared_tf_stats)

    @testset "plot_topography_stats - analytic result" begin
        out = EegFun.plot_topography_stats(tf_analytic_result, freq_range = (8.0, 13.0), n_topos = 4, display_plot = false)
        @test out.fig isa Figure
        @test length(out.axes) == 4
    end

    @testset "plot_topography_stats - data modes" begin
        for mode in (:tvalues, :difference)
            out = EegFun.plot_topography_stats(
                tf_analytic_result,
                freq_range = (8.0, 13.0),
                n_topos = 3,
                topo_data = mode,
                display_plot = false,
            )
            @test out.fig isa Figure
        end

        # Invalid mode
        @test_throws Exception EegFun.plot_topography_stats(
            tf_analytic_result,
            freq_range = (8.0, 13.0),
            topo_data = :invalid,
            display_plot = false,
        )
    end

    @testset "plot_topography_stats - invalid freq range" begin
        @test_throws Exception EegFun.plot_topography_stats(tf_analytic_result, freq_range = (100.0, 200.0), display_plot = false)
    end

    # ────────────────────────────────────────────────────────────
    # 4. plot_tf — single TimeFreqData & Vector{TimeFreqData}
    # ────────────────────────────────────────────────────────────
    @testset "plot_tf (single TF)" begin
        tf = tf_participant1[1]

        # Single layout
        out_single = EegFun.plot_tf(tf; layout = :single, display_plot = false)
        @test out_single.fig isa Figure
        @test length(out_single.axes) == 1

        # Grid layout with channel subset
        out_grid = EegFun.plot_tf(tf; channel_selection = EegFun.channels(:Oz, :Pz, :Cz), layout = :grid, display_plot = false)
        @test out_grid.fig isa Figure
        @test length(out_grid.axes) == 3

        # Baseline and custom colorrange
        out_bl = EegFun.plot_tf(
            tf;
            channel_selection = EegFun.channels(:Oz),
            baseline_interval = (-0.5, 0.0),
            baseline_method = :db,
            colorrange = (-2.0, 2.0),
            display_plot = false,
        )
        @test out_bl.fig isa Figure
    end

    @testset "plot_tf (multi TF)" begin
        # Vector of conditions
        out_multi = EegFun.plot_tf(tf_participant1; channel_selection = EegFun.channels(:Oz), display_plot = false)
        @test out_multi.fig isa Figure
        @test length(out_multi.axes) == length(tf_participant1)

        # Multi with baseline
        out_multi_bl = EegFun.plot_tf(
            tf_participant1;
            channel_selection = EegFun.channels(:Oz),
            baseline_interval = (-0.5, 0.0),
            baseline_method = :percent,
            display_plot = false,
        )
        @test out_multi_bl.fig isa Figure
    end

    @testset "plot_tf - invalid" begin
        @test_throws Exception EegFun.plot_tf(EegFun.TimeFreqData[])
    end

    @testset "plot_tf and topography ergonomics" begin
        tf = tf_participant1[1]
        # Direct channel symbol and string plot order
        out = EegFun.plot_tf(tf; channel_selection = :Oz, channel_plot_order = ["Oz"], display_plot = false)
        @test out.fig isa Figure

        # String channel selection in multi TF
        out_multi = EegFun.plot_tf(tf_participant1; channel_selection = "Oz", display_plot = false)
        @test out_multi.fig isa Figure

        # plot_topography with Vector freq_range and baseline_interval
        out_topo = EegFun.plot_topography(tf; freq_range = [8.0, 12.0], baseline_interval = [-0.5, 0.0], display_plot = false)
        @test out_topo.fig isa Figure

        # plot_topography for AbstractVector{<:TimeFreqData}
        out_topo_vec = EegFun.plot_topography(@view(tf_participant1[1:2]); freq_range = [8.0, 12.0], display_plot = false)
        @test out_topo_vec.fig isa Figure

        # Return shape includes .axes and .ax alias on single TF topography, plus normalized reversed freq_range
        out_rev = EegFun.plot_topography(
            tf;
            freq_range = (13.0, 8.0),
            baseline_interval = (-0.5, 0.0),
            baseline_method = :db,
            colormap = :coolwarm,
            display_plot = false,
        )
        @test out_rev.axes == [out_rev.ax]
        @test occursin("8-13 Hz", out_rev.ax.title[])
        @test !occursin("13-8", out_rev.ax.title[])
        cf = only(filter(p -> p isa Makie.Contourf, out_rev.ax.scene.plots))
        @test cf.colorrange[][1] ≈ -cf.colorrange[][2]
    end

    @testset "plot_tf_stats shared symmetric colorrange across channels" begin
        for mode in (:tvalues, :difference)
            out_stats = EegFun.plot_tf_stats(
                tf_analytic_result;
                content = mode,
                channel_selection = EegFun.channels(:Oz, :Pz, :Cz),
                layout = :grid,
                colormap = :coolwarm,
                display_plot = false,
            )
            @test out_stats.fig isa Figure
            @test length(out_stats.axes) == 3
            crs = [first(filter(p->p isa Makie.Heatmap, ax.scene.plots)).colorrange[] for ax in out_stats.axes]
            @test crs[1][1] ≈ -crs[1][2]
            @test all(cr -> cr == crs[1], crs)
        end

        # Non-significant channel under mask_style = :alpha is dimmed (alpha = 0.15), not drawn at 1.0
        tf_res_nosig = deepcopy(tf_analytic_result)
        fill!(tf_res_nosig.masks.positive, false)
        fill!(tf_res_nosig.masks.negative, false)
        out_nosig = EegFun.plot_tf_stats(
            tf_res_nosig;
            channel_selection = EegFun.channels(:Oz),
            mask_style = :alpha,
            display_plot = false,
        )
        hm_nosig = only(filter(p -> p isa Makie.Heatmap, out_nosig.axes[1].scene.plots))
        @test hm_nosig.alpha[] ≈ 0.15
    end

    @testset "plot_topography_stats (TF) subsets unanalyzed layout channels" begin
        tf_res_sub = deepcopy(tf_analytic_result)
        n_orig_layout = size(tf_res_sub.data[1].layout.data, 1)
        # Drop the last 2 electrodes from the stats result while keeping them in tf_res_sub.data[1].layout
        n_keep = length(tf_res_sub.electrodes) - 2
        tf_res_sub = EegFun.TFAnalyticResult(
            tf_res_sub.test_info,
            tf_res_sub.data,
            EegFun.TFStatMatrix(
                tf_res_sub.stat_matrix.t[1:n_keep, :, :],
                isnothing(tf_res_sub.stat_matrix.p) ? nothing : tf_res_sub.stat_matrix.p[1:n_keep, :, :],
            ),
            EegFun.TFMasks(tf_res_sub.masks.positive[1:n_keep, :, :], tf_res_sub.masks.negative[1:n_keep, :, :]),
            tf_res_sub.electrodes[1:n_keep],
            tf_res_sub.frequencies,
            tf_res_sub.time_points,
            tf_res_sub.critical_t,
        )
        out_sub = EegFun.plot_topography_stats(tf_res_sub; freq_range = (13.0, 8.0), n_topos = 2, point_plot = true, display_plot = false)
        @test out_sub.fig isa Figure
        @test length(out_sub.axes) == 2
        lbls = filter(b -> b isa Label, out_sub.fig.content)
        @test any(l -> occursin("8–13 Hz", l.text[]), lbls)
        # Point scatter on each topo axis should only have n_keep points, NOT n_orig_layout
        scatters = filter(p -> p isa Makie.Scatter, out_sub.axes[1].scene.plots)
        @test any(sc -> length(sc[1][]) == n_keep, scatters)
        @test !any(sc -> length(sc[1][]) == n_orig_layout, scatters)
    end

end # @testset "TF Plots"
