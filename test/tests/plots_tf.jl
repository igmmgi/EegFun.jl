using Test
using CairoMakie

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

end # @testset "TF Plots"
