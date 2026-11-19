using Test
using DataFrames
using Makie
using EegFun

@testset "plot_erp" begin
    dat = EegFun.create_test_erp_data(n_channels = 3)
    dat2 = EegFun.create_test_erp_data(n_channels = 3, condition = 2)

    @testset "single ErpData" begin
        result = EegFun.plot_erp(dat; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "vector of ErpData" begin
        result = EegFun.plot_erp([dat, dat2]; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp(dat; channel_selection = EegFun.channels([1, 2]), display_plot = false)
        @test result.fig isa Figure
    end

    @testset "baseline interval" begin
        result = EegFun.plot_erp(dat; baseline_interval = (-0.1, 0.0), display_plot = false)
        @test result.fig isa Figure
    end

    @testset "legend kwargs and nbanks" begin
        leg_kw = Dict(:framevisible => false)
        result = EegFun.plot_erp([dat, dat2]; legend_kwargs = leg_kw, legend_nbanks = 2, layout = :grid, display_plot = false)
        @test result.fig isa Figure
        @test haskey(leg_kw, :framevisible)
    end

    @testset "empty data check" begin
        @test_throws EegFun.EegFunError EegFun.plot_erp(dat; condition_selection = EegFun.conditions("nonexistent"), display_plot = false)
    end

    @testset "channel_plot_order and AbstractVector views" begin
        # String channel_plot_order
        res = EegFun.plot_erp(dat; channel_plot_order = ["Ch2", "Ch1"], display_plot = false)
        @test res.fig isa Figure

        # AbstractVector view
        vec_dats = [dat, dat2]
        res_view = EegFun.plot_erp(@view(vec_dats[1:2]); channel_plot_order = [:Ch1, :Ch2], display_plot = false)
        @test res_view.fig isa Figure
    end

    @testset "plot_erp! and plot_erp_errorbar! kwargs & average_channels" begin
        p_erps = EegFun.ErpData[]
        for p = 1:3, c = 1:2
            e = EegFun.create_test_erp_data(participant = p, condition = c, fs = 100, n_channels = 3, seed = 10 * p + c)
            # Pure between-participant intercept + deterministic condition effect -> zero within-subject variance
            for ch in (:Ch1, :Ch2, :Ch3)
                e.data[!, ch] .= 50.0 * p + (c == 1 ? 2.0 : -2.0)
            end
            push!(p_erps, e)
        end

        fig = Figure()
        ax1 = Axis(fig[1, 1])
        EegFun.plot_erp!(
            fig,
            ax1,
            [dat, dat2];
            condition_selection = EegFun.conditions([2]),
            interval_selection = EegFun.times(0.0, 0.5),
            sample_selection = EegFun.samples(),
            average_channels = true,
        )
        @test !isempty(ax1.scene.plots)

        ax2 = Axis(fig[1, 2])
        EegFun.plot_erp_errorbar!(
            fig,
            ax2,
            p_erps;
            condition_selection = EegFun.conditions([1, 2]),
            interval_selection = EegFun.times(0.0, 0.5),
            average_channels = true,
            error_bars = :within_sem,
        )
        @test any(p -> p isa Makie.Lines, ax2.scene.plots)
        @test any(p -> p isa Makie.Band, ax2.scene.plots)

        # In-memory Cousineau-Morey :within_sem removes pure participant intercepts via _participant_group_key
        _, vars_within = EegFun._compute_dynamic_grand_averages(p_erps, :within_sem)
        _, vars_sem = EegFun._compute_dynamic_grand_averages(p_erps, :sem)
        @test maximum(vars_within[1][1].Ch1) < 1e-20
        @test minimum(vars_sem[1][1].Ch1) > 1000.0

        # Multi-participant in-memory plot_erp splits into one figure per participant
        multi_figs = EegFun.plot_erp(p_erps; display_plot = false)
        @test multi_figs isa Vector{NamedTuple}
        @test length(multi_figs) == 3

        # plot_channel_spectrum excludes deleted channels still absent from DataFrame
        cont = EegFun.create_test_continuous_data(n = 256, fs = 256, n_channels = 3)
        select!(cont.data, Not(:Ch3))
        spec_fig = EegFun.plot_channel_spectrum(cont; display_plot = false)
        @test spec_fig.fig isa Figure
    end

    @testset "plot_erp_filter_gui ContinuousData time origin & relative .jld2 input_dir" begin
        cont_full = EegFun.create_test_continuous_data(n = 500, fs = 200, n_channels = 2)
        cont_sub = EegFun.subset(cont_full; interval_selection = EegFun.times(0.5, 1.5))
        gui_res = EegFun.plot_erp_filter_gui(cont_sub; display_plot = false)
        line_plots = filter(p -> p isa Makie.Lines, gui_res.axes[1].scene.plots)
        @test !isempty(line_plots)
        pts = line_plots[1][1][]
        @test first(pts)[1] ≈ first(cont_sub.data.time) atol = 1e-5
        @test last(pts)[1] ≈ last(cont_sub.data.time) atol = 1e-5

        mktempdir() do tmpdir
            EegFun.save_data(joinpath(tmpdir, "single_erp.jld2"), [dat, dat2])
            p_erps = [
                EegFun.create_test_erp_data(participant = p, condition = c, fs = 100, n_channels = 3, seed = 20 * p + c) for p = 1:3 for
                c = 1:2
            ]
            EegFun.save_data(joinpath(tmpdir, "group_erps.jld2"), p_erps)

            res_erp = EegFun.plot_erp("single_erp.jld2"; input_dir = tmpdir, display_plot = false)
            @test res_erp.fig isa Figure

            res_eb = EegFun.plot_erp_errorbar("group_erps.jld2"; input_dir = tmpdir, display_plot = false)
            @test res_eb.fig isa Figure

            res_gfp = EegFun.plot_gfp("single_erp.jld2"; input_dir = tmpdir, display_plot = false)
            @test res_gfp.fig isa Figure

            res_topo = EegFun.plot_topography("single_erp.jld2"; input_dir = tmpdir, display_plot = false)
            @test res_topo.fig isa Figure

            fake_gui = (filename = Observable("single_erp.jld2"), directory = Observable(tmpdir))
            @test EegFun._gui_filepath(fake_gui) == joinpath(tmpdir, "single_erp.jld2")
        end
    end
end
