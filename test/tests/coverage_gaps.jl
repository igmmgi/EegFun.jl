using Test
using EegFun
using CairoMakie
using DataFrames
using JLD2
using Random

@testset "Targeted Core & Headless Plotting Coverage" begin
    # 1. Pipeline raw data inspection (src/pipelines/utils/utils.jl: check_raw_data)
    @testset "check_raw_data" begin
        bdf_file = EegFun.example_path("data/bdf/example1.bdf")
        redirect_stdio(stdout = devnull, stderr = devnull) do
            mktempdir() do raw_dir
                cp(bdf_file, joinpath(raw_dir, "example1.bdf"))
                res = EegFun.check_raw_data(raw_dir)
                @test res isa DataFrame
                @test nrow(res) == 1
                @test isfile(joinpath(raw_dir, "check_raw_data.log"))
            end

            mktempdir() do empty_dir
                @test_throws Exception EegFun.check_raw_data(empty_dir)
            end
        end
    end

    # 2. TriggerInfo show, trigger_count & batch trigger_info / subset (triggers.jl, subset_files.jl)
    @testset "TriggerInfo show, trigger_count, batch trigger_info & subset" begin
        c_dat = EegFun.create_test_continuous_data(n = 200, fs = 200, n_channels = 2, seed = 10)
        c_dat.data.trigger .= 0
        empty_shown = sprint(show, EegFun.trigger_count(c_dat))
        @test occursin("No non-zero triggers", empty_shown)

        c_dat.data.trigger[20] = 1
        c_dat.data.trigger[60] = 2
        c_dat.data.trigger[120] = 1
        t_cnt = EegFun.trigger_count(c_dat)
        shown = sprint(show, t_cnt)
        @test occursin("1", shown)
        @test occursin("2", shown)

        raw_bdf = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
        bdf_cnt = EegFun.trigger_count(raw_bdf)
        @test occursin("Raw vs Cleaned", sprint(show, bdf_cnt))

        ep_vec = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 3, n_channels = 2, fs = 100, n = 40, seed = 15)
        for (ci, ed) in enumerate(ep_vec), df in ed.data
            df[!, :trigger] = zeros(Int, nrow(df))
            z_idx = argmin(abs.(df.time))
            df.trigger[z_idx] = ci
        end
        ti_vec = EegFun.trigger_info(ep_vec)
        @test length(ti_vec) == 2

        mktempdir() do tmp
            p1_path = joinpath(tmp, "1_epochs.jld2")
            jldsave(p1_path; data = ep_vec)
            redirect_stdio(stdout = devnull, stderr = devnull) do
                ti_direct = EegFun.trigger_info(p1_path)
                @test length(ti_direct) == 2
                ti_batch = EegFun.trigger_info("epochs"; input_dir = tmp)
                @test length(ti_batch) == 2

                # Batch file subsetting (src/utils/data/subset_files.jl)
                EegFun.subset("epochs"; input_dir = tmp, channel_selection = EegFun.channels([:Ch1]))
                sub_file = joinpath(tmp, "subset_epochs", "1_epochs.jld2")
                @test isfile(sub_file)
                sub_loaded = EegFun.read_data(sub_file)
                @test EegFun.channel_labels(sub_loaded[1]) == [:Ch1]
            end
        end
    end

    # 3. Direct LIBSVM / LIBLINEAR binary & multiclass classifier (src/analysis/decoding/liblinear.jl)
    @testset "libsvm_classifier binary, multiclass & edge cases" begin
        rng = MersenneTwister(99)
        # Binary classification
        X_bin = vcat(randn(rng, 20, 4) .- 2.0, randn(rng, 20, 4) .+ 2.0)
        y_bin = vcat(fill(1, 20), fill(2, 20))
        preds_bin = EegFun.libsvm_classifier(X_bin, y_bin, X_bin; cost = 1.0)
        @test sum(preds_bin .== y_bin) / length(y_bin) > 0.9

        # Multiclass (3-class) classification
        X_mc = vcat(randn(rng, 15, 3) .- 3.0, randn(rng, 15, 3), randn(rng, 15, 3) .+ 3.0)
        y_mc = vcat(fill(1, 15), fill(2, 15), fill(3, 15))
        preds_mc = EegFun.libsvm_classifier(X_mc, y_mc, X_mc; cost = 1.0)
        @test sum(preds_mc .== y_mc) / length(y_mc) > 0.85

        # Empty & validation branches
        @test isempty(EegFun.libsvm_classifier(zeros(0, 3), Int[], zeros(2, 3)))
        @test_throws Exception EegFun.libsvm_classifier(X_bin, y_bin, X_bin; cost = -1.0)
        X_bad = copy(X_bin)
        X_bad[1, 1] = NaN
        @test_throws Exception EegFun.libsvm_classifier(X_bad, y_bin, X_bin)
    end

    # 4. Epoched ICA & subsampling (src/analysis/processing/ica/api.jl)
    @testset "run_ica on EpochData and subsampling" begin
        ep_vec = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 4, n_channels = 3, fs = 100, n = 50, seed = 21)
        prms = EegFun.IcaPrms(max_iter = 10)
        ica_ep = EegFun.run_ica(ep_vec[1]; percentage_of_data = 60.0, params = prms, rng = MersenneTwister(3))
        @test ica_ep isa EegFun.InfoIca
        @test size(ica_ep.unmixing, 1) == 3

        # Duplicate sample warning branch when overlapping epochs share :samples column
        for ed in ep_vec, df in ed.data
            df[!, :samples] = collect(1:nrow(df))
        end
        ica_vec = EegFun.run_ica(ep_vec; remove_duplicates = true, percentage_of_data = 80.0, params = prms, rng = MersenneTwister(4))
        @test ica_vec isa EegFun.InfoIca
    end

    # 5. Time-Frequency baseline methods & non-mutating dispatches (src/analysis/time_frequency/baseline.jl)
    @testset "TF baseline methods & dispatches" begin
        ep = EegFun.create_test_epoch_data(n = 64, fs = 128, n_epochs = 3, n_channels = 2, seed = 5)
        tf = EegFun.tf_morlet(ep; frequencies = [10.0, 20.0], cycles = 3)
        t_min, t_max = extrema(tf.data_power.time)
        b_win = (t_min, 0.5 * (t_min + t_max))

        for m in (:db, :percent, :relchange, :zscore, :absolute)
            tf_b = EegFun.tf_baseline(tf, b_win; method = m)
            @test tf_b isa EegFun.TimeFreqData
            @test !isnothing(tf_b.baseline)
        end

        # Full-interval default baseline & alias
        tf_copy = copy(tf)
        EegFun.baseline!(tf_copy; method = :db)
        @test !isnothing(tf_copy.baseline)
        # Double-baseline guard throws
        @test_throws Exception EegFun.tf_baseline!(tf_copy, b_win; method = :db)
        # Invalid method throws
        @test_throws Exception EegFun.tf_baseline(tf, b_win; method = :not_a_method)

        # Single-trial TimeFreqEpochData baseline
        tf_trials = EegFun.tf_morlet(ep; frequencies = [10.0, 20.0], cycles = 3, return_trials = true)
        tf_trials_b = EegFun.tf_baseline(tf_trials, b_win; method = :percent)
        @test tf_trials_b isa EegFun.TimeFreqEpochData
    end

    # 6. Headless ERP errorbar, stat heatmap, spherical-spline topography & ERP measurement plots
    @testset "Headless ERP errorbar, stat heatmap, spherical-spline topo & measurement plots" begin
        lay = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi32.csv"))
        EegFun.polar_to_cartesian_xy!(lay)
        EegFun.polar_to_cartesian_xyz!(lay)
        chs = lay.data.label[1:6]
        lay6 = EegFun.subset_layout(lay; channel_selection = EegFun.channels(chs))

        t_vec = collect(-0.1:0.02:0.3)
        n_t = length(t_vec)
        erps = EegFun.ErpData[]
        for p = 1:4, (cond, amp) in ((1, 4.0), (2, 1.0))
            df = DataFrame(
                time = t_vec,
                sample = 1:n_t,
                condition = fill(cond, n_t),
                condition_name = fill("Cond$cond", n_t),
                participant = fill(p, n_t),
            )
            rng_p = MersenneTwister(50 * p + cond)
            for (ci, ch) in enumerate(chs)
                df[!, ch] = (amp + 0.3 * ci) .* exp.(-((t_vec .- 0.12) ./ 0.05) .^ 2) .+ 0.2 .* randn(rng_p, n_t)
            end
            push!(erps, EegFun.ErpData("p$p", cond, "Cond$cond", df, copy(lay6), 50, EegFun.AnalysisInfo(), 20))
        end

        # plot_erp_errorbar (:single, :grid, single ErpData, mutating plot_erp! / plot_erp_errorbar!)
        r_eb = EegFun.plot_erp_errorbar(erps; layout = :single, channel_selection = EegFun.channels(chs[1]), display_plot = false)
        @test r_eb.fig isa Makie.Figure
        r_eb_grid = EegFun.plot_erp_errorbar(
            erps;
            layout = :grid,
            error_bars = :within_sem,
            channel_selection = EegFun.channels(chs[1:2]),
            display_plot = false,
        )
        @test r_eb_grid.fig isa Makie.Figure
        r_eb_single = EegFun.plot_erp_errorbar(erps[1]; channel_selection = EegFun.channels(chs[1]), display_plot = false)
        @test r_eb_single.fig isa Makie.Figure

        fig_mut = Makie.Figure()
        ax_mut = Makie.Axis(fig_mut[1, 1])
        EegFun.plot_erp!(fig_mut, ax_mut, erps[1]; channel_selection = EegFun.channels(chs[1]))
        ax_mut2 = Makie.Axis(fig_mut[1, 2])
        EegFun.plot_erp_errorbar!(fig_mut, ax_mut2, erps; channel_selection = EegFun.channels(chs[1]))

        # plot_stat_heatmap & plot_stat_heatmap!
        prep = EegFun.prepare_stats(erps; design = :paired)
        stat_res = EegFun.analytic_test(prep)
        hm_res = EegFun.plot_stat_heatmap(stat_res; figure_title = "Analytic T-Map", display_plot = false)
        @test hm_res.fig isa Makie.Figure

        perm_res = EegFun.permutation_test(prep; n_permutations = 15, show_progress = false, rng = MersenneTwister(11))
        hm_perm = EegFun.plot_stat_heatmap(perm_res; time_unit = :ms, display_plot = false)
        @test hm_perm.fig isa Makie.Figure

        # Spherical spline & alternative interpolation methods + highlight_channels in plot_topography
        hl_group = [(channels = chs[1:2], color = :yellow, size = 10, label = true)]
        for interp in (:spherical_spline, :nearest, :shepard, :multiquadratic)
            rt = EegFun.plot_topography(
                erps[1];
                interval_selection = EegFun.times(0.08, 0.16),
                method = interp,
                gridscale = 24,
                highlight_channels = hl_group,
                display_plot = false,
            )
            @test rt.fig isa Makie.Figure
        end

        # Scrubber headless construction
        r_scrub = EegFun.plot_topography(erps[1]; scrubber = true, gridscale = 20, display_plot = false, interactive = false)
        @test r_scrub.fig isa Makie.Figure

        # plot_erp_measurements (single & average_channels = true & .jld2 file) + file-batch plot_erp_errorbar
        for mtype in ("mean_amplitude", "max_peak_latency")
            rm = EegFun.plot_erp_measurements(
                erps[1:2],
                mtype;
                analysis_interval = (0.05, 0.2),
                baseline_interval = (-0.1, 0.0),
                channel_selection = EegFun.channels(chs[1:2]),
                average_channels = (mtype == "max_peak_latency"),
                display_plot = false,
            )
            @test rm.fig isa Makie.Figure
        end

        mktempdir() do tmp
            jld_path = joinpath(tmp, "1_erps.jld2")
            jld_path2 = joinpath(tmp, "2_erps.jld2")
            jldsave(jld_path; data = erps[1:2])
            jldsave(jld_path2; data = erps[3:4])
            rm_file = EegFun.plot_erp_measurements(
                jld_path,
                "mean_amplitude";
                analysis_interval = (0.05, 0.2),
                channel_selection = EegFun.channels(chs[1]),
                display_plot = false,
            )
            @test rm_file.fig isa Makie.Figure

            reb_file = EegFun.plot_erp_errorbar("erps"; input_dir = tmp, channel_selection = EegFun.channels(chs[1]), display_plot = false)
            @test reb_file.fig isa Makie.Figure
        end
    end
end
