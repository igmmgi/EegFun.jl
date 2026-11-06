using Test
using EegFun
using CairoMakie
using DataFrames
using JLD2
using Random

@testset "Core utilities and headless plotting" begin
    # 1. Pipeline raw data inspection (src/pipelines/utils/utils.jl: check_raw_data)
    @testset "check_raw_data" begin
        bdf_file = EegFun.example_path("data/bdf/example1.bdf")
        hdr = EegFun.read_raw_data(bdf_file; header_only = true)
        n_chan = length(hdr.channel_labels)
        sr = Int(hdr.sample_rate[1])
        n_rec = 40
        short_buf = open(io -> read(io, 256 * (n_chan + 1) + n_rec * n_chan * sr * 3), bdf_file, "r")
        short_buf[237:244] .= codeunits(rpad(string(n_rec), 8))

        redirect_stdio(stdout = devnull, stderr = devnull) do
            mktempdir() do raw_dir
                short_path = joinpath(raw_dir, "short1.bdf")
                write(short_path, short_buf)
                res = EegFun.check_raw_data(raw_dir)
                @test res isa DataFrame
                @test nrow(res) == 1
                @test isfile(joinpath(raw_dir, "check_raw_data.log"))

                raw_bdf = EegFun.read_raw_data(short_path)
                bdf_cnt = EegFun.trigger_count(raw_bdf)
                @test occursin("Raw vs Cleaned", sprint(show, bdf_cnt))
                @test nrow(EegFun.all_data(bdf_cnt)) >= 1
            end

            mktempdir() do empty_dir
                @test_throws Exception EegFun.check_raw_data(empty_dir)
            end
        end
    end

    # 2. TriggerInfo show, trigger_count, batch trigger_info, subset & allcomb/epoch helpers
    @testset "TriggerInfo, subset_files, allcomb & epoch helpers" begin
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

        # mark_epoch_intervals! on Symbol and Vector{Symbol} Bool columns
        c_dat.data.ev_a = falses(nrow(c_dat.data))
        c_dat.data.ev_b = falses(nrow(c_dat.data))
        c_dat.data.ev_a[40] = true
        c_dat.data.ev_b[120] = true
        EegFun.mark_epoch_intervals!(c_dat, :ev_a, (-0.05, 0.05))
        @test any(c_dat.data.ev_a_window)
        EegFun.mark_epoch_intervals!(c_dat, [:ev_a, :ev_b], (-0.05, 0.05); channel_out = :combined_win)
        @test any(c_dat.data.combined_win)

        # condition_parse_epoch with allcomb(...) string expressions + epochs_count
        cfg_dict = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict("name" => "CombCond", "trigger_sequences" => "allcomb([1, 2], 10:11, :any, [[3, 4]], 99)"),
                    Dict("name" => "SimpleTrigger", "trigger_sequences" => [[1]]),
                ],
            ),
        )
        conds = EegFun.condition_parse_epoch(cfg_dict)
        @test length(conds) == 2
        @test length(conds[1].trigger_sequences) == 4
        @test_throws Exception EegFun.condition_parse_epoch(
            Dict("epochs" => Dict("conditions" => [Dict("name" => "Bad", "trigger_sequences" => "not_allcomb")])),
        )
        @test_throws Exception EegFun.condition_parse_epoch(
            Dict("epochs" => Dict("conditions" => [Dict("name" => "Bad2", "trigger_sequences" => "allcomb()")])),
        )

        cnt_df = EegFun.epochs_count(c_dat, conds; print_table = false)
        @test cnt_df isa DataFrame
        @test cnt_df.n_epochs[2] == 2

        ep_vec = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 4, n_channels = 2, fs = 100, n = 40, seed = 15)
        for (ci, ed) in enumerate(ep_vec), (ei, df) in enumerate(ed.data)
            df[!, :trigger] = zeros(Int, nrow(df))
            df.trigger[argmin(abs.(df.time))] = ci
            df[!, :bad_a] = fill(ei == 1, nrow(df))
            df[!, :bad_b] = fill(ei == 2, nrow(df))
        end
        ti_vec = EegFun.trigger_info(ep_vec)
        @test length(ti_vec) == 2

        # epochs_table & reject_epochs / reject_epochs! with Vector{Symbol}
        et_df = EegFun.epochs_table(ep_vec; print_table = false)
        @test et_df.n_epochs == [4, 4]
        ep_clean = EegFun.reject_epochs(ep_vec, [:bad_a, :bad_b])
        @test EegFun.n_epochs(ep_clean[1]) == 2
        ep_mut = deepcopy(ep_vec)
        EegFun.reject_epochs!(ep_mut, [:bad_a, :bad_b])
        @test EegFun.n_epochs(ep_mut[1]) == 2

        mktempdir() do tmp
            toml_path = joinpath(tmp, "conds.toml")
            write(toml_path, "[epochs]\n[[epochs.conditions]]\nname = \"C1\"\ntrigger_sequences = [[1]]\n")
            @test EegFun.epochs_count(c_dat, toml_path; print_table = false).n_epochs == [2]

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

    # 3. Config parameter documentation & data access helpers (config.jl, access.jl)
    @testset "show_parameter_info & data access helpers" begin
        redirect_stdio(stdout = devnull, stderr = devnull) do
            EegFun.show_parameter_info()
            EegFun.show_parameter_info(parameter_name = "preprocess")
            EegFun.show_parameter_info(parameter_name = "preprocess.ica.percentage_of_data")
            EegFun.show_parameter_info(parameter_name = "nonexistent.parameter")
        end

        ep_vec = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 3, n_channels = 2, fs = 100, n = 20, seed = 8)
        erps_v = EegFun.average_epochs(ep_vec)
        @test nrow(EegFun.all_data(DataFrame[])) == 0
        @test nrow(EegFun.all_data([DataFrame(a = [1, 2]), DataFrame(a = [3])])) == 3
        @test nrow(EegFun.all_data(erps_v)) == 40
        @test nrow(EegFun.all_data(ep_vec[1]; epoch_selection = EegFun.epochs(1:2))) == 40
        @test ncol(EegFun.meta_data(ep_vec[1])) >= 1
        @test ncol(EegFun.channel_data(ep_vec[1])) == 2
        for df in ep_vec[1].data
            df[!, :extra_col] = ones(20)
        end
        @test ncol(EegFun.extra_data(ep_vec[1]; epoch_selection = EegFun.epochs(1:1))) >= 1
        @test :Ch1 in EegFun.all_labels(ep_vec[1], 1)
        @test :Ch1 in EegFun.all_labels(ep_vec[1].data[1])
        @test length(EegFun.time_vector(ep_vec[1], 2)) == 20
        @test length(EegFun.time_vector(erps_v)) == 20
        @test EegFun.channel_labels(erps_v[1], [1:2]) == [:Ch1, :Ch2]
    end

    # 4. Design-aware block permutation matrix generator & TF nonparametric permutations (permutations.jl, clustering.jl)
    @testset "generate_permutation_matrix & TF nonparametric permutations" begin
        rng = MersenneTwister(123)
        df_bal = DataFrame(Subject = repeat(1:3, inner = 4), Item = repeat(1:4, outer = 3), Cond = repeat([1, 2], 6))

        pm_within = EegFun.generate_permutation_matrix(df_bal, :Subject; n_perms = 10, type = :within, rng = rng)
        @test size(pm_within) == (12, 10)
        @test sort(pm_within[:, 1]) == collect(1:12)
        @test sort(pm_within[1:4, 1]) == [1, 2, 3, 4]

        pm_sync = EegFun.generate_permutation_matrix(df_bal, :Item; sync_col = :Subject, n_perms = 10, type = :synchronized, rng = rng)
        @test size(pm_sync) == (12, 10)
        @test sort(pm_sync[:, 1]) == collect(1:12)

        pm_glob = EegFun.generate_permutation_matrix(df_bal, :Subject; n_perms = 10, type = :global, rng = rng)
        @test size(pm_glob) == (12, 10)
        @test sort(pm_glob[:, 1]) == collect(1:12)

        @test size(EegFun.generate_permutation_matrix(df_bal, :Subject; n_perms = 0), 2) == 0
        @test_throws Exception EegFun.generate_permutation_matrix(df_bal, :MissingCol; n_perms = 5)
        @test_throws Exception EegFun.generate_permutation_matrix(df_bal, :Item; n_perms = 5, type = :synchronized)
        @test_throws Exception EegFun.generate_permutation_matrix(
            df_bal[1:11, :],
            :Item;
            sync_col = :Subject,
            n_perms = 5,
            type = :synchronized,
        )
        @test_throws Exception EegFun.generate_permutation_matrix(df_bal, :Subject; n_perms = 5, type = :invalid_type)

        # TF nonparametric_common & nonparametric_individual streaming permutation tests
        tfs = EegFun.TimeFreqData[]
        for p = 1:5, (cond, off) in ((1, 5.0), (2, 1.0))
            push!(
                tfs,
                EegFun.create_test_tf_data(
                    file = "p$p",
                    condition = cond,
                    condition_name = "c$cond",
                    n_channels = 2,
                    frequencies = [8.0, 12.0],
                    time_points = [0.0, 0.1, 0.2],
                    power_offset = off,
                    noise = 0.2,
                ),
            )
        end
        EegFun.get_neighbors_xy!(tfs[1].layout, 1.5)
        for tf in tfs
            tf.layout = copy(tfs[1].layout)
        end
        tf_prep = EegFun.prepare_stats(tfs; design = :paired)
        for tm in (:nonparametric_common, :nonparametric_individual)
            tf_perm = EegFun.permutation_test(
                tf_prep;
                n_permutations = 10,
                threshold_method = tm,
                cluster_type = :spatiotemporal,
                min_num_neighbors = 1,
                show_progress = false,
                rng = MersenneTwister(33),
            )
            @test tf_perm isa EegFun.TFStatsResult
            @test EegFun.frequency_vector(tf_perm) == [8.0, 12.0]
        end
    end

    # 5. Direct LIBSVM / LIBLINEAR binary & multiclass classifier + erp_measurements! (liblinear.jl, erp_measurements.jl)
    @testset "libsvm_classifier & mutating erp_measurements!" begin
        rng = MersenneTwister(99)
        X_bin = vcat(randn(rng, 20, 4) .- 2.0, randn(rng, 20, 4) .+ 2.0)
        y_bin = vcat(fill(1, 20), fill(2, 20))
        preds_bin = EegFun.libsvm_classifier(X_bin, y_bin, X_bin; cost = 1.0)
        @test sum(preds_bin .== y_bin) / length(y_bin) > 0.9

        X_mc = vcat(randn(rng, 15, 3) .- 3.0, randn(rng, 15, 3), randn(rng, 15, 3) .+ 3.0)
        y_mc = vcat(fill(1, 15), fill(2, 15), fill(3, 15))
        preds_mc = EegFun.libsvm_classifier(X_mc, y_mc, X_mc; cost = 1.0)
        @test sum(preds_mc .== y_mc) / length(y_mc) > 0.85

        @test isempty(EegFun.libsvm_classifier(zeros(0, 3), Int[], zeros(2, 3)))
        @test_throws Exception EegFun.libsvm_classifier(X_bin, y_bin, X_bin; cost = -1.0)
        X_bad = copy(X_bin)
        X_bad[1, 1] = NaN
        @test_throws Exception EegFun.libsvm_classifier(X_bad, y_bin, X_bin)

        # Mutating erp_measurements! on ErpData, EpochData, and Vector{ErpData}
        erp1 = EegFun.create_test_erp_data(participant = 1, condition = 1)
        erp2 = EegFun.create_test_erp_data(participant = 1, condition = 2)
        nt1 = EegFun.erp_measurements!(copy(erp1), :mean_amplitude; analysis_interval = (0.05, 0.2), baseline_interval = (-0.1, 0.0))
        @test nt1 isa NamedTuple
        nt_vec =
            EegFun.erp_measurements!([copy(erp1), copy(erp2)], "mean_amplitude"; analysis_interval = (0.05, 0.2), average_channels = true)
        @test length(nt_vec) == 2
        @test isnothing(EegFun.erp_measurements!(EegFun.ErpData[], "mean_amplitude"; analysis_interval = (0.05, 0.2)))
        @test isnothing(
            EegFun.erp_measurements!(
                copy(erp1),
                "mean_amplitude";
                analysis_interval = (0.05, 0.2),
                channel_selection = EegFun.channels(Symbol[]),
            ),
        )

        ep_m = EegFun.create_test_epoch_data(n = 60, fs = 200, n_epochs = 3, n_channels = 2)
        ep_rows = EegFun.erp_measurements!(ep_m, "mean_amplitude"; analysis_interval = (0.05, 0.2), baseline_interval = (0.0, 0.05))
        @test length(ep_rows) == 3
    end

    # 6. Epoched ICA, subsampling & headless ICA artifact / topography plots (ica/api.jl, plot_ica_artifacts.jl, plot_ica.jl)
    @testset "Epoched ICA & headless ICA artifact component plots" begin
        ep_vec = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 4, n_channels = 3, fs = 100, n = 50, seed = 21)
        lay32 = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi32.csv"))
        EegFun.polar_to_cartesian_xy!(lay32)
        EegFun.polar_to_cartesian_xyz!(lay32)
        lay3 = EegFun.subset_layout(lay32; channel_selection = EegFun.channels(lay32.data.label[1:3]))
        for ed in ep_vec
            for (ci, old_ch) in enumerate([:Ch1, :Ch2, :Ch3])
                new_ch = lay3.data.label[ci]
                for df in ed.data
                    rename!(df, old_ch => new_ch)
                end
            end
            ed.layout = copy(lay3)
        end

        prms = EegFun.IcaPrms(max_iter = 10)
        ica_ep = EegFun.run_ica(ep_vec[1]; percentage_of_data = 60.0, params = prms, rng = MersenneTwister(3))
        @test ica_ep isa EegFun.InfoIca
        @test size(ica_ep.unmixing, 1) == 3
        @test nrow(EegFun.all_data(ica_ep)) == 3

        ica_vec = EegFun.run_ica(ep_vec; remove_duplicates = true, percentage_of_data = 80.0, params = prms, rng = MersenneTwister(4))
        @test ica_vec isa EegFun.InfoIca

        metrics_df = DataFrame(
            Component = [1, 2, 3],
            vEOG_zscore = [3.5, 0.2, -0.4],
            hEOG_zscore = [0.1, -3.2, 0.3],
            z_spatial_kurtosis = [0.5, 0.2, 4.1],
            power_ratio_zscore = [0.2, 3.8, -0.1],
            harmonic_ratio = [1.0, 6.0, 1.2],
            mean_ibi_s = [0.85, NaN, 0.90],
            std_ibi_s = [0.04, NaN, 0.25],
            peak_ratio = [0.92, 0.1, 0.4],
        )
        @test EegFun.plot_eog_component_features(Dict(:vEOG => [1], :hEOG => [2]), metrics_df; display_plot = false).fig isa Makie.Figure
        @test EegFun.plot_eog_component_features(Dict(:vEOG => Int[], :hEOG => Int[]), DataFrame(); display_plot = false).fig isa
              Makie.Figure
        @test EegFun.plot_spatial_kurtosis_components([3], metrics_df; display_plot = false).fig isa Makie.Figure
        @test EegFun.plot_line_noise_components([2], metrics_df; legend = true, display_plot = false).fig isa Makie.Figure
        @test EegFun.plot_ecg_component_features([1], metrics_df; legend = true, display_plot = false).fig isa Makie.Figure

        art_comps = EegFun.ArtifactComponents(Dict(:vEOG => [1], :hEOG => [2]), [1], [2], [3], Int[])
        @test EegFun.plot_artifact_components(ica_ep, art_comps; gridscale = 12, display_plot = false).fig isa Makie.Figure
        empty_art = EegFun.ArtifactComponents(Dict(:vEOG => Int[], :hEOG => Int[]), Int[], Int[], Int[], Int[])
        @test EegFun.plot_artifact_components(ica_ep, empty_art; display_plot = false).fig isa Makie.Figure

        @test EegFun.plot_topography(
            ica_ep;
            component_selection = EegFun.components(1:2),
            use_global_scale = true,
            gridscale = 12,
            display_plot = false,
        ).fig isa Makie.Figure
    end

    # 7. Time-Frequency baseline methods, batch TF baseline & non-mutating dispatches (time_frequency/baseline.jl)
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

        # Non-standard / shuffled TF row order fallback in baseline
        tf_shuf = copy(tf)
        tf_shuf.data_power = tf_shuf.data_power[reverse(1:nrow(tf_shuf.data_power)), :]
        for m in (:db, :zscore)
            @test !isnothing(EegFun.tf_baseline(tf_shuf, b_win; method = m).baseline)
        end

        tf_copy = copy(tf)
        EegFun.baseline!(tf_copy; method = :db)
        @test !isnothing(tf_copy.baseline)
        @test_throws Exception EegFun.tf_baseline!(tf_copy, b_win; method = :db)
        @test_throws Exception EegFun.tf_baseline(tf, b_win; method = :not_a_method)

        tf_vec_copy = [copy(tf)]
        EegFun.baseline!(tf_vec_copy; method = :percent)
        @test !isnothing(tf_vec_copy[1].baseline)

        tf_trials = EegFun.tf_morlet(ep; frequencies = [10.0, 20.0], cycles = 3, return_trials = true)
        tf_trials_b = EegFun.tf_baseline(tf_trials; method = :percent)
        @test tf_trials_b isa EegFun.TimeFreqEpochData

        mktempdir() do tmp
            jldsave(joinpath(tmp, "1_tf.jld2"); data = [copy(tf)])
            redirect_stdio(stdout = devnull, stderr = devnull) do
                EegFun.tf_baseline("tf", b_win; input_dir = tmp, method = :db)
                @test isfile(joinpath(tmp, "tf_baseline_tf", "1_tf.jld2"))
            end
        end
    end

    # 8. Headless ERP errorbar, stat waveforms/heatmap, spherical-spline topography & ERP measurement plots
    @testset "Headless ERP errorbar, stat waveforms/heatmap, topo & measurement plots" begin
        lay = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi32.csv"))
        EegFun.polar_to_cartesian_xy!(lay)
        EegFun.polar_to_cartesian_xyz!(lay)
        chs = lay.data.label[1:6]
        lay6 = EegFun.subset_layout(lay; channel_selection = EegFun.channels(chs))

        t_vec = collect(-0.1:0.02:0.3)
        n_t = length(t_vec)
        erps = EegFun.ErpData[]
        for p = 1:4, (cond, amp) in ((1, 5.0), (2, 0.5))
            df = DataFrame(
                time = t_vec,
                sample = 1:n_t,
                condition = fill(cond, n_t),
                condition_name = fill("Cond$cond", n_t),
                participant = fill(p, n_t),
            )
            rng_p = MersenneTwister(50 * p + cond)
            for (ci, ch) in enumerate(chs)
                df[!, ch] = (amp + 0.3 * ci) .* exp.(-((t_vec .- 0.12) ./ 0.05) .^ 2) .+ 0.15 .* randn(rng_p, n_t)
            end
            push!(erps, EegFun.ErpData("p$p", cond, "Cond$cond", df, copy(lay6), 50, EegFun.AnalysisInfo(), 20))
        end

        r_eb = EegFun.plot_erp_errorbar(erps; layout = :single, channel_selection = EegFun.channels(chs[1]), display_plot = false)
        @test r_eb.fig isa Makie.Figure
        r_eb_grid = EegFun.plot_erp_errorbar(
            erps;
            layout = :grid,
            error_bars = :within_sem,
            average_channels = true,
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

        # plot_erp_stats (difference wave, offset, SE bands, t-values, critical-t, significance positions) & plot_stat_heatmap
        prep = EegFun.prepare_stats(erps; design = :paired)
        stat_res = EegFun.analytic_test(prep)
        @test EegFun.plot_erp_stats(
            stat_res;
            channel_selection = EegFun.channels(chs[1:2]),
            layout = :grid,
            plot_erp = true,
            plot_difference = true,
            difference_offset = -3.0,
            plot_se = true,
            plot_significance = true,
            significance_position = :bottom,
            display_plot = false,
        ).fig isa Makie.Figure
        @test EegFun.plot_erp_stats(
            stat_res;
            channel_selection = EegFun.channels(chs[1]),
            plot_erp = false,
            plot_difference = false,
            plot_tvalues = true,
            plot_critical_t = true,
            plot_significance = true,
            significance_position = :zero,
            display_plot = false,
        ).fig isa Makie.Figure

        hm_res = EegFun.plot_stat_heatmap(stat_res; figure_title = "Analytic T-Map", mask_style = :hide, display_plot = false)
        @test hm_res.fig isa Makie.Figure

        perm_res = EegFun.permutation_test(prep; n_permutations = 10, show_progress = false, rng = MersenneTwister(11))
        hm_perm = EegFun.plot_stat_heatmap(perm_res; time_unit = :ms, display_plot = false)
        @test hm_perm.fig isa Makie.Figure

        # Spherical spline & alternative interpolation methods + highlight_channels in plot_topography
        hl_group = [(channels = chs[1:2], color = :yellow, size = 10, label = true)]
        for interp in (:spherical_spline, :nearest, :shepard, :multiquadratic)
            rt = EegFun.plot_topography(
                erps[1];
                interval_selection = EegFun.times(0.08, 0.16),
                method = interp,
                gridscale = 14,
                highlight_channels = hl_group,
                display_plot = false,
            )
            @test rt.fig isa Makie.Figure
        end

        r_scrub = EegFun.plot_topography(erps[1]; scrubber = true, gridscale = 12, display_plot = false, interactive = false)
        @test r_scrub.fig isa Makie.Figure

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

            r_erp_file = EegFun.plot_erp(jld_path; channel_selection = EegFun.channels(chs[1]), display_plot = false)
            @test r_erp_file.fig isa Makie.Figure

            r_topo_file = EegFun.plot_topography(
                "erps";
                input_dir = tmp,
                interval_selection = EegFun.times(0.08, 0.16),
                gridscale = 12,
                display_plot = false,
            )
            @test r_topo_file[1].fig isa Makie.Figure
        end
    end

    # 9. Tables interface, repair_info, show methods, and channel_metrics edge cases
    @testset "Tables, repair_info, show methods, and channel_metrics" begin
        erp_s = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 12)
        EegFun.polar_to_cartesian_xyz!(erp_s.layout)
        EegFun.get_neighbors_xyz!(erp_s.layout, 1.2)
        ep_m = EegFun.create_test_epoch_data(n = 40, n_epochs = 3, n_channels = 4)
        tf_s = EegFun.create_test_tf_data(n_channels = 2)

        # Tables.jl + DataFrames trait/forwarding methods (src/utils/tables.jl)
        @test EegFun.Tables.istable(typeof(erp_s))
        @test EegFun.Tables.columnaccess(typeof(erp_s))
        @test EegFun.Tables.rowaccess(typeof(erp_s))
        @test nrow(erp_s) == nrow(erp_s.data)
        @test ncol(erp_s) == ncol(erp_s.data)
        @test names(erp_s) == names(erp_s.data)
        @test describe(erp_s) isa DataFrame
        @test nrow(ep_m) == 3 * 40
        @test ncol(ep_m) == ncol(ep_m.data[1])
        @test names(ep_m) == names(ep_m.data[1])
        @test describe(ep_m) isa DataFrame
        @test length(first(ep_m, 2)) == 2
        @test length(last(ep_m, 2)) == 2
        @test DataFrame(tf_s) isa DataFrame
        @test EegFun.Tables.istable(Vector{EegFun.ErpData})
        @test EegFun.Tables.columnaccess(Vector{EegFun.ErpData})
        @test EegFun.Tables.rowaccess(Vector{EegFun.ErpData})
        @test !isempty(collect(EegFun.Tables.rows([erp_s])))
        @test !isnothing(EegFun.Tables.schema([erp_s]))

        # Rich show methods (src/types/core.jl & src/pipelines/utils/repair_info.jl)
        meas_res = EegFun.erp_measurements(erp_s, "mean_amplitude"; analysis_interval = (0.05, 0.2), baseline_interval = (-0.1, 0.0))
        @test contains(sprint(show, meas_res), "ErpMeasurementsResult")

        cond_rich = EegFun.EpochCondition(
            name = "RichCond",
            trigger_sequences = [[1, 2], [1, 3]],
            reference_index = 2,
            timing_pairs = [(1, 2)],
            min_interval = 0.1,
            max_interval = 1.5,
            mask_before_trigger = 99,
            mask_after_trigger = 98,
        )
        @test contains(sprint(show, cond_rich), "Timing pairs")
        @test contains(sprint(show, erp_s.layout), "showing first 3 and last 3")
        nb_plain = Dict{Symbol,EegFun.Neighbors}(pairs(erp_s.layout.neighbors))
        @test contains(sprint(show, MIME("text/plain"), nb_plain), "showing first 3 and last 3")

        # EpochRepairInfo, ChannelRepairInfo, ArtifactInfo show with populated data
        EegFun.polar_to_cartesian_xyz!(ep_m.layout)
        EegFun.get_neighbors_xyz!(ep_m.layout, 1.5)
        ep_rep = EegFun.create_epoch_repair_info(ep_m, Dict(1 => [:Ch1]), Dict(2 => [:Ch2]), :neighbor_interpolation, ep_m.layout.neighbors)
        @test contains(sprint(show, ep_rep), "Neighbors used")
        cont_rep = EegFun.ContinuousRepairInfo("c_rep", :neighbor_interpolation, [:Ch1], [:Ch2])
        ch_rep_info = EegFun.ChannelRepairInfo(cont_rep, [ep_rep])
        @test contains(sprint(show, ch_rep_info), "Continuous level")
        art_comp = EegFun.ArtifactComponents(Dict(:vEOG => [1]), Int[], Int[], Int[], Int[])
        rej_item = EegFun.detect_bad_epochs_automatic(ep_m; abs_criterion = 1e6)
        art_full = EegFun.ArtifactInfo([cont_rep], [rej_item], art_comp)
        @test contains(sprint(show, art_full), "ICA components removed")

        # Non-contiguous sample selection & Float32 columns in channel_metrics
        c_f32 = EegFun.create_test_continuous_data(n = 60, fs = 100, n_channels = 3)
        c_f32.data[!, :Ch1] = Float32.(c_f32.data.Ch1)
        c_f32.data[!, :Ch3] .= 5.0 # constant channel -> ss == 0 branch
        non_contig_sel = x -> [i in (2, 5, 9, 14, 20, 27) for i = 1:nrow(x)]
        cm_nc = EegFun.correlation_matrix(c_f32; sample_selection = non_contig_sel)
        @test nrow(cm_nc) == 3
        cm_contig_f32 = EegFun.correlation_matrix(c_f32)
        @test nrow(cm_contig_f32) == 3
        X_nc = EegFun._extract_channel_sample_matrix(c_f32.data, [2, 5, 9, 14], [:Ch1, :Ch2])
        @test size(X_nc) == (4, 2)
        X_contig = EegFun._extract_channel_sample_matrix(c_f32.data, collect(1:10), [:Ch1, :Ch2])
        @test size(X_contig) == (10, 2)
    end
end
