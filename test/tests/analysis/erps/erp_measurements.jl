# Tests for ERP amplitude, peak latency, area, and fractional latency measurements.

using Test
using JLD2
using DataFrames
using CSV
using CairoMakie
using EegFun

@testset "Batch ERP Measurements" begin
    function create_known_signal_erp(
        participant::Int,
        condition::Int,
        signal_func::Function;
        fs::Int = 1000,
        t_start::Float64 = 0.0,
        t_end::Float64 = 1.0,
    )
        time = collect(range(t_start, t_end, length = Int((t_end - t_start) * fs) + 1))
        df = DataFrame(
            time = time,
            sample = 1:length(time),
            condition = fill(condition, length(time)),
            condition_name = fill("condition_$condition", length(time)),
            participant = fill(participant, length(time)),
            Ch1 = signal_func.(time),
        )
        layout = EegFun.Layout(DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]), nothing, nothing, nothing)
        return EegFun.ErpData("test", condition, "condition_$condition", df, layout, fs, EegFun.AnalysisInfo(), 1)
    end

    mktempdir() do test_dir
        for p = 1:3
            jldsave(
                joinpath(test_dir, "$(p)_erps_unrejected.jld2");
                data = [EegFun.create_test_erp_data(participant = p, condition = c) for c in (1, 2)],
            )
        end
        for p = 1:2
            jldsave(
                joinpath(test_dir, "$(p)_epochs_unrejected.jld2");
                data = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_channels = 3),
            )
        end

        @testset "CSV output, logging, auto-dir & all measurement types" begin
            out_dir = joinpath(test_dir, "measurements")
            res = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = out_dir,
            )
            @test res isa EegFun.ErpMeasurementsResult && nrow(res.data) == 6 && ncol(res.data) >= 5
            @test names(res.data)[1:4] == ["file", "participant", "condition", "condition_name"]
            csv_df = CSV.read(joinpath(out_dir, "erp_measurements.csv"), DataFrame)
            @test nrow(csv_df) == 6 && all(c -> c in names(csv_df), ("participant", "condition", "Ch1"))
            log_txt = read(joinpath(out_dir, "erp_measurements.log"), String)
            @test occursin("ERP measurement analysis started", log_txt) && occursin("Analysis complete", log_txt)

            # Custom output file & auto-generated directory
            cust_dir = joinpath(test_dir, "custom_out")
            EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = cust_dir,
                output_file = "custom_meas",
            )
            @test isfile(joinpath(cust_dir, "custom_meas.csv"))
            EegFun.erp_measurements("erps_unrejected", "mean_amplitude"; analysis_interval = (0.1, 0.2), input_dir = test_dir)
            @test isdir(joinpath(test_dir, "measurements_mean_amplitude"))

            # All measurement types on batch JLD2
            types_out = joinpath(test_dir, "types_out")
            for (atype, check_fn) in (
                ("max_peak_amplitude", v -> all(v .> 0)),
                ("min_peak_amplitude", v -> all(isfinite, v)),
                ("max_peak_latency", v -> all(0.1 .<= v .<= 0.2)),
                ("min_peak_latency", v -> all(0.1 .<= v .<= 0.2)),
                ("peak_to_peak_amplitude", v -> all(v .> 0)),
                ("peak_to_peak_latency", v -> all(0.0 .<= v .<= 0.1)),
                ("rectified_area", v -> all(v .>= 0)),
                ("integral", v -> all(isfinite, v)),
                ("positive_area", v -> all(v .>= 0)),
                ("negative_area", v -> all(v .>= 0)),
                ("fractional_area_latency", v -> all(0.1 .<= v .<= 0.2)),
                ("fractional_peak_latency", v -> all(0.1 .<= v .<= 0.2)),
            )
                r = EegFun.erp_measurements(
                    "erps_unrejected",
                    atype;
                    analysis_interval = (0.1, 0.2),
                    input_dir = test_dir,
                    output_dir = types_out,
                )
                @test r isa EegFun.ErpMeasurementsResult && nrow(r.data) == 6
                for ch in (:Ch1, :Ch2, :Ch3)
                    vals = r.data[!, ch]
                    @test all(isfinite, vals) && check_fn(vals)
                end
            end
        end

        @testset "Kwargs, EpochData, selection filters, baseline & error/edge cases" begin
            kw_out = joinpath(test_dir, "kw_out")
            @test EegFun.erp_measurements(
                "erps_unrejected",
                "max_peak_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = kw_out,
                local_interval = 5,
            ) isa EegFun.ErpMeasurementsResult
            @test EegFun.erp_measurements(
                "erps_unrejected",
                "fractional_area_latency";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = kw_out,
                fractional_area_fraction = 0.3,
            ) isa EegFun.ErpMeasurementsResult
            @test EegFun.erp_measurements(
                "erps_unrejected",
                "fractional_peak_latency";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = kw_out,
                fractional_peak_fraction = 0.7,
                fractional_peak_direction = :offset,
            ) isa EegFun.ErpMeasurementsResult

            # Invalid kwargs & error handling
            for (atype, bad_kw) in (
                ("max_peak_amplitude", (; local_interval = 0)),
                ("fractional_area_latency", (; fractional_area_fraction = 1.5)),
                ("fractional_peak_latency", (; fractional_peak_fraction = -0.1)),
                ("fractional_peak_latency", (; fractional_peak_direction = :invalid)),
            )
                @test_throws Exception EegFun.erp_measurements(
                    "erps_unrejected",
                    atype;
                    analysis_interval = (0.1, 0.2),
                    input_dir = test_dir,
                    output_dir = kw_out,
                    bad_kw...,
                )
            end
            @test_throws Exception EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = joinpath(test_dir, "nonexistent"),
            )
            @test_throws Exception EegFun.erp_measurements(
                "erps_unrejected",
                "invalid_type";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
            )
            @test_throws Exception EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.2, 0.1),
                input_dir = test_dir,
            )

            # EpochData batch processing
            res_ep = EegFun.erp_measurements(
                "epochs_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = joinpath(test_dir, "ep_out"),
            )
            @test res_ep isa EegFun.ErpMeasurementsResult && nrow(res_ep.data) == 40 && "epoch" in names(res_ep.data)

            # Participant, condition & channel selection + baseline correction
            res_p = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                participant_selection = EegFun.participants([1, 2]),
                output_dir = kw_out,
            )
            @test nrow(res_p.data) == 4 && Set(res_p.data.participant) == Set([1, 2])

            res_c = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                condition_selection = EegFun.conditions([1]),
                output_dir = kw_out,
            )
            @test nrow(res_c.data) == 3 && all(==(1), res_c.data.condition)

            res_ch = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                channel_selection = EegFun.channels([:Ch1, :Ch2]),
                output_dir = kw_out,
            )
            @test "Ch1" in names(res_ch.data) && "Ch2" in names(res_ch.data) && "Ch3" ∉ names(res_ch.data)

            res_bl = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = [0.1, 0.2],
                baseline_interval = [-0.2, 0.0],
                input_dir = test_dir,
                output_dir = kw_out,
            )
            res_nobl = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                input_dir = test_dir,
                output_dir = kw_out,
            )
            @test res_bl.data.Ch1 != res_nobl.data.Ch1

            # Edge cases returning nothing or gracefully skipping empty baseline
            @test isnothing(
                EegFun.erp_measurements(
                    "erps_unrejected",
                    "mean_amplitude";
                    analysis_interval = (2.1, 3.0),
                    input_dir = test_dir,
                    output_dir = kw_out,
                ),
            )
            @test isnothing(
                EegFun.erp_measurements(
                    "nonexistent_pattern",
                    "mean_amplitude";
                    analysis_interval = (0.1, 0.2),
                    input_dir = test_dir,
                    output_dir = kw_out,
                ),
            )
            @test isnothing(
                EegFun.erp_measurements(
                    "erps_unrejected",
                    "mean_amplitude";
                    analysis_interval = (0.1, 0.2),
                    channel_selection = EegFun.channels(Symbol[]),
                    input_dir = test_dir,
                    output_dir = kw_out,
                ),
            )
            @test isnothing(
                EegFun.erp_measurements(
                    "erps_unrejected",
                    "mean_amplitude";
                    analysis_interval = (0.1, 0.2),
                    participant_selection = EegFun.participants(999),
                    input_dir = test_dir,
                    output_dir = kw_out,
                ),
            )
            @test EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.2),
                baseline_interval = (10.0, 11.0),
                input_dir = test_dir,
                output_dir = kw_out,
            ) isa EegFun.ErpMeasurementsResult
            res_narrow = EegFun.erp_measurements(
                "erps_unrejected",
                "mean_amplitude";
                analysis_interval = (0.1, 0.1001),
                input_dir = test_dir,
                output_dir = kw_out,
            )
            @test isnothing(res_narrow) || res_narrow isa EegFun.ErpMeasurementsResult
        end
    end

    @testset "Known analytical signal validation" begin
        erp_const5 = create_known_signal_erp(1, 1, t -> 5.0)
        erp_const3 = create_known_signal_erp(1, 1, t -> 3.0)
        erp_neg2 = create_known_signal_erp(1, 1, t -> -2.0)
        erp_max_pk = create_known_signal_erp(1, 1, t -> -(t - 0.5)^2 + 1.0)
        erp_min_pk = create_known_signal_erp(1, 1, t -> (t - 0.5)^2 - 1.0)
        erp_sin = create_known_signal_erp(1, 1, t -> sin(2π * t))
        erp_mixed = create_known_signal_erp(1, 1, t -> t - 0.5)
        erp_robust = create_known_signal_erp(1, 1, t -> -2 * (t - 0.5)^2 + 1.0 + 0.1 * sin(20π * t))
        erp_ramp = create_known_signal_erp(1, 1, t -> t)
        erp_tri = create_known_signal_erp(1, 1, t -> t < 0.5 ? 2.0 * t : 1.0 - 2.0 * (t - 0.5))
        erp_bl5 = create_known_signal_erp(1, 1, t -> 5.0; t_start = -0.2, t_end = 1.0)

        for (erp, atype, kw, expected, tol) in (
            (erp_const5, "mean_amplitude", (;), 5.0, 1e-6),
            (erp_max_pk, "max_peak_amplitude", (;), 1.0, 0.01),
            (erp_max_pk, "max_peak_latency", (;), 0.5, 0.01),
            (erp_min_pk, "min_peak_amplitude", (;), -1.0, 0.01),
            (erp_min_pk, "min_peak_latency", (;), 0.5, 0.01),
            (erp_sin, "peak_to_peak_amplitude", (;), 2.0, 0.1),
            (erp_sin, "peak_to_peak_latency", (;), 0.5, 0.03),
            (erp_const3, "rectified_area", (;), 3.0, 0.03),
            (erp_const3, "integral", (;), 3.0, 0.03),
            (erp_const3, "positive_area", (;), 3.0, 0.03),
            (erp_const3, "negative_area", (;), 0.0, 0.01),
            (erp_neg2, "negative_area", (;), 2.0, 0.02),
            (erp_neg2, "rectified_area", (;), 2.0, 0.02),
            (erp_neg2, "integral", (;), -2.0, 0.02),
            (erp_neg2, "positive_area", (;), 0.0, 0.01),
            (erp_mixed, "integral", (;), 0.0, 0.01),
            (erp_mixed, "positive_area", (;), 0.125, 0.005),
            (erp_mixed, "negative_area", (;), 0.125, 0.005),
            (erp_mixed, "rectified_area", (;), 0.25, 0.005),
            (erp_robust, "max_peak_latency", (; local_interval = 5), 0.5, 0.03),
            (erp_ramp, "fractional_area_latency", (; fractional_area_fraction = 0.0), 0.0, 0.01),
            (erp_ramp, "fractional_area_latency", (; fractional_area_fraction = 0.25), 0.5, 0.03),
            (erp_ramp, "fractional_area_latency", (; fractional_area_fraction = 0.5), sqrt(0.5), 0.04),
            (erp_ramp, "fractional_area_latency", (; fractional_area_fraction = 1.0), 1.0, 0.01),
            (erp_tri, "fractional_peak_latency", (; fractional_peak_fraction = 0.0, fractional_peak_direction = :onset), 0.0, 0.01),
            (erp_tri, "fractional_peak_latency", (; fractional_peak_fraction = 0.0, fractional_peak_direction = :offset), 1.0, 0.01),
            (erp_tri, "fractional_peak_latency", (; fractional_peak_fraction = 1.0, fractional_peak_direction = :onset), 0.5, 0.01),
            (erp_tri, "fractional_peak_latency", (; fractional_peak_fraction = 0.5, fractional_peak_direction = :onset), 0.25, 0.02),
            (erp_tri, "fractional_peak_latency", (; fractional_peak_fraction = 0.5, fractional_peak_direction = :offset), 0.75, 0.02),
            (erp_bl5, "mean_amplitude", (; analysis_interval = (0.1, 0.2), baseline_interval = (-0.2, 0.0)), 0.0, 0.01),
        )
            iv = haskey(kw, :analysis_interval) ? kw.analysis_interval : (0.0, 1.0)
            res = EegFun.erp_measurements(erp, atype; analysis_interval = iv, kw...)
            @test isapprox(res.data[1, :Ch1], expected; atol = tol)
        end
        @test EegFun.erp_measurements(erp_robust, "max_peak_amplitude"; analysis_interval = (0.0, 1.0), local_interval = 5).data[1, :Ch1] >
              0.9

        plot_res = EegFun.plot_erp_measurements(
            [erp_bl5],
            "mean_amplitude";
            analysis_interval = [0.1, 0.2],
            baseline_interval = [-0.2, 0.0],
            channel_selection = :Ch1,
            display_plot = false,
        )
        @test (plot_res isa Figure) || (plot_res isa NamedTuple && plot_res.fig isa Figure)
    end

    @testset "In-memory erp_measurements, area_boundary, peak_polarity, and underscore channels" begin
        t_pts = collect(range(-0.2, 1.0, step = 0.01))
        raw_sig = [t < 0.0 ? 3.0 : 3.0 + (-4.0 + 8.0 * t) for t in t_pts]
        df_inmem = DataFrame(time = t_pts, sample = 1:length(t_pts), Ch1_bip = copy(raw_sig), vEOG = fill(999.0, length(t_pts)))
        lay_inmem = EegFun.Layout(DataFrame(label = [:Ch1_bip], inc = [0.0], azi = [0.0]), nothing, nothing, nothing)
        erp_inmem = EegFun.ErpData("sub1", 1, "condA", df_inmem, lay_inmem, 100, EegFun.AnalysisInfo(), 10)
        orig_col = copy(erp_inmem.data.Ch1_bip)

        res_rect = EegFun.erp_measurements(erp_inmem, "rectified_area"; analysis_interval = (0.0, 1.0), baseline_interval = (-0.2, -0.01))
        @test res_rect isa EegFun.ErpMeasurementsResult && erp_inmem.data.Ch1_bip == orig_col
        @test :Ch1_bip in propertynames(res_rect.data) && :vEOG ∉ propertynames(res_rect.data)
        @test isapprox(res_rect.data[1, :Ch1_bip], 2.0; atol = 0.05)

        for (atype, exp_val) in (("positive_area", 1.0), ("negative_area", 1.0), ("integral", 0.0))
            r = EegFun.erp_measurements(erp_inmem, atype; analysis_interval = (0.0, 1.0), baseline_interval = (-0.2, -0.01))
            @test isapprox(r.data[1, :Ch1_bip], exp_val; atol = 0.05)
        end

        for (ab, exp_lat) in ((:positive, 0.5 + 0.5 / sqrt(2)), (:rectified, 0.5))
            r = EegFun.erp_measurements(
                erp_inmem,
                "fractional_area_latency";
                analysis_interval = (0.0, 1.0),
                baseline_interval = (-0.2, -0.01),
                area_boundary = ab,
                fractional_area_fraction = 0.5,
            )
            @test isapprox(r.data[1, :Ch1_bip], exp_lat; atol = 0.015)
        end

        asym_sig = [-10.0 * exp(-((t - 0.2) / 0.05)^2) + 5.0 * exp(-((t - 0.8) / 0.05)^2) for t in t_pts]
        erp_asym = EegFun.ErpData(
            "asym",
            1,
            "c1",
            DataFrame(time = t_pts, sample = 1:length(t_pts), Ch1_bip = asym_sig),
            lay_inmem,
            100,
            EegFun.AnalysisInfo(),
            10,
        )
        lat_abs = EegFun.erp_measurements(
            erp_asym,
            "fractional_peak_latency";
            analysis_interval = (0.0, 1.0),
            fractional_peak_fraction = 0.5,
            peak_polarity = :abs,
        )
        lat_pos = EegFun.erp_measurements(
            erp_asym,
            "fractional_peak_latency";
            analysis_interval = (0.0, 1.0),
            fractional_peak_fraction = 0.5,
            peak_polarity = :positive,
        )
        lat_neg = EegFun.erp_measurements(
            erp_asym,
            "fractional_peak_latency";
            analysis_interval = (0.0, 1.0),
            fractional_peak_fraction = 0.5,
            peak_polarity = :negative,
        )
        @test isapprox(lat_abs.data[1, :Ch1_bip], lat_neg.data[1, :Ch1_bip]; atol = 1e-12) && lat_abs.data[1, :Ch1_bip] < 0.2
        @test 0.7 < lat_pos.data[1, :Ch1_bip] < 0.8

        ep_dfs = [
            let d = copy(erp_asym.data)
                d.epoch = fill(k, nrow(d))
                d
            end for k in (1, 2)
        ]
        ep_inmem = EegFun.EpochData("ep_sub", 1, "c1", ep_dfs, lay_inmem, 100, EegFun.AnalysisInfo())
        res_ep = EegFun.erp_measurements(ep_inmem, "mean_amplitude"; analysis_interval = (0.15, 0.25))
        @test nrow(res_ep.data) == 2 && :epoch in propertynames(res_ep.data)
        @test nrow(EegFun.erp_measurements([erp_inmem, erp_asym], "mean_amplitude"; analysis_interval = (0.1, 0.3)).data) == 2

        mp_erps = [create_known_signal_erp(p, c, t -> 10.0^c * p; t_start = 0.0, t_end = 0.5) for c in (1, 2) for p in (1, 2)]
        res_c2 =
            EegFun.erp_measurements(mp_erps, "mean_amplitude"; analysis_interval = (0.1, 0.4), condition_selection = EegFun.conditions(2))
        @test nrow(res_c2.data) == 2 && all(==(2), res_c2.data.condition) && sort(res_c2.data.participant) == [1, 2]

        zc_lay = EegFun.Layout(DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]), nothing, nothing, nothing)
        zc_erp = EegFun.ErpData(
            "zc",
            1,
            "c1",
            DataFrame(time = [0.0, 1.0], sample = [1, 2], Ch1 = [-2.0, 10.0]),
            zc_lay,
            1,
            EegFun.AnalysisInfo(),
            1,
        )
        @test EegFun.erp_measurements(zc_erp, "positive_area").data[1, :Ch1] ≈ 25.0 / 6.0
        @test EegFun.erp_measurements(zc_erp, "negative_area").data[1, :Ch1] ≈ 1.0 / 6.0
        @test EegFun.erp_measurements(zc_erp, "rectified_area").data[1, :Ch1] ≈ 26.0 / 6.0
        @test EegFun.erp_measurements(zc_erp, "integral").data[1, :Ch1] ≈ 4.0

        neg_only_erp = create_known_signal_erp(1, 1, t -> -5.0 - sin(π * t); fs = 100, t_start = 0.0, t_end = 1.0)
        @test isnan(
            EegFun.erp_measurements(neg_only_erp, "fractional_area_latency"; area_boundary = :positive, fractional_area_fraction = 0.5).data[
                1,
                :Ch1,
            ],
        )
        @test isnan(
            EegFun.erp_measurements(neg_only_erp, "fractional_peak_latency"; peak_polarity = :positive, fractional_peak_fraction = 0.5).data[
                1,
                :Ch1,
            ],
        )

        mut_c1 = EegFun.create_test_erp_data(participant = 1, condition = 1, fs = 100)
        mut_c2 = EegFun.create_test_erp_data(participant = 1, condition = 2, fs = 100)
        EegFun.erp_measurements!([mut_c1, mut_c2], "mean_amplitude")
        diff_after_mut = EegFun.condition_difference([mut_c1, mut_c2], (1, 2))
        @test EegFun.erp_measurements(diff_after_mut, "mean_amplitude").data[1, :condition_name] == "difference_1_2"
        @test EegFun.erp_measurements(
            mut_c1,
            "mean_amplitude";
            analysis_interval = df -> 0.1 .<= df.time .<= 0.3,
            baseline_interval = EegFun.times(),
        ) isa EegFun.ErpMeasurementsResult
        @test EegFun.erp_measurements(mut_c1, "mean_amplitude"; analysis_interval = EegFun.times(0.1234)) isa EegFun.ErpMeasurementsResult
    end
end
