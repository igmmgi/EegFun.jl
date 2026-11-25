# Tests for the standalone HTML preprocessing pipeline QC report generator (`pipeline_report`).

using Test
using CairoMakie
using CSV
using DataFrames
using JLD2
using OrderedCollections
using EegFun

@testset "Pipeline HTML QC Report (pipeline_report)" begin
    @testset "Validation and error guards" begin
        @test_throws EegFun.EegFunError EegFun.pipeline_report("/nonexistent/eegfun_dir_xyz_999")
        mktempdir() do empty_dir
            @test_throws EegFun.EegFunError EegFun.pipeline_report(empty_dir)
            @test_throws EegFun.EegFunError EegFun.pipeline_report(empty_dir; low_epoch_threshold = -5.0)
            @test_throws EegFun.EegFunError EegFun.pipeline_report(empty_dir; low_epoch_threshold = 105.0)
            @test_throws EegFun.EegFunError EegFun.pipeline_report(empty_dir; max_ica_pct = -1.0)
            @test_throws EegFun.EegFunError EegFun.pipeline_report(empty_dir; max_ica_pct = 120.0)
        end
    end

    @testset "Full multi-participant report with plots, flags, ICA, and config" begin
        mktempdir() do out_dir
            # Write realistic config.toml
            EegFun.generate_config_template(; filename = joinpath(out_dir, "config.toml"))

            # Participant 1: Clean dataset (100% / 95% epochs kept, 1 repaired channel, 2 removed ICs)
            topo_layout = EegFun.Layout(
                DataFrame(
                    label = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5, :Ch6],
                    inc = [65.0, 65.0, 0.0, 45.0, 80.0, 80.0],
                    azi = [-35.0, 35.0, 0.0, 180.0, -145.0, 145.0],
                ),
                nothing,
                nothing,
                nothing,
            )
            EegFun.polar_to_cartesian_xy!(topo_layout)
            EegFun.polar_to_cartesian_xyz!(topo_layout)

            cdat1 = EegFun.create_test_continuous_data(n = 600, fs = 250, n_channels = 6, seed = 101)
            cdat1.layout = copy(topo_layout)
            ica1 = EegFun.run_ica(cdat1; params = EegFun.IcaPrms(max_iter = 5))
            erps1 = [
                EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 6),
                EegFun.create_test_erp_data(participant = 1, condition = 2, n_channels = 6),
            ]
            erps1[1].layout = copy(topo_layout)
            erps1[2].layout = copy(topo_layout)

            rep1 = EegFun.ContinuousRepairInfo("continuous_repair", :spherical_spline, [:Ch2], Symbol[])
            step1_s1 = [
                EegFun.EpochRejectionInfo(
                    "rejection_step1",
                    EegFun.EpochInfo(1, "Standard", 40),
                    2,
                    100.0,
                    [EegFun.Rejection(:Ch1, 3), EegFun.Rejection(:Ch2, 7)],
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch1, 3), EegFun.Rejection(:Ch2, 7)],
                    OrderedDict(3 => [:Ch1]),
                    OrderedDict{Int,Vector{Symbol}}(),
                ),
                EegFun.EpochRejectionInfo(
                    "rejection_step1",
                    EegFun.EpochInfo(2, "Target", 40),
                    2,
                    100.0,
                    [EegFun.Rejection(:Ch3, 5), EegFun.Rejection(:Ch4, 12)],
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch3, 5), EegFun.Rejection(:Ch4, 12)],
                    OrderedDict{Int,Vector{Symbol}}(),
                    OrderedDict{Int,Vector{Symbol}}(),
                ),
            ]
            step2_s1 = [
                EegFun.EpochRejectionInfo(
                    "rejection_step2",
                    EegFun.EpochInfo(1, "Standard", 40),
                    1,
                    100.0,
                    [EegFun.Rejection(:Ch2, 7)],
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch2, 7)],
                    nothing,
                    nothing,
                ),
                EegFun.EpochRejectionInfo(
                    "rejection_step2",
                    EegFun.EpochInfo(2, "Target", 40),
                    2,
                    100.0,
                    [EegFun.Rejection(:Ch3, 5), EegFun.Rejection(:Ch4, 12)],
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch3, 5), EegFun.Rejection(:Ch4, 12)],
                    nothing,
                    nothing,
                ),
            ]
            ics1 = EegFun.ArtifactComponents(Dict(:vEOG => [1], :hEOG => Int[]), Int[], Int[], Int[], [2])
            art1 = EegFun.ArtifactInfo([rep1], vcat(step1_s1, step2_s1), ics1)

            erps1_raw = [copy(erps1[1]), copy(erps1[2])]
            erps1_raw[1].data[!, :Ch1] .+= 25.0 # simulate pre-ICA frontal artifact
            EegFun.save_data(joinpath(out_dir, "sub01_artifact_info.jld2"), art1)
            EegFun.save_data(joinpath(out_dir, "sub01_ica.jld2"), ica1)
            EegFun.save_data(joinpath(out_dir, "sub01_erps.jld2"), erps1)
            EegFun.save_data(joinpath(out_dir, "sub01_erps_raw.jld2"), erps1_raw)

            # Participant 2: Flagged dataset (skipped channel + <70% retention in Target + HTML-special condition name)
            erps2 = [
                EegFun.create_test_erp_data(participant = 2, condition = 1, n_channels = 6),
                EegFun.create_test_erp_data(participant = 2, condition = 2, n_channels = 6),
            ]
            erps2[1].layout = copy(cdat1.layout)
            erps2[2].layout = copy(cdat1.layout)
            rep2 = EegFun.ContinuousRepairInfo("continuous_repair", :neighbor_interpolation, [:Ch2], [:Ch5])
            rej2 = [
                EegFun.EpochRejectionInfo(
                    "autoreject",
                    EegFun.EpochInfo(1, "Standard", 20),
                    2,
                    100.0,
                    nothing,
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch1, i) for i = 1:2],
                    OrderedDict(4 => [:Ch2]),
                    OrderedDict{Int,Vector{Symbol}}(),
                ),
                EegFun.EpochRejectionInfo(
                    "autoreject",
                    EegFun.EpochInfo(2, "Target <A&B>", 20),
                    10,
                    100.0,
                    nothing,
                    0.0,
                    nothing,
                    [EegFun.Rejection(:Ch1, i) for i = 1:10],
                    OrderedDict{Int,Vector{Symbol}}(),
                    OrderedDict{Int,Vector{Symbol}}(),
                ),
            ]
            art2 = EegFun.ArtifactInfo([rep2], rej2, nothing)
            EegFun.save_data(joinpath(out_dir, "sub02_artifact_info.jld2"), art2)
            EegFun.save_data(joinpath(out_dir, "sub02_erps.jld2"), erps2)

            report_path = EegFun.pipeline_report(out_dir; title = "Study <QC> & Summary")
            @test report_path == abspath(joinpath(out_dir, "pipeline_report.html"))
            @test isfile(report_path)

            html = read(report_path, String)
            @test occursin("Study &lt;QC&gt; &amp; Summary", html)
            @test occursin("Target &lt;A&amp;B&gt;", html)
            @test occursin("sub01", html)
            @test occursin("sub02", html)
            @test occursin("1 / 2", html) # 1 of 2 flagged
            @test occursin("Ch5", html)   # skipped channel badge
            @test occursin("Manual:1", html)
            @test occursin("Cohort Grand-Average &amp; Spectral Preview", html)
            @test occursin("Draft Methods &amp; Attrition Boilerplate", html)
            @test occursin("avg neighbors =", html)
            @test occursin("download=\"pipeline_qc_summary.csv\"", html)
            @test occursin("Flagged Only (1)", html)
            @test occursin("Cohort sensor repair and interpolation heatmap", html)
            @test occursin("Cohort pre- vs post-cleaning ERPs and GFP", html)
            @test occursin("Cohort cleaned power spectral density", html)
            @test occursin("Pre- vs post-cleaning overlay for sub01", html)
            @test !occursin("sub01_erps_raw</a>", html) # raw files must not spawn duplicate participant stems
            @test occursin("Expand All", html)
            @test occursin("Collapse All", html)
            @test occursin("77 / 80", html) # sub01 total kept / original epochs
            @test count("data:image/png;base64,", html) >= 10

            # Cohort-level ICA % and truncated trial count review flags (≥3 participants)
            erps3 = [
                EegFun.create_test_erp_data(participant = 3, condition = 1, n_channels = 6),
                EegFun.create_test_erp_data(participant = 3, condition = 2, n_channels = 6),
            ]
            erps3[1].layout = copy(topo_layout)
            erps3[2].layout = copy(topo_layout)
            ics3 = EegFun.ArtifactComponents(Dict(:vEOG => [1, 2], :hEOG => [2, 3]), Int[], Int[], Int[], Int[])
            res_dup = EegFun.plot_artifact_components(ica1, ics3; display_plot = false)
            @test length(res_dup.axes) == 3
            @test res_dup.axes[2].title[] == "vEOG + hEOG 2"

            rej3 = [
                EegFun.EpochRejectionInfo(
                    "autoreject",
                    EegFun.EpochInfo(1, "Standard", 40),
                    0,
                    100.0,
                    nothing,
                    0.0,
                    nothing,
                    EegFun.Rejection[],
                    OrderedDict(2 => [:Ch1, :Ch2, :Ch3], 5 => [:Ch1, :Ch4, :Ch5, :Ch6]),
                    nothing,
                ),
                EegFun.EpochRejectionInfo(
                    "autoreject",
                    EegFun.EpochInfo(2, "Target", 40),
                    0,
                    100.0,
                    nothing,
                    0.0,
                    nothing,
                    EegFun.Rejection[],
                    nothing,
                    nothing,
                ),
            ]
            art3 = EegFun.ArtifactInfo(EegFun.ContinuousRepairInfo[], rej3, ics3)
            EegFun.save_data(joinpath(out_dir, "sub03_artifact_info.jld2"), art3)
            EegFun.save_data(joinpath(out_dir, "sub03_ica.jld2"), ica1)
            EegFun.save_data(joinpath(out_dir, "sub03_erps.jld2"), erps3)

            cohort3_path = EegFun.pipeline_report(out_dir; output_file = "cohort3_report.html", include_plots = false)
            html_c3 = read(cohort3_path, String)
            @test occursin("high ICA removal: 3/6 components", html_c3)
            @test occursin("1 shared", html_c3)
            @test occursin("low trial count: 40 original epochs", html_c3) # sub02 has 40 vs cohort median 80
            @test occursin("id=\"cohort-summary\"", html_c3)
            @test occursin("↑ Summary", html_c3)
            @test occursin("openHashTarget", html_c3)
            @test !occursin("Continuous channel status (—)", html)

            # Verify tidy on-disk `pipeline_qc_summary.csv` schema, values, and top-5 channel truncation
            csv_path = joinpath(out_dir, "pipeline_qc_summary.csv")
            @test isfile(csv_path)
            qc_df = CSV.read(csv_path, DataFrame)
            @test nrow(qc_df) == 3
            @test names(qc_df)[1:7] ==
                  ["participant", "status", "epochs_original", "epochs_repaired", "epochs_rejected", "epochs_kept", "epochs_kept_pct"]
            @test all(
                col -> col in names(qc_df),
                [
                    "Standard_original",
                    "Standard_kept",
                    "Standard_kept_pct",
                    "Target_original",
                    "Target_kept",
                    "Target_kept_pct",
                    "n_channels",
                    "n_repaired_channels",
                    "n_skipped_channels",
                    "n_epoch_interp_channels",
                    "n_epoch_channel_interps",
                    "ic_available",
                    "ic_removed",
                    "ic_removed_pct",
                    "repaired_channels",
                    "skipped_channels",
                    "top_epoch_interp_channels",
                    "flag_reasons",
                ],
            )
            @test qc_df.epochs_original == [80, 40, 80]
            @test qc_df.epochs_repaired == [1, 1, 2]
            @test qc_df.epochs_rejected == [3, 12, 0]
            @test qc_df.epochs_kept == [77, 28, 80]
            @test qc_df.epochs_kept_pct == [96.2, 70.0, 100.0]
            @test qc_df.Standard_kept == [39, 18, 40]
            @test qc_df.Target_kept == [38, 10, 40]
            @test qc_df.n_channels == [6, 6, 6]
            @test qc_df.n_repaired_channels == [1, 1, 0]
            @test qc_df.n_skipped_channels == [0, 1, 0]
            @test qc_df.n_epoch_interp_channels == [1, 1, 6]
            @test qc_df.n_epoch_channel_interps == [1, 1, 7]
            @test qc_df.ic_available[1] == 6
            @test ismissing(qc_df.ic_available[2])
            @test qc_df.ic_removed == [2, 0, 3]
            @test qc_df.ic_removed_pct[1] == 33.3
            @test qc_df.ic_removed_pct[3] == 50.0
            @test occursin("(+1 more)", qc_df.top_epoch_interp_channels[3])
            @test names(qc_df)[end] == "flag_reasons"

            # Tabular-only mode (`include_plots = false`, `save_csv = false`) and custom output path
            custom_html_path = joinpath(out_dir, "subreports", "tabular_only.html")
            ret_custom = EegFun.pipeline_report(out_dir; output_file = custom_html_path, save_csv = false, include_plots = false)
            @test ret_custom == abspath(custom_html_path)
            @test isfile(ret_custom)
            @test !isfile(joinpath(out_dir, "subreports", "pipeline_qc_summary.csv"))
            html_no_plots = read(ret_custom, String)
            @test !occursin("data:image/png;base64,", html_no_plots)
            @test occursin("sub01", html_no_plots)
            @test occursin("sub02", html_no_plots)
        end
    end

    @testset "Partial output directory (epochs-only and coordinate-less layout fallback)" begin
        mktempdir() do out_dir
            ep = EegFun.create_test_epoch_data_vector(n = 50, n_epochs = 6, n_channels = 4)
            EegFun.save_data(joinpath(out_dir, "p1_epochs.jld2"), ep)
            p = EegFun.pipeline_report(out_dir; include_plots = true)
            @test isfile(p)
            html_partial = read(p, String)
            @test occursin("p1", html_partial)
            @test count("data:image/png;base64,", html_partial) >= 1
        end
    end
end
