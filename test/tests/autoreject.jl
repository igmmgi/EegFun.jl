using Test
using DataFrames
using OrderedCollections
using Statistics
using EegFun
using EegFun:
    AutoRejectConfig,
    AutoRejectResult,
    autoreject,
    autoreject!,
    PreprocessConfig,
    EpochData,
    EpochRejectionInfo,
    n_epochs,
    times,
    create_test_epoch_data,
    create_test_epoch_data_vector,
    create_test_layout

@testset "autoreject" begin

    @testset "AutoRejectConfig and PipelineConfig integration" begin
        cfg = AutoRejectConfig()
        @test cfg.apply == false
        @test cfg.mode == :local
        @test cfg.n_folds == 10
        @test cfg.n_candidates == 30
        @test cfg.consensus_perc == 0.2
        @test cfg.max_repair_perc == 0.1
        @test cfg.repair_method == :spherical_spline
        @test cfg.seed == 42

        # Dictionary constructor
        dict = Dict(
            "apply" => true,
            "mode" => "global",
            "n_folds" => 5,
            "n_candidates" => 20,
            "consensus_perc" => 0.3,
            "max_repair_perc" => 0.15,
            "repair_method" => "neighbor_interpolation",
            "seed" => 123,
        )
        parsed_cfg = AutoRejectConfig(dict)
        @test parsed_cfg.apply == true
        @test parsed_cfg.mode == :global
        @test parsed_cfg.n_folds == 5
        @test parsed_cfg.n_candidates == 20
        @test parsed_cfg.consensus_perc == 0.3
        @test parsed_cfg.max_repair_perc == 0.15
        @test parsed_cfg.repair_method == :neighbor_interpolation
        @test parsed_cfg.seed == 123

        # Nested in PreprocessConfig
        prep_cfg = PreprocessConfig(Dict("autoreject" => dict))
        @test prep_cfg.autoreject.apply == true
        @test prep_cfg.autoreject.mode == :global
    end

    @testset "autoreject on clean synthetic epochs" begin
        # 12 clean epochs with identical sinusoidal signals across 4 channels
        time = collect(0:99) ./ 250.0
        n_eps = 12
        layout = create_test_layout(layout_type = :topo)
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        dfs = [
            DataFrame(
                time = time,
                sample = 1:100,
                epoch = fill(ep, 100),
                participant = fill(1, 100),
                condition = fill(1, 100),
                condition_name = fill("cond1", 100),
                Fp1 = sin.(2 * π * 5 .* time),
                Fp2 = sin.(2 * π * 5 .* time),
                Cz = cos.(2 * π * 5 .* time),
                Pz = cos.(2 * π * 5 .* time),
                O1 = sin.(2 * π * 10 .* time),
                O2 = sin.(2 * π * 10 .* time),
            ) for ep = 1:n_eps
        ]
        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        # Non-mutating run
        original_n = n_epochs(ep_data)
        cleaned_epochs, result = autoreject(ep_data; n_folds = 3, n_candidates = 15, seed = 42, verbose = false)

        @test result isa AutoRejectResult
        @test n_epochs(ep_data) == original_n  # original untouched
        @test n_epochs(cleaned_epochs) == original_n  # all clean epochs kept
        @test length(result.epochs_rejected) == 0
        @test length(result.epochs_clean) == original_n
        @test length(result.thresholds) == 6
        @test all(v -> v > 0, values(result.thresholds))

        # Check summary table
        @test result.summary_table isa DataFrame
        @test "Retention Rate (%)" in result.summary_table.Metric
        @test result.channel_table isa DataFrame
        @test nrow(result.channel_table) == 6

        # Check Base.show
        io = IOBuffer()
        show(io, result)
        str = String(take!(io))
        @test occursin("AutoRejectResult", str)
        @test occursin("kept", str)
    end

    @testset "autoreject local repair and rejection" begin
        # Create 15 epochs, 6 channels with 3D coordinates (topo layout)
        layout = create_test_layout(layout_type = :topo)
        n_eps = 15
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        dfs = DataFrame[]
        time = collect(0:99) ./ 250.0

        for ep = 1:n_eps
            df = DataFrame(
                time = time,
                sample = 1:100,
                epoch = fill(ep, 100),
                participant = fill(1, 100),
                condition = fill(1, 100),
                condition_name = fill("cond1", 100),
            )
            for ch in channels
                df[!, ch] = sin.(2 * π * 5 .* time)
            end

            # Inject a single bad channel artifact in epoch 2 (1 bad ch out of 6 -> repaired)
            if ep == 2
                df[1:20, :Fp1] .+= 250.0
            end

            # Inject all channels bad in epoch 5 (6 bad chs out of 6 -> rejected)
            if ep == 5
                for ch in channels
                    df[1:20, ch] .+= 300.0
                end
            end

            push!(dfs, df)
        end

        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        # Run local autoreject with max_repair_perc = 0.2 (allows 1 bad ch out of 6 to repair)
        # and consensus_perc = 0.4 (rejects if > 2 bad chs)
        cleaned_epochs, res = autoreject(
            ep_data;
            mode = :local,
            n_folds = 3,
            n_candidates = 20,
            consensus_perc = 0.4,
            max_repair_perc = 0.2,
            seed = 42,
            verbose = false,
        )

        @test 2 in res.epochs_repaired
        @test 5 in res.epochs_rejected
        @test haskey(res.repaired_channels, 2)
        @test :Fp1 in res.repaired_channels[2]

        # Cleaned epoch 2 should have had Fp1 repaired
        f_p1_repaired = cleaned_epochs.data[2][!, :Fp1]
        @test maximum(abs, f_p1_repaired) < 150.0

        # Epoch 5 should have been pruned
        @test n_epochs(cleaned_epochs) == 14

        # Test reject = false: returns all 15 epochs, but res still records rejection
        cleaned_no_prune, res_no_prune = autoreject(
            ep_data;
            mode = :local,
            n_folds = 3,
            n_candidates = 20,
            consensus_perc = 0.4,
            max_repair_perc = 0.2,
            seed = 42,
            reject = false,
            verbose = false,
        )
        @test n_epochs(cleaned_no_prune) == 15
        @test 5 in res_no_prune.epochs_rejected
    end

    @testset "autoreject global mode" begin
        layout = create_test_layout(layout_type = :topo)
        n_eps = 10
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        dfs = DataFrame[]
        time = collect(0:99) ./ 250.0

        for ep = 1:n_eps
            df = DataFrame(
                time = time,
                sample = 1:100,
                epoch = fill(ep, 100),
                participant = fill(1, 100),
                condition = fill(1, 100),
                condition_name = fill("cond1", 100),
            )
            for ch in channels
                df[!, ch] = sin.(2 * π * 5 .* time)
            end
            if ep == 3
                df[1:20, :Cz] .+= 350.0  # pulse artifact
            end
            push!(dfs, df)
        end

        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        cleaned_epochs, res = autoreject(ep_data; mode = :global, n_folds = 3, seed = 42, verbose = false)

        @test res.mode == :global
        # In global mode: single threshold across all channels, no repair
        all_th = collect(values(res.thresholds))
        @test all(th -> th == all_th[1], all_th)
        @test length(res.epochs_repaired) == 0
        @test 3 in res.epochs_rejected
        @test n_epochs(cleaned_epochs) == 9
    end

    @testset "autoreject on Vector{EpochData} (multiple conditions)" begin
        # 2 conditions, 8 epochs each
        ep_vec = create_test_epoch_data_vector(conditions = 1:2, n_epochs = 8, n_channels = 4, fs = 250, n = 100, seed = 42)

        # Inject pulse artifact in condition 1, epoch 3
        ep_vec[1].data[3][1:20, :Ch1] .+= 400.0

        cleaned_vec, res = autoreject(ep_vec; n_folds = 3, consensus_perc = 0.2, seed = 42, verbose = false)

        @test length(cleaned_vec) == 2
        @test res isa AutoRejectResult
        @test length(res.rejection_info) == 2
        @test res.rejection_info[1] isa EpochRejectionInfo
        @test res.rejection_info[2] isa EpochRejectionInfo

        # Condition 1 epoch 3 was rejected
        @test 3 in res.rejection_info[1].rejected[1].epoch
    end

    @testset "autoreject edge cases" begin
        # Very few epochs (e.g. 2 epochs)
        small_ep = create_test_epoch_data(n_epochs = 2, n_channels = 3, fs = 250, n = 50, seed = 42)
        cleaned_small, res_small = autoreject(small_ep; verbose = false)
        @test res_small isa AutoRejectResult
        @test length(res_small.thresholds) == 3

        # Invalid mode throws error
        @test_throws Exception autoreject(small_ep; mode = :invalid_mode)

        # Invalid consensus_perc and max_repair_perc throw error
        @test_throws Exception autoreject(small_ep; consensus_perc = 1.5)
        @test_throws Exception autoreject(small_ep; consensus_perc = -0.2)
        @test_throws Exception autoreject(small_ep; max_repair_perc = -0.1)
        @test_throws Exception autoreject(small_ep; max_repair_perc = 1.2)

        # Interval selection
        ep_data = create_test_epoch_data(n_epochs = 6, n_channels = 3, fs = 1000, n = 1000, seed = 42)
        cleaned_win, res_win = autoreject(ep_data; interval_selection = times(0.1, 0.4), n_folds = 2, verbose = false)
        @test res_win isa AutoRejectResult
    end

    @testset "pipeline preprocess with autoreject" begin
        data_dir = joinpath(dirname(dirname(@__DIR__)), "resources", "data")
        test_bdf = joinpath(data_dir, "bdf", "example1.bdf")

        if isfile(test_bdf)
            mktempdir() do temp_dir
                epoch_file = joinpath(temp_dir, "epochs.toml")
                write(
                    epoch_file,
                    """
[epochs]
[[epochs.conditions]]
name = "Test Condition"
trigger_sequences = [1]
""",
                )

                pipeline_file = joinpath(temp_dir, "pipeline.toml")
                bdf_dir = replace(joinpath(data_dir, "bdf"), '\\' => '/')
                write(
                    pipeline_file,
                    """
[files.input]
directory = "$bdf_dir"
raw_data_files = "example1\\\\.bdf"
layout_file = "biosemi72.csv"
epoch_condition_file = "epochs.toml"

[files.output]
directory = "output"
save_continuous_data_raw = false
save_continuous_data_corrected = false

[preprocess]
epoch_start = -0.1
epoch_end = 0.5

[preprocess.autoreject]
apply = true
mode = "local"
n_folds = 3
n_candidates = 15
consensus_perc = 0.3
max_repair_perc = 0.1
""",
                )

                redirect_stdio(stdout = devnull, stderr = devnull) do
                    EegFun.preprocess(pipeline_file; base_dir = temp_dir)
                end

                out_dir = joinpath(temp_dir, "output")
                @test isdir(out_dir)
                @test isfile(joinpath(out_dir, "example1_epochs.jld2"))
                @test isfile(joinpath(out_dir, "example1_erps.jld2"))
                @test isfile(joinpath(out_dir, "example1_artifact_info.jld2"))
            end
        end
    end

end
