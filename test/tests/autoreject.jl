using Test
using DataFrames
using OrderedCollections
using Statistics
using LinearAlgebra
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

        # Default PreprocessConfig has autoreject.apply == false
        prep_cfg_default = PreprocessConfig()
        @test prep_cfg_default.autoreject.apply == false
        @test PreprocessConfig(Dict{String,Any}()).autoreject.apply == false

        # Default config generation has preprocess.autoreject.apply == false
        default_dict = EegFun._generate_default_config()
        @test haskey(default_dict, "preprocess")
        @test haskey(default_dict["preprocess"], "autoreject")
        @test default_dict["preprocess"]["autoreject"]["apply"] == false
    end

    @testset "autoreject on clean synthetic epochs" begin
        # 12 clean epochs with identical sinusoidal signals across 6 channels
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

    @testset "autoreject non-mutation of input EpochData" begin
        dat = create_test_epoch_data(n_epochs = 8, n_channels = 4, fs = 250, n = 80, seed = 42)
        orig_dfs = [copy(df) for df in dat.data]
        orig_n = n_epochs(dat)

        cleaned, res = autoreject(dat; n_folds = 2, seed = 42, verbose = false)

        @test n_epochs(dat) == orig_n
        for (ep_i, df_orig) in enumerate(orig_dfs)
            @test Matrix(df_orig[:, [:Ch1, :Ch2, :Ch3, :Ch4]]) == Matrix(dat.data[ep_i][:, [:Ch1, :Ch2, :Ch3, :Ch4]])
        end
    end

    @testset "autoreject determinism with fixed seed" begin
        dat1 = create_test_epoch_data(n_epochs = 10, n_channels = 4, fs = 250, n = 80, seed = 10)
        dat2 = create_test_epoch_data(n_epochs = 10, n_channels = 4, fs = 250, n = 80, seed = 10)

        # Inject same artifact
        dat1.data[2][1:15, :Ch1] .+= 250.0
        dat2.data[2][1:15, :Ch1] .+= 250.0

        _, res1 = autoreject(dat1; n_folds = 3, seed = 777, verbose = false)
        _, res2 = autoreject(dat2; n_folds = 3, seed = 777, verbose = false)

        @test res1.thresholds == res2.thresholds
        @test res1.epochs_clean == res2.epochs_clean
        @test res1.epochs_repaired == res2.epochs_repaired
        @test res1.epochs_rejected == res2.epochs_rejected
    end

    @testset "autoreject synthetic recovery and RMSE reduction" begin
        layout = create_test_layout(layout_type = :topo)
        n_eps = 15
        channels = [:Fp1, :Fp2, :Cz, :Pz, :O1, :O2]
        time = collect(0:99) ./ 250.0
        dfs = DataFrame[]

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
            push!(dfs, df)
        end

        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        # Ground truth clean signal for Cz in epoch 3
        clean_cz = copy(ep_data.data[3][!, :Cz])

        # Inject transient artifact into Cz on epoch 3
        ep_data.data[3][20:40, :Cz] .+= 300.0
        corrupted_rmse = sqrt(mean((ep_data.data[3][!, :Cz] .- clean_cz) .^ 2))

        cleaned_data, res =
            autoreject(ep_data; mode = :local, n_folds = 3, consensus_perc = 0.5, max_repair_perc = 0.2, seed = 42, verbose = false)

        @test 3 in res.epochs_repaired
        @test haskey(res.repaired_channels, 3)
        @test :Cz in res.repaired_channels[3]

        repaired_cz = cleaned_data.data[3][!, :Cz]
        repaired_rmse = sqrt(mean((repaired_cz .- clean_cz) .^ 2))

        # Repaired signal should significantly improve RMSE against original clean signal
        @test repaired_rmse < corrupted_rmse * 0.3
    end

    @testset "autoreject middle band (consensus vs repair threshold)" begin
        # 6 channels: consensus_perc = 0.5 (k_max = 3), max_repair_perc = 0.17 (r_max = 1)
        # Epoch with 2 bad channels (r_max < 2 <= k_max):
        # Epoch should be RETAINED (not rejected).
        # The worst 1 channel should be repaired; the 2nd should be tracked in skipped_channels.
        layout = create_test_layout(layout_type = :topo)
        channels = layout.data.label
        n_eps = 15
        time = collect(0:99) ./ 250.0
        dfs = DataFrame[]

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
            push!(dfs, df)
        end

        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        # Inject 2 bad channels into epoch 4 with different severity
        ep_data.data[4][1:20, channels[1]] .+= 200.0  # mild
        ep_data.data[4][1:20, channels[2]] .+= 450.0  # severe

        cleaned, res = autoreject(
            ep_data;
            channel_selection = EegFun.channels(channels),
            consensus_perc = 0.5,
            max_repair_perc = 0.17,
            n_folds = 3,
            seed = 42,
            verbose = false,
        )

        @test 4 in res.epochs_repaired
        @test 4 ∉ res.epochs_rejected
        @test haskey(res.repaired_channels, 4)
        @test length(res.repaired_channels[4]) == 1
        # Worst channel (channels[2]) should be repaired
        @test channels[2] in res.repaired_channels[4]
        # Mild channel should be in skipped_channels
        @test haskey(res.skipped_channels, 4)
        @test channels[1] in res.skipped_channels[4]
    end

    @testset "autoreject failed repair moves to rejected" begin
        # Create epoch data with an empty layout (no coordinates, no neighbours)
        channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        layout = EegFun.Layout(DataFrame(label = channels), nothing, nothing, nothing)
        ep_data = create_test_epoch_data(n_epochs = 10, n_channels = 4, fs = 250, n = 50, seed = 42)
        ep_data.layout = layout

        # Inject artifact in epoch 2
        ep_data.data[2][1:20, :Ch1] .+= 300.0

        cleaned, res = autoreject(ep_data; n_folds = 3, consensus_perc = 0.5, max_repair_perc = 0.3, seed = 42, verbose = false)

        # Because layout has no coordinates and no neighbours, repair fails
        # Epoch 2 must be moved to rejected, NOT kept as clean/repaired!
        @test 2 in res.epochs_rejected
        @test 2 ∉ res.epochs_repaired
    end

    @testset "autoreject fallback for < 4 epochs" begin
        small_ep = create_test_epoch_data(n_epochs = 2, n_channels = 3, fs = 250, n = 50, seed = 42)
        cleaned_small, res_small = autoreject(small_ep; verbose = false)

        @test res_small isa AutoRejectResult
        @test isempty(res_small.epochs_rejected)
        @test isempty(res_small.epochs_repaired)
        @test length(res_small.epochs_clean) == 2
        @test n_epochs(cleaned_small) == 2
    end

    @testset "autoreject edge cases and validation" begin
        ep_data = create_test_epoch_data(n_epochs = 6, n_channels = 3, fs = 1000, n = 1000, seed = 42)

        # Invalid mode throws error
        @test_throws Exception autoreject(ep_data; mode = :invalid_mode)

        # Invalid repair method throws error
        @test_throws Exception autoreject(ep_data; repair_method = :invalid_method)

        # Invalid n_folds < 2 throws error
        @test_throws Exception autoreject(ep_data; n_folds = 1)

        # Invalid consensus_perc and max_repair_perc throw error
        @test_throws Exception autoreject(ep_data; consensus_perc = 1.5)
        @test_throws Exception autoreject(ep_data; consensus_perc = -0.2)
        @test_throws Exception autoreject(ep_data; max_repair_perc = -0.1)
        @test_throws Exception autoreject(ep_data; max_repair_perc = 1.2)

        # Interval selection
        cleaned_win, res_win = autoreject(ep_data; interval_selection = times(0.1, 0.4), n_folds = 2, verbose = false)
        @test res_win isa AutoRejectResult
    end

    @testset "autoreject global mode" begin
        layout = create_test_layout(layout_type = :topo)
        n_eps = 12
        channels = [:Fp1, :Fp2, :Cz, :Pz]
        time = collect(0:99) ./ 250.0
        dfs = DataFrame[]

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
                df[1:20, :Fp1] .+= 300.0
            end
            push!(dfs, df)
        end

        ep_data = EpochData("test", 1, "cond1", dfs, layout, 250, EegFun.AnalysisInfo())

        cleaned_epochs, res = autoreject(ep_data; mode = :global, n_folds = 3, seed = 42, verbose = false)

        @test res.mode == :global
        @test length(res.epochs_repaired) == 0  # Global mode does not repair
        @test 3 in res.epochs_rejected
        @test length(unique(values(res.thresholds))) == 1  # single global threshold
        @test n_epochs(cleaned_epochs) == 11
    end

    @testset "autoreject on Vector{EpochData} (multiple conditions)" begin
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
        @test any(r -> r.epoch == 3, res.rejection_info[1].rejected)
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
