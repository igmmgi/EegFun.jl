# Tests for between-condition difference waveforms and time-frequency contrasts.

using Test
using JLD2
using DataFrames
using CSV
using EegFun

@testset "Batch Difference Conditions" begin
    isdefined(@__MODULE__, :run_batch_jld2_contract_tests!) ||
        Base.include(@__MODULE__, joinpath(pkgdir(EegFun), "test", "tests", "test_helpers.jl"))

    run_batch_jld2_contract_tests!(
        (pat; kwargs...) -> EegFun.condition_difference(pat, [(1, 2), (3, 4)]; kwargs...);
        make_data = p -> [EegFun.create_test_erp_data(participant = p, condition = c) for c = 1:4],
        pattern = "erps_unrejected",
        n_participants = 3,
        supports_condition_selection = false,
        log_file = "condition_difference.log",
        log_snippets = ["Batch condition differencing started", "Batch operation complete"],
        expected_auto_dir = "differences_erps_unrejected_1-2_3-4",
    )

    mktempdir() do test_dir
        orig1 = [EegFun.create_test_erp_data(participant = 1, condition = c) for c = 1:4]
        jldsave(joinpath(test_dir, "1_erps_unrejected.jld2"); data = orig1)
        jldsave(
            joinpath(test_dir, "99_erps_unrejected.jld2");
            data = [EegFun.create_test_erp_data(participant = 99, condition = c) for c = 1:2],
        )

        # Vector pairs syntax, subtraction integrity, metadata & partial missing-condition handling
        out_dir = joinpath(test_dir, "out")
        EegFun.condition_difference("erps_unrejected", [[1, 2], [3, 4]]; input_dir = test_dir, output_dir = out_dir)
        diffs1 = load(joinpath(out_dir, "1_erps_unrejected.jld2"), "data")
        @test [d.condition_name for d in diffs1] == ["difference_1_2", "difference_3_4"]
        @test [d.condition for d in diffs1] == [1, 2]
        @test diffs1[1].data.Ch1 ≈ orig1[1].data.Ch1 .- orig1[2].data.Ch1
        @test diffs1[1].sample_rate == orig1[1].sample_rate
        @test diffs1[1].n_epochs == min(orig1[1].n_epochs, orig1[2].n_epochs)
        @test length(load(joinpath(out_dir, "99_erps_unrejected.jld2"), "data")) == 1

        @test_throws Exception EegFun.condition_difference("erps_unrejected", []; input_dir = test_dir)
        @test_throws Exception EegFun.condition_difference("erps_unrejected", [(1, "invalid")]; input_dir = test_dir)
    end
end

@testset "In-Memory Difference Conditions" begin

    @testset "Basic difference wave" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        diff = EegFun.condition_difference(erps, [(1, 2)])

        @test length(diff) == 1
        @test diff[1].condition == 1  # pair index
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "Multiple difference pairs" begin
        erps = [
            EegFun.create_test_erp_data(condition = 1),
            EegFun.create_test_erp_data(condition = 2),
            EegFun.create_test_erp_data(condition = 3),
            EegFun.create_test_erp_data(condition = 4),
        ]

        diff = EegFun.condition_difference(erps, [(1, 2), (3, 4)])

        @test length(diff) == 2
        @test diff[1].condition_name == "difference_1_2"
        @test diff[2].condition_name == "difference_3_4"
    end

    @testset "Vector syntax" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        diff = EegFun.condition_difference(erps, [[1, 2]])

        @test length(diff) == 1
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "Data integrity" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        diff = EegFun.condition_difference(erps, [(1, 2)])

        # Verify difference is actually erp1 - erp2 for EEG channels
        for ch in [:Fz, :Cz, :Pz]
            if hasproperty(erps[1].data, ch) && hasproperty(erps[2].data, ch) && hasproperty(diff[1].data, ch)
                expected = erps[1].data[!, ch] .- erps[2].data[!, ch]
                @test all(abs.(diff[1].data[!, ch] .- expected) .< 1e-10)
            end
        end
    end

    @testset "Metadata preservation" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        diff = EegFun.condition_difference(erps, [(1, 2)])

        @test diff[1].sample_rate == erps[1].sample_rate
        @test nrow(diff[1].data) == nrow(erps[1].data)
        @test diff[1].n_epochs == min(erps[1].n_epochs, erps[2].n_epochs)
    end

    @testset "Missing condition skipping" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        # Pair (3, 4) should be skipped, but (1, 2) should succeed
        diff = EegFun.condition_difference(erps, [(1, 2), (3, 4)])

        @test length(diff) == 1  # Only (1, 2) succeeded
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "All conditions missing" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        # No valid pairs → should throw
        @test_throws Exception EegFun.condition_difference(erps, [(5, 6)])
    end

    @testset "Non-contiguous .condition numbers matched by .condition field" begin
        e2 = EegFun.create_test_erp_data(condition = 2)
        e5 = EegFun.create_test_erp_data(condition = 5)
        e2.data.Ch1 .= 30.0
        e5.data.Ch1 .= 10.0
        diff = EegFun.condition_difference([e2, e5], [(2, 5)])
        @test length(diff) == 1
        @test diff[1].data.Ch1[1] ≈ 20.0
    end
end

@testset "TF Condition Difference" begin

    @testset "Basic difference" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        @test length(diff) == 1
        @test diff[1].condition == 1  # pair index
        @test diff[1].condition_name == "difference_1_2"
        @test diff[1] isa EegFun.TimeFreqData
    end

    @testset "Multiple difference pairs" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1),
            EegFun.create_test_tf_data(condition = 2),
            EegFun.create_test_tf_data(condition = 3),
            EegFun.create_test_tf_data(condition = 4),
        ]

        diff = EegFun.condition_difference(tfs, [(1, 2), (3, 4)])

        @test length(diff) == 2
        @test diff[1].condition_name == "difference_1_2"
        @test diff[2].condition_name == "difference_3_4"
    end

    @testset "Vector syntax" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [[1, 2]])

        @test length(diff) == 1
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "Power data integrity" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1, power_offset = 6.0, noise = 0.0),
            EegFun.create_test_tf_data(condition = 2, power_offset = 4.0, noise = 0.0),
        ]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        # Difference of 6.0 - 4.0 should be 2.0
        for ch in [:Ch1, :Ch2, :Ch3]
            @test all(diff[1].data_power[!, ch] .≈ 2.0)
        end
    end

    @testset "Phase data integrity" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        # Phase should be subtracted too
        for ch in [:Ch1, :Ch2, :Ch3]
            expected_phase = tfs[1].data_phase[!, ch] .- tfs[2].data_phase[!, ch]
            @test all(abs.(diff[1].data_phase[!, ch] .- expected_phase) .< 1e-10)
        end
    end

    @testset "Metadata preservation" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        @test diff[1].sample_rate == tfs[1].sample_rate
        @test diff[1].method == tfs[1].method
        @test nrow(diff[1].data_power) == nrow(tfs[1].data_power)
        @test nrow(diff[1].data_phase) == nrow(tfs[1].data_phase)

        # time and freq columns should be preserved
        @test diff[1].data_power.time == tfs[1].data_power.time
        @test diff[1].data_power.freq == tfs[1].data_power.freq
    end

    @testset "Missing condition skipping" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # Pair (3, 4) should be skipped, but (1, 2) should succeed
        diff = EegFun.condition_difference(tfs, [(1, 2), (3, 4)])

        @test length(diff) == 1
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "All conditions missing" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # No valid pairs → should throw
        @test_throws Exception EegFun.condition_difference(tfs, [(5, 6)])
    end
end

@testset "Batch TF Condition Difference" begin
    test_dir = mktempdir()

    @testset "Basic batch TF difference" begin
        # Create test TF files
        for participant = 1:2
            tfs = [
                EegFun.create_test_tf_data(participant = participant, condition = 1),
                EegFun.create_test_tf_data(participant = participant, condition = 2),
                EegFun.create_test_tf_data(participant = participant, condition = 3),
                EegFun.create_test_tf_data(participant = participant, condition = 4),
            ]

            file_path = joinpath(test_dir, "$(participant)_tf_morlet.jld2")
            jldsave(file_path; data = tfs)
        end

        output_dir = joinpath(test_dir, "tf_differences")

        result = EegFun.condition_difference("tf_morlet", [(1, 2), (3, 4)], input_dir = test_dir, output_dir = output_dir)

        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "1_tf_morlet.jld2" in output_files
        @test "2_tf_morlet.jld2" in output_files

        # Load and verify output
        differences = load(joinpath(output_dir, "1_tf_morlet.jld2"), "data")
        @test length(differences) == 2

        # Verify data type is TimeFreqData
        @test differences[1] isa EegFun.TimeFreqData
        @test differences[2] isa EegFun.TimeFreqData

        # Verify naming
        @test differences[1].condition_name == "difference_1_2"
        @test differences[2].condition_name == "difference_3_4"

        # Verify both power and phase data exist
        @test nrow(differences[1].data_power) > 0
        @test nrow(differences[1].data_phase) > 0
    end

    @testset "No matching files" begin
        output_dir = joinpath(test_dir, "tf_diff_none")

        result = EegFun.condition_difference("nonexistent_pattern", [(1, 2)], input_dir = test_dir, output_dir = output_dir)
        @test isnothing(result)
    end

    @testset "Multi-participant in-memory Vector{ErpData}" begin
        mp = [
            EegFun.create_test_erp_data(participant = 1, condition = 1),
            EegFun.create_test_erp_data(participant = 1, condition = 2),
            EegFun.create_test_erp_data(participant = 2, condition = 1),
            EegFun.create_test_erp_data(participant = 2, condition = 2),
        ]
        diffs = EegFun.condition_difference(mp, [(1, 2)])
        @test length(diffs) == 2
        @test [d.file for d in diffs] == [mp[1].file, mp[3].file]

        # Derived extra float channels (:vEOG / channel_difference!) are subtracted, single-pair batch works, invalid pairs error
        d1 = EegFun.create_test_erp_data(participant = 1, condition = 1)
        d2 = EegFun.create_test_erp_data(participant = 1, condition = 2)
        d1.data[!, :vEOG] = fill(25.0, nrow(d1.data))
        d2.data[!, :vEOG] = fill(10.0, nrow(d2.data))
        diff_ex = EegFun.condition_difference([d1, d2], (1, 2))
        @test all(≈(15.0), diff_ex[1].data.vEOG)
        @test_throws EegFun.EegFunError EegFun.condition_difference([d1, d2], [[1]])
        @test_throws EegFun.EegFunError EegFun.condition_difference([d1, d2], [[1, 2, 3]])

        mktempdir() do tmp_sp
            jldsave(joinpath(tmp_sp, "1_erps.jld2"); data = [d1, d2])
            res_sp = EegFun.condition_difference("erps", (1, 2); input_dir = tmp_sp)
            @test res_sp.success == 1
        end
    end

    rm(test_dir, recursive = true)
end
