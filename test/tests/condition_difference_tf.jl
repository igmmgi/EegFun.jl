using Test
using DataFrames
using JLD2

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
        @test result.success == 0
    end

    rm(test_dir, recursive = true)
end
