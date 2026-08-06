using Test
using JLD2
using DataFrames
using CSV

@testset "Batch Difference Conditions" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic difference wave creation" begin
        # Create test ERP files
        for participant = 1:3
            erps = [
                EegFun.create_test_erp_data(participant = participant, condition = 1),
                EegFun.create_test_erp_data(participant = participant, condition = 2),
                EegFun.create_test_erp_data(participant = participant, condition = 3),
                EegFun.create_test_erp_data(participant = participant, condition = 4),
            ]

            file_path = joinpath(test_dir, "$(participant)_erps_unrejected.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "differences")

        # Test basic difference creation
        result = EegFun.condition_difference("erps_unrejected", [(1, 2), (3, 4)], input_dir = test_dir, output_dir = output_dir)

        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        # Should have 3 files: 1, 2, 3 (participants from this test)
        @test length(output_files) >= 3  # At least 3 files from this test
        @test "1_erps_unrejected.jld2" in output_files
        @test "2_erps_unrejected.jld2" in output_files
        @test "3_erps_unrejected.jld2" in output_files

        # Load and verify one output file
        output_file = joinpath(output_dir, "1_erps_unrejected.jld2")
        @test isfile(output_file)

        differences = load(output_file, "data")
        @test length(differences) == 2  # Two difference waves

        # Verify difference wave structure
        diff1 = differences[1]
        @test diff1.condition == 1  # First difference labeled as condition 1
        @test diff1.condition_name == "difference_1_2"

        diff2 = differences[2]
        @test diff2.condition == 2  # Second difference labeled as condition 2
        @test diff2.condition_name == "difference_3_4"
    end

    @testset "Single participant processing" begin
        output_dir = joinpath(test_dir, "differences_single")

        result = EegFun.condition_difference(
            "erps_unrejected",
            [(1, 2)],
            input_dir = test_dir,
            participant_selection = EegFun.participants(2),
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "2_erps_unrejected.jld2" in output_files
        # Should have at least 1 file (participant 2), but might have more from other tests
        @test length(output_files) >= 1
    end

    @testset "Vector condition pairs" begin
        output_dir = joinpath(test_dir, "differences_vector")

        result = EegFun.condition_difference("erps_unrejected", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        # Should have at least 3 files (participants 1, 2, 3), but might have more from other tests
        @test length(output_files) >= 3
    end

    @testset "Missing conditions handling" begin
        # Create file with only some conditions
        erps = [
            EegFun.create_test_erp_data(participant = 99, condition = 1),
            EegFun.create_test_erp_data(participant = 99, condition = 2),
            # Missing conditions 3 and 4
        ]

        file_path = joinpath(test_dir, "99_erps_unrejected.jld2")
        jldsave(file_path; data = erps)

        output_dir = joinpath(test_dir, "differences_missing")

        result = EegFun.condition_difference(
            "erps_unrejected",
            [(1, 2), (3, 4)],
            input_dir = test_dir,
            participant_selection = EegFun.participants(99),
            output_dir = output_dir,
        )

        # Should still create file but only with available pairs
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "99_erps_unrejected.jld2" in output_files
        # Should have at least 1 file (participant 99), but might have more from other tests
        @test length(output_files) >= 1

        differences = load(joinpath(output_dir, "99_erps_unrejected.jld2"), "data")
        @test length(differences) == 1  # Only one difference wave (1-2)
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception EegFun.condition_difference("erps_unrejected", [(1, 2)], input_dir = "/nonexistent/dir")
        end




        @testset "Empty condition pairs" begin
            @test_throws Exception EegFun.condition_difference("erps_unrejected", [], input_dir = test_dir)
        end

        @testset "Invalid condition pairs" begin
            @test_throws Exception EegFun.condition_difference("erps_unrejected", [(1, "invalid")], input_dir = test_dir)
        end
    end

    @testset "Data integrity" begin
        output_dir = joinpath(test_dir, "differences_integrity")

        result = EegFun.condition_difference(
            "erps_unrejected",
            [(1, 2)],
            input_dir = test_dir,
            participant_selection = EegFun.participants(1),
            output_dir = output_dir,
        )

        # Load original and difference data
        original_file = joinpath(test_dir, "1_erps_unrejected.jld2")
        original_erps = load(original_file, "data")

        diff_file = joinpath(output_dir, "1_erps_unrejected.jld2")
        differences = load(diff_file, "data")

        # Verify difference calculation
        erp1 = original_erps[1]  # Condition 1
        erp2 = original_erps[2]  # Condition 2
        diff = differences[1]    # Difference 1-2

        # Check that difference is actually erp1 - erp2
        for ch in [:Fz, :Cz, :Pz]
            if hasproperty(erp1.data, ch) && hasproperty(erp2.data, ch) && hasproperty(diff.data, ch)
                expected_diff = erp1.data[!, ch] .- erp2.data[!, ch]
                @test all(abs.(diff.data[!, ch] .- expected_diff) .< 1e-10)
            end
        end

        # Verify metadata preservation
        @test diff.sample_rate == erp1.sample_rate
        @test nrow(diff.data) == nrow(erp1.data)
        @test diff.n_epochs == min(erp1.n_epochs, erp2.n_epochs)
    end

    @testset "Edge cases" begin

        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "differences_none")

            # No files match so the function returns nothing
            result = EegFun.condition_difference("nonexistent_pattern", [(1, 2)], input_dir = test_dir, output_dir = output_dir)
            @test isnothing(result)
        end

        @testset "Empty ERP data" begin
            # Create file with empty ERP list
            empty_file = joinpath(test_dir, "empty_erps_unrejected.jld2")
            jldsave(empty_file; data = EegFun.ErpData[])

            output_dir = joinpath(test_dir, "differences_empty")

            result = EegFun.condition_difference(
                "erps_unrejected",
                [(1, 2)],
                input_dir = test_dir,
                participant_selection = EegFun.participants(999),  # Non-existent participant
                output_dir = output_dir,
            )

            # The function returns nothing when no files match
            @test isnothing(result)
        end
    end

    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_output")

            result = EegFun.condition_difference("erps_unrejected", [(1, 2)], input_dir = test_dir, output_dir = custom_dir)

            @test isdir(custom_dir)
            # Expect 5 files: 1, 2, 3, 99, and empty (but empty will have 0 differences)
            @test length(readdir(custom_dir)) == 5
        end

        @testset "Auto-generated output directory" begin
            result = EegFun.condition_difference("erps_unrejected", [(1, 2)], input_dir = test_dir)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "differences_erps_unrejected_1-2")
            @test isdir(expected_dir)
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "differences_logging")

        result = EegFun.condition_difference("erps_unrejected", [(1, 2)], input_dir = test_dir, output_dir = output_dir)

        # Check that log file was created
        log_file = joinpath(output_dir, "condition_difference.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch condition differencing started", log_content)
        @test occursin("Found 5 JLD2 files", log_content)
        @test occursin("Batch operation complete", log_content)
    end

    # Cleanup
    rm(test_dir, recursive = true)
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

    rm(test_dir, recursive = true)
end
