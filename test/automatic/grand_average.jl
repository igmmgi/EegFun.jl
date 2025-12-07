"""
Test suite for src/analysis/grand_average.jl
"""

using Test
using JLD2
using DataFrames


@testset "Batch Grand Average" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic grand averaging" begin
        # Create test ERP files for multiple participants
        for participant = 1:4
            erps = [
                create_test_erp_data(participant, 1),
                create_test_erp_data(participant, 2),
                create_test_erp_data(participant, 3),
            ]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "grand_average")

        # Test basic grand averaging
        result = eegfun.grand_average("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

        # Verify output file was created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "grand_average_erps_cleaned.jld2" in output_files

        # Load and verify grand averages
        grand_avg_file = joinpath(output_dir, "grand_average_erps_cleaned.jld2")
        @test isfile(grand_avg_file)

        grand_averages = load(grand_avg_file, "data")
        @test length(grand_averages) == 3  # One for each condition

        # Verify grand average structure
        for grand_avg in grand_averages
            condition_num = grand_avg.condition
            @test condition_num in [1, 2, 3]
            @test grand_avg.condition_name == "grand_avg_condition_$condition_num"
            @test grand_avg.sample_rate == 1000.0
            @test nrow(grand_avg.data) == 2501  # n_timepoints
        end
    end

    @testset "Participant filtering" begin
        output_dir = joinpath(test_dir, "grand_average_filtered")

        # Test with specific participants
        result =
            eegfun.grand_average("erps_cleaned", input_dir = test_dir, participant_selection = eegfun.participants([1, 2, 3]), output_dir = output_dir)

        @test isdir(output_dir)
        grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")
        @test length(grand_averages) == 3

        # Verify that only 3 participants were used (check n_epochs)
        for grand_avg in grand_averages
            @test grand_avg.n_epochs == 30  # 3 participants × 10 epochs each
        end
    end

    @testset "Condition filtering" begin
        output_dir = joinpath(test_dir, "grand_average_conditions")

        # Test with specific conditions
        result = eegfun.grand_average("erps_cleaned", input_dir = test_dir, condition_selection = eegfun.conditions([1, 2]), output_dir = output_dir)

        @test isdir(output_dir)
        grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")
        @test length(grand_averages) == 2  # Only conditions 1 and 2

        # Verify condition numbers
        condition_nums = [grand_avg.condition for grand_avg in grand_averages]
        @test 1 in condition_nums
        @test 2 in condition_nums
    end

    @testset "Grand average calculation verification" begin
        output_dir = joinpath(test_dir, "grand_average_verification")

        result = eegfun.grand_average(
            "erps_cleaned",
            input_dir = test_dir,
            participant_selection = eegfun.participants([1, 2]),
            condition_selection = eegfun.conditions([1]),
            output_dir = output_dir,
        )

        # Load original data for verification
        erp1_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        erp1_data = load(erp1_file, "data")
        erp1_cond1 = erp1_data[1]  # Condition 1 from participant 1

        erp2_file = joinpath(test_dir, "2_erps_cleaned.jld2")
        erp2_data = load(erp2_file, "data")
        erp2_cond1 = erp2_data[1]  # Condition 1 from participant 2

        # Load grand average
        grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")
        grand_avg_cond1 = grand_averages[1]

        # Verify grand average is actually the mean of the two ERPs
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(erp1_cond1.data, ch) &&
               hasproperty(erp2_cond1.data, ch) &&
               hasproperty(grand_avg_cond1.data, ch)
                expected_grand_avg = (erp1_cond1.data[!, ch] .+ erp2_cond1.data[!, ch]) ./ 2
                @test all(abs.(grand_avg_cond1.data[!, ch] .- expected_grand_avg) .< 1e-10)
            end
        end

        # Verify total epochs
        @test grand_avg_cond1.n_epochs == erp1_cond1.n_epochs + erp2_cond1.n_epochs
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.grand_average("erps_cleaned", input_dir = "/nonexistent/dir")
        end

        @testset "Non-ERP pattern" begin
            @test_throws Exception eegfun.grand_average("epochs_cleaned", input_dir = test_dir)
        end

        @testset "No matching files" begin
            @test_throws Exception eegfun.grand_average("nonexistent_pattern", input_dir = test_dir)
        end

        @testset "Files with no ERP data" begin
            # Create file with no 'erps' variable
            no_erps_file = joinpath(test_dir, "no_erps_cleaned.jld2")
            jldsave(no_erps_file; other_data = "test")

            output_dir = joinpath(test_dir, "grand_average_no_erps")

            result = eegfun.grand_average(
                "erps_cleaned",
                input_dir = test_dir,
                participant_selection = eegfun.participants(999),  # Non-existent participant
                output_dir = output_dir,
            )

            @test result === nothing
        end
    end

    @testset "Edge cases" begin
        @testset "Single participant (should skip grand average)" begin
            output_dir = joinpath(test_dir, "grand_average_single")

            result =
                eegfun.grand_average("erps_cleaned", input_dir = test_dir, participant_selection = eegfun.participants([1]), output_dir = output_dir)

            @test result === nothing
        end

        @testset "Insufficient participants for some conditions" begin
            # Create file with only some conditions for one participant
            erps = [
                create_test_erp_data(5, 1),
                create_test_erp_data(5, 2),
                # Missing condition 3
            ]

            file_path = joinpath(test_dir, "5_erps_cleaned.jld2")
            jldsave(file_path; data = erps)

            output_dir = joinpath(test_dir, "grand_average_insufficient")

            result = eegfun.grand_average(
                "erps_cleaned",
                input_dir = test_dir,
                participant_selection = eegfun.participants([1, 2, 5]),
                output_dir = output_dir,
            )

            # Should create grand averages for all conditions
            @test isdir(output_dir)
            grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")
            @test length(grand_averages) == 3  # All conditions
        end

        @testset "Empty ERP data" begin
            # Create file with empty ERP list
            empty_file = joinpath(test_dir, "empty_erps_cleaned.jld2")
            jldsave(empty_file; data = eegfun.ErpData[])

            output_dir = joinpath(test_dir, "grand_average_empty")

            result = eegfun.grand_average(
                "erps_cleaned",
                input_dir = test_dir,
                participant_selection = eegfun.participants(999),  # Non-existent participant
                output_dir = output_dir,
            )

            @test result === nothing
        end
    end

    @testset "Data structure validation" begin
        output_dir = joinpath(test_dir, "grand_average_structure")

        result = eegfun.grand_average(
            "erps_cleaned",
            input_dir = test_dir,
            participant_selection = eegfun.participants([1, 2]),
            condition_selection = eegfun.conditions([1]),
            output_dir = output_dir,
        )

        grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")
        grand_avg = grand_averages[1]

        # Verify ErpData structure
        @test grand_avg isa eegfun.ErpData
        @test grand_avg.sample_rate == 1000.0
        @test grand_avg.layout isa eegfun.Layout
        @test grand_avg.analysis_info isa eegfun.AnalysisInfo

        # Verify DataFrame structure
        @test grand_avg.data isa DataFrame
        @test nrow(grand_avg.data) == 2501  # n_timepoints
        @test "time" in names(grand_avg.data)
        # condition and condition_name are now struct fields, not DataFrame columns
        @test "condition" ∉ names(grand_avg.data)
        @test "condition_name" ∉ names(grand_avg.data)
        @test "Ch1" in names(grand_avg.data)
        @test "Ch2" in names(grand_avg.data)
        @test "Ch3" in names(grand_avg.data)

        # Verify metadata (now in struct fields)
        @test grand_avg.condition == 1
        @test grand_avg.condition_name == "grand_avg_condition_1"
        @test hasproperty(grand_avg, :condition)
        @test hasproperty(grand_avg, :condition_name)
    end

    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_grand_average")

            result = eegfun.grand_average("erps_cleaned", input_dir = test_dir, output_dir = custom_dir)

            @test isdir(custom_dir)
            @test "grand_average_erps_cleaned.jld2" in readdir(custom_dir)
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.grand_average("erps_cleaned", input_dir = test_dir)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "grand_average_erps_cleaned")
            @test isdir(expected_dir)
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "grand_average_logging")

        result = eegfun.grand_average("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

        # Check that log file was created
        log_file = joinpath(output_dir, "grand_average.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch grand averaging started", log_content)
        @test occursin("Found 7 JLD2 files", log_content)
        @test occursin("Grand averaging complete", log_content)
    end

    @testset "Multiple conditions with different participant counts" begin
        # Use a separate directory to avoid matching files from earlier tests
        separate_test_dir = joinpath(test_dir, "different_counts")
        mkpath(separate_test_dir)
        
        # Create scenario where some conditions have more participants than others
        # 3 participants have condition 1
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1)]
            file_path = joinpath(separate_test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        # Only 2 participants have condition 2 (append to existing files)
        for participant = 1:2
            erp2 = create_test_erp_data(participant, 2)
            file_path = joinpath(separate_test_dir, "$(participant)_erps_cleaned.jld2")
            # Load existing data and append condition 2
            existing_data = load(file_path, "data")
            push!(existing_data, erp2)
            jldsave(file_path; data = existing_data)
        end

        output_dir = joinpath(separate_test_dir, "grand_average_different_counts")

        result = eegfun.grand_average("erps_cleaned", input_dir = separate_test_dir, output_dir = output_dir)

        @test isdir(output_dir)
        grand_averages = load(joinpath(output_dir, "grand_average_erps_cleaned.jld2"), "data")

        # Should have grand averages for conditions 1 and 2
        @test length(grand_averages) == 2
        # Condition 1: 3 participants × 10 epochs = 30 epochs
        # Condition 2: 2 participants × 10 epochs = 20 epochs
        # Find which condition is which (they're sorted by condition number)
        sorted_grand_averages = sort(grand_averages, by = x -> x.condition)
        @test sorted_grand_averages[1].condition == 1
        @test sorted_grand_averages[1].n_epochs == 30
        @test sorted_grand_averages[2].condition == 2
        @test sorted_grand_averages[2].n_epochs == 20
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
