using Test
using JLD2
using DataFrames

@testset "Batch Average Conditions" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic average wave creation" begin
        # Create test ERP files
        for participant = 1:3
            erps = [
                EegFun.create_test_erp_data(participant = participant, condition = 1),
                EegFun.create_test_erp_data(participant = participant, condition = 2),
                EegFun.create_test_erp_data(participant = participant, condition = 3),
                EegFun.create_test_erp_data(participant = participant, condition = 4),
            ]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "averages")

        # Test basic average creation
        result = EegFun.condition_average("erps_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) >= 3
        @test "1_erps_cleaned.jld2" in output_files
        @test "2_erps_cleaned.jld2" in output_files
        @test "3_erps_cleaned.jld2" in output_files

        # Load and verify one output file
        output_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        @test isfile(output_file)

        averages = load(output_file, "data")
        @test length(averages) == 2  # Two average waves

        # Verify average wave structure
        avg1 = averages[1]
        @test avg1.condition == 1  # First average labeled as condition 1
        @test avg1.condition_name == "average_1_2"

        avg2 = averages[2]
        @test avg2.condition == 2  # Second average labeled as condition 2
        @test avg2.condition_name == "average_3_4"
    end

    @testset "Single participant processing" begin
        output_dir = joinpath(test_dir, "averages_single")

        result = EegFun.condition_average(
            "erps_cleaned",
            [[1, 2]],
            input_dir = test_dir,
            participant_selection = EegFun.participants(2),
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "2_erps_cleaned.jld2" in output_files
        @test length(output_files) >= 1
    end

    @testset "Missing conditions handling" begin
        # Create file with only some conditions
        erps = [
            EegFun.create_test_erp_data(participant = 99, condition = 1),
            EegFun.create_test_erp_data(participant = 99, condition = 2),
            # Missing conditions 3 and 4
        ]

        file_path = joinpath(test_dir, "99_erps_cleaned.jld2")
        jldsave(file_path; data = erps)

        output_dir = joinpath(test_dir, "averages_missing")

        result = EegFun.condition_average(
            "erps_cleaned",
            [[1, 2], [3, 4]],
            input_dir = test_dir,
            participant_selection = EegFun.participants(99),
            output_dir = output_dir,
        )

        # Should still create file but only with available groups
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "99_erps_cleaned.jld2" in output_files
        @test length(output_files) >= 1

        averages = load(joinpath(output_dir, "99_erps_cleaned.jld2"), "data")
        @test length(averages) == 1  # Only one average wave (1, 2)
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception EegFun.condition_average("erps_cleaned", [[1, 2]], input_dir = "/nonexistent/dir")
        end

        @testset "Non-ERP pattern" begin
            @test_throws Exception EegFun.condition_average("epochs_cleaned", [[1, 2]], input_dir = test_dir)
        end

        @testset "Empty condition groups" begin
            @test_throws Exception EegFun.condition_average("erps_cleaned", Vector{Int}[], input_dir = test_dir)
        end
    end

    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_avg_output")

            result = EegFun.condition_average("erps_cleaned", [[1, 2]], input_dir = test_dir, output_dir = custom_dir)

            @test isdir(custom_dir)
            @test length(readdir(custom_dir)) >= 3
        end

        @testset "Auto-generated output directory" begin
            result = EegFun.condition_average("erps_cleaned", [[1, 2]], input_dir = test_dir)

            expected_dir = joinpath(test_dir, "averages_erps_cleaned_1-2")
            @test isdir(expected_dir)
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "averages_logging")

        result = EegFun.condition_average("erps_cleaned", [[1, 2]], input_dir = test_dir, output_dir = output_dir)

        # Check that log file was created
        log_file = joinpath(output_dir, "condition_average.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch condition averaging started", log_content)
        @test occursin("Found", log_content)
        @test occursin("Batch operation complete", log_content)
    end

    # Cleanup
    rm(test_dir, recursive = true)
end

@testset "In-Memory Average Conditions" begin

    @testset "Basic average wave" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        avg = EegFun.condition_average(erps, [[1, 2]])

        @test length(avg) == 1
        @test avg[1].condition == 1  # group index
        @test avg[1].condition_name == "average_1_2"
    end

    @testset "Multiple average groups" begin
        erps = [
            EegFun.create_test_erp_data(condition = 1),
            EegFun.create_test_erp_data(condition = 2),
            EegFun.create_test_erp_data(condition = 3),
            EegFun.create_test_erp_data(condition = 4),
        ]

        avg = EegFun.condition_average(erps, [[1, 2], [3, 4]])

        @test length(avg) == 2
        @test avg[1].condition_name == "average_1_2"
        @test avg[2].condition_name == "average_3_4"
    end

    @testset "Average all conditions" begin
        erps = [
            EegFun.create_test_erp_data(condition = 1),
            EegFun.create_test_erp_data(condition = 2),
            EegFun.create_test_erp_data(condition = 3),
        ]

        avg = EegFun.condition_average(erps, [[1, 2, 3]])

        @test length(avg) == 1
        @test avg[1].condition_name == "average_1_2_3"
    end

    @testset "Data integrity" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        avg = EegFun.condition_average(erps, [[1, 2]])

        # Verify average is actually (erp1 + erp2) / 2 for EEG channels
        for ch in [:Fz, :Cz, :Pz]
            if hasproperty(erps[1].data, ch) && hasproperty(erps[2].data, ch) && hasproperty(avg[1].data, ch)
                expected = (erps[1].data[!, ch] .+ erps[2].data[!, ch]) ./ 2
                @test all(abs.(avg[1].data[!, ch] .- expected) .< 1e-10)
            end
        end
    end

    @testset "Metadata preservation" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        avg = EegFun.condition_average(erps, [[1, 2]])

        @test avg[1].sample_rate == erps[1].sample_rate
        @test nrow(avg[1].data) == nrow(erps[1].data)
        @test avg[1].n_epochs == erps[1].n_epochs + erps[2].n_epochs
    end

    @testset "Missing condition skipping" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        # Group [3, 4] should be skipped, but [1, 2] should succeed
        avg = EegFun.condition_average(erps, [[1, 2], [3, 4]])

        @test length(avg) == 1  # Only [1, 2] succeeded
        @test avg[1].condition_name == "average_1_2"
    end

    @testset "All conditions missing" begin
        erps = [EegFun.create_test_erp_data(condition = 1), EegFun.create_test_erp_data(condition = 2)]

        # No valid groups → should throw
        @test_throws Exception EegFun.condition_average(erps, [[5, 6]])
    end
end
