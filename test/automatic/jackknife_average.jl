"""
Test suite for src/analysis/jackknife_average.jl
"""

using Test
using JLD2
using DataFrames


# Use generic create_test_erp_data from test_utils.jl
# create_test_erp_data(participant, condition, n_timepoints, n_channels)

@testset "Jackknife Average" begin
    @testset "In-memory jackknife_average function" begin
        @testset "Basic jackknife averaging" begin
            # Create test ERP data for 4 participants
            erps = [create_test_erp_data(i, 1) for i = 1:4]

            # Create jackknife averages
            jackknife_results = eegfun.jackknife_average(erps)

            # Should have 4 jackknife averages (one per participant)
            @test length(jackknife_results) == 4

            # Each jackknife should be an ErpData object
            for jk in jackknife_results
                @test jk isa eegfun.ErpData
                @test jk.sample_rate == 1000.0
                @test nrow(jk.data) == 2501
            end
        end

        @testset "Jackknife calculation verification" begin
            # Create simple test data with known values
            erps = [create_test_erp_data(i, 1, fs = 50) for i = 1:3]

            jackknife_results = eegfun.jackknife_average(erps)

            # Jackknife 1 should be average of participants 2 and 3
            # Jackknife 2 should be average of participants 1 and 3
            # Jackknife 3 should be average of participants 1 and 2

            # Verify first jackknife (excluding participant 1)
            jk1 = jackknife_results[1]
            expected_ch = (erps[2].data.Ch1 .+ erps[3].data.Ch1) ./ 2
            @test all(abs.(jk1.data.Ch1 .- expected_ch) .< 1e-10)

            # Verify second jackknife (excluding participant 2)
            jk2 = jackknife_results[2]
            expected_ch = (erps[1].data.Ch1 .+ erps[3].data.Ch1) ./ 2
            @test all(abs.(jk2.data.Ch1 .- expected_ch) .< 1e-10)

            # Verify third jackknife (excluding participant 3)
            jk3 = jackknife_results[3]
            expected_ch = (erps[1].data.Ch1 .+ erps[2].data.Ch1) ./ 2
            @test all(abs.(jk3.data.Ch1 .- expected_ch) .< 1e-10)
        end

        @testset "Jackknife with different channels" begin
            # Test that all channels are processed correctly
            erps = [create_test_erp_data(i, 1, fs = 50, n_channels = 3) for i = 1:3]

            jackknife_results = eegfun.jackknife_average(erps)

            # Verify all channels are present in jackknife results
            for jk in jackknife_results
                @test hasproperty(jk.data, :Ch1)
                @test hasproperty(jk.data, :Ch2)
                @test hasproperty(jk.data, :Ch3)
            end
        end

        @testset "Error handling: insufficient participants" begin
            # Only 1 participant - should throw error
            erps = [create_test_erp_data(1, 1)]
            @test_throws Exception eegfun.jackknife_average(erps)
        end

        @testset "Error handling: mismatched sample rates" begin
            # Create ERPs with different sample rates
            erp1 = create_test_erp_data(1, 1)
            erp2 = create_test_erp_data(2, 1)
            erp2 = eegfun.ErpData(erp2.data, erp2.layout, 500.0, erp2.analysis_info, erp2.n_epochs)  # Different sample rate

            @test_throws Exception eegfun.jackknife_average([erp1, erp2])
        end

        @testset "Error handling: mismatched time points" begin
            # Create ERPs with different numbers of time points
            erp1 = create_test_erp_data(1, 1, fs = 100)
            erp2 = create_test_erp_data(2, 1, fs = 50)

            @test_throws Exception eegfun.jackknife_average([erp1, erp2])
        end

        @testset "Metadata preservation" begin
            erps = [create_test_erp_data(i, 1) for i = 1:3]

            jackknife_results = eegfun.jackknife_average(erps)

            # Check that metadata is preserved
            for (i, jk) in enumerate(jackknife_results)
                @test hasproperty(jk.data, :time)
                @test hasproperty(jk.data, :condition)
                @test hasproperty(jk.data, :condition_name)
                @test occursin("jackknife", jk.data.condition_name[1])
            end
        end

        @testset "n_epochs calculation" begin
            # Create ERPs with specific n_epochs
            erps = [create_test_erp_data(i, 1) for i = 1:4]
            for (i, erp) in enumerate(erps)
                erps[i] = eegfun.ErpData(erp.data, erp.layout, erp.sample_rate, erp.analysis_info, i * 10)
            end

            jackknife_results = eegfun.jackknife_average(erps)

            # Jackknife 1 excludes participant 1 (10 epochs), includes 2 (20) + 3 (30) + 4 (40) = 90
            @test jackknife_results[1].n_epochs == 90

            # Jackknife 2 excludes participant 2 (20 epochs), includes 1 (10) + 3 (30) + 4 (40) = 80
            @test jackknife_results[2].n_epochs == 80

            # Jackknife 3 excludes participant 3 (30 epochs), includes 1 (10) + 2 (20) + 4 (40) = 70
            @test jackknife_results[3].n_epochs == 70

            # Jackknife 4 excludes participant 4 (40 epochs), includes 1 (10) + 2 (20) + 3 (30) = 60
            @test jackknife_results[4].n_epochs == 60
        end
    end

    # Create temporary test directory for batch processing tests
    test_dir = mktempdir()

    @testset "Batch jackknife_average function" begin
        @testset "Basic batch processing" begin
            # Create test LRP files for multiple participants
            for participant = 1:4
                lrp_data = create_test_erp_data(participant, 1)
                file_path = joinpath(test_dir, "$(participant)_lrp.jld2")
                save(file_path, "lrp", lrp_data)
            end

            output_dir = joinpath(test_dir, "jackknife_test")

            # Test basic jackknife averaging
            result = eegfun.jackknife_average("lrp", input_dir = test_dir, output_dir = output_dir)

            # Verify output directory was created
            @test isdir(output_dir)

            # Verify output files (one per participant)
            output_files = readdir(output_dir)
            @test "1_lrp.jld2" in output_files
            @test "2_lrp.jld2" in output_files
            @test "3_lrp.jld2" in output_files
            @test "4_lrp.jld2" in output_files

            # Load and verify jackknife data
            jk1 = load(joinpath(output_dir, "1_lrp.jld2"), "jackknife")
            @test jk1 isa eegfun.ErpData
            @test nrow(jk1.data) == 2501
        end

        @testset "Multiple conditions" begin
            # Create test LRP files with multiple conditions
            for participant = 1:4
                lrp_data = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]
                file_path = joinpath(test_dir, "$(participant)_multi_lrp.jld2")
                save(file_path, "lrp", lrp_data)
            end

            output_dir = joinpath(test_dir, "jackknife_multi")

            result = eegfun.jackknife_average("multi_lrp", input_dir = test_dir, output_dir = output_dir)

            @test isdir(output_dir)

            # Load and verify - should have vector of ErpData for multiple conditions
            jk1 = load(joinpath(output_dir, "1_multi_lrp.jld2"), "jackknife")
            @test jk1 isa Vector{eegfun.ErpData}
            @test length(jk1) == 2

            # Verify each condition
            @test jk1[1].data.condition[1] == 1
            @test jk1[2].data.condition[1] == 2
        end

        @testset "Participant filtering" begin
            output_dir = joinpath(test_dir, "jackknife_filtered")

            # Test with specific participants
            # Note: pattern "lrp" matches both "_lrp" and "_multi_lrp" files
            # So we need to check the actual output corresponds to input files
            result =
                eegfun.jackknife_average("lrp", input_dir = test_dir, participants = [1, 2, 3], output_dir = output_dir)

            @test isdir(output_dir)

            # Should have output files for the participants requested
            output_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir))

            # Files with pattern "lrp" include both "_lrp.jld2" and "_multi_lrp.jld2"
            # Check that we have the participants we requested
            @test "1_lrp.jld2" in output_files
            @test "2_lrp.jld2" in output_files
            @test "3_lrp.jld2" in output_files
            @test "1_multi_lrp.jld2" in output_files
            @test "2_multi_lrp.jld2" in output_files
            @test "3_multi_lrp.jld2" in output_files

            # And not participant 4
            @test !("4_lrp.jld2" in output_files)
            @test !("4_multi_lrp.jld2" in output_files)
        end

        @testset "Condition filtering" begin
            output_dir = joinpath(test_dir, "jackknife_cond_filter")

            # Test with specific conditions
            result =
                eegfun.jackknife_average("multi_lrp", input_dir = test_dir, conditions = [1], output_dir = output_dir)

            @test isdir(output_dir)

            # Load and verify - should only have condition 1
            jk1 = load(joinpath(output_dir, "1_multi_lrp.jld2"), "jackknife")
            @test jk1 isa eegfun.ErpData  # Single ErpData, not vector
            @test jk1.data.condition[1] == 1
        end

        @testset "Custom data variable" begin
            # Create test ERP files (not LRP)
            for participant = 1:3
                erp_data = [create_test_erp_data(participant, 1)]
                file_path = joinpath(test_dir, "$(participant)_erps.jld2")
                save(file_path, "erps", erp_data)
            end

            output_dir = joinpath(test_dir, "jackknife_erps")

            # Test with different data variable name
            result = eegfun.jackknife_average("erps", input_dir = test_dir, output_dir = output_dir, data_var = "erps")

            @test isdir(output_dir)

            # Verify output files
            @test "1_erps.jld2" in readdir(output_dir)
            @test "2_erps.jld2" in readdir(output_dir)
            @test "3_erps.jld2" in readdir(output_dir)
        end

        @testset "Jackknife calculation verification in batch" begin
            # Create a separate test directory with only single-condition files to avoid ambiguity
            verify_dir = mktempdir()

            # Create test files with specific known values
            for participant = 1:3
                lrp_data = create_test_erp_data(participant, 1)
                file_path = joinpath(verify_dir, "$(participant)_verify.jld2")
                save(file_path, "lrp", lrp_data)
            end

            output_dir = joinpath(verify_dir, "jackknife_output")

            result = eegfun.jackknife_average("verify", input_dir = verify_dir, output_dir = output_dir)

            # Load original data
            lrp1 = load(joinpath(verify_dir, "1_verify.jld2"), "lrp")
            lrp2 = load(joinpath(verify_dir, "2_verify.jld2"), "lrp")
            lrp3 = load(joinpath(verify_dir, "3_verify.jld2"), "lrp")

            # Load jackknife for participant 1 (should exclude participant 1)
            jk1 = load(joinpath(output_dir, "1_verify.jld2"), "jackknife")

            # Since there's only one condition, should be single ErpData
            @test jk1 isa eegfun.ErpData

            # Verify: jackknife 1 should be average of participants 2 and 3
            expected_ch = (lrp2.data.Ch1 .+ lrp3.data.Ch1) ./ 2
            @test all(abs.(jk1.data.Ch1 .- expected_ch) .< 1e-10)

            # Cleanup
            rm(verify_dir, recursive = true)
        end

        @testset "Error handling: insufficient files" begin
            # Create directory with only 1 file
            single_dir = mktempdir()
            save(joinpath(single_dir, "1_lrp.jld2"), "lrp", create_test_erp_data(1, 1))

            output_dir = joinpath(single_dir, "jackknife_insufficient")

            result = eegfun.jackknife_average("lrp", input_dir = single_dir, output_dir = output_dir)

            # Should return nothing due to insufficient participants
            @test result === nothing

            rm(single_dir, recursive = true)
        end

        @testset "Error handling: no matching files" begin
            output_dir = joinpath(test_dir, "jackknife_nomatch")

            result = eegfun.jackknife_average("nonexistent", input_dir = test_dir, output_dir = output_dir)

            # Should return nothing when no files found
            @test result === nothing
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "jackknife_logging")

            result = eegfun.jackknife_average("lrp", input_dir = test_dir, output_dir = output_dir)

            # Check that log file was created
            log_file = joinpath(output_dir, "jackknife.log")
            @test isfile(log_file)

            # Verify log content
            log_content = read(log_file, String)
            @test occursin("Batch jackknife averaging started", log_content)
            @test occursin("Found", log_content)
            @test occursin("Jackknife averaging complete", log_content)
        end

        @testset "Output structure" begin
            output_dir = joinpath(test_dir, "jackknife_structure")

            result =
                eegfun.jackknife_average("lrp", input_dir = test_dir, participants = [1, 2], output_dir = output_dir)

            jk1 = load(joinpath(output_dir, "1_lrp.jld2"), "jackknife")

            # When multiple conditions exist, result is Vector{ErpData}, otherwise single ErpData
            if jk1 isa Vector
                @test length(jk1) > 0
                @test all(x -> x isa eegfun.ErpData, jk1)

                # Test first condition
                jk1_cond = jk1[1]
                @test jk1_cond.sample_rate == 1000.0
                @test jk1_cond.layout isa eegfun.Layout
                @test jk1_cond.analysis_info isa eegfun.AnalysisInfo
                @test jk1_cond.data isa DataFrame
                @test nrow(jk1_cond.data) == 2501
                @test "time" in names(jk1_cond.data)
                @test "condition" in names(jk1_cond.data)
                @test "condition_name" in names(jk1_cond.data)
                @test "Ch1" in names(jk1_cond.data)
                @test "Ch2" in names(jk1_cond.data)
                @test "Ch3" in names(jk1_cond.data)
            else
                # Single condition case
                @test jk1 isa eegfun.ErpData
                @test jk1.sample_rate == 1000.0
                @test jk1.layout isa eegfun.Layout
                @test jk1.analysis_info isa eegfun.AnalysisInfo
                @test jk1.data isa DataFrame
                @test nrow(jk1.data) == 2501
                @test "time" in names(jk1.data)
                @test "condition" in names(jk1.data)
                @test "condition_name" in names(jk1.data)
                @test "Ch1" in names(jk1.data)
                @test "Ch2" in names(jk1.data)
                @test "Ch3" in names(jk1.data)
            end
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.jackknife_average("lrp", input_dir = test_dir)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "jackknife_lrp")
            @test isdir(expected_dir)
        end
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
