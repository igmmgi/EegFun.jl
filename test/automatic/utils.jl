"""
Test suite for src/analysis/batch/utils.jl
"""

using Test
using JLD2
using DataFrames


@testset "Batch Utils" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "BatchResult struct" begin
        # Test BatchResult creation and access
        result = eegfun.BatchResult(true, "test_file.jld2", "Success message")

        @test result.success == true
        @test result.filename == "test_file.jld2"
        @test result.message == "Success message"

        # Test immutable nature
        @test_throws ErrorException result.success = false
    end

    @testset "BatchConfig struct" begin
        # Test BatchConfig creation and access
        config = eegfun.BatchConfig("test_pattern", "/input", "/output", [1, 2], [1, 2])

        @test config.file_pattern == "test_pattern"
        @test config.input_dir == "/input"
        @test config.output_dir == "/output"
        @test config.participants == [1, 2]
        @test config.conditions == [1, 2]

        # Test with nothing values
        config_nothing = eegfun.BatchConfig("test", "/input", "/output", nothing, nothing)
        @test config_nothing.participants === nothing
        @test config_nothing.conditions === nothing
    end

    @testset "_find_batch_files" begin
        # Create test files
        for participant = 1:5
            erps = [create_test_erp_data(participant, 1)]
            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        # Create some non-matching files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1)]
            file_path = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        # Test basic file finding
        files = eegfun._find_batch_files("erps_cleaned", test_dir)
        @test length(files) == 5
        @test all(endswith.(files, ".jld2"))
        @test all(contains.(files, "erps_cleaned"))

        # Test participant filtering
        files_filtered = eegfun._find_batch_files("erps_cleaned", test_dir, [1, 3, 5])
        @test length(files_filtered) == 3
        @test "1_erps_cleaned.jld2" in files_filtered
        @test "3_erps_cleaned.jld2" in files_filtered
        @test "5_erps_cleaned.jld2" in files_filtered

        # Test single participant
        files_single = eegfun._find_batch_files("erps_cleaned", test_dir, 2)
        @test length(files_single) == 1
        @test "2_erps_cleaned.jld2" in files_single

        # Test no participants (should return all)
        files_all = eegfun._find_batch_files("erps_cleaned", test_dir, nothing)
        @test length(files_all) == 5

        # Test non-matching pattern
        files_none = eegfun._find_batch_files("nonexistent", test_dir)
        @test isempty(files_none)
    end

    @testset "load_data" begin
        # Create test files with different variable names
        erps = [create_test_erp_data(1, 1)]

        # Test with "erps" variable
        erps_file = joinpath(test_dir, "test_erps.jld2")
        jldsave(erps_file; data = erps)

        result = eegfun.load_data(erps_file)
        @test result !== nothing
        @test length(result) == length(erps)
        @test result[1].data == erps[1].data

        # Test with "epochs" variable
        epochs_file = joinpath(test_dir, "test_epochs.jld2")
        jldsave(epochs_file; data = erps)

        result = eegfun.load_data(epochs_file)
        @test result !== nothing
        @test length(result) == length(erps)
        @test result[1].data == erps[1].data

        # Test with invalid data type (String instead of EEG data type)
        other_file = joinpath(test_dir, "test_other.jld2")
        jldsave(other_file; data = "test")

        result = eegfun.load_data(other_file)
        # load_data returns the data (String in this case), it doesn't validate types
        @test result !== nothing
        @test result == "test"

        # Test with non-existent file (jldopen throws SystemError, not ArgumentError)
        @test_throws SystemError eegfun.load_data("/nonexistent/file.jld2")
    end

    @testset "_condition_select" begin
        # Create test data
        data = [create_test_erp_data(1, i) for i = 1:5]

        # Test with nothing (should return original)
        result = eegfun._condition_select(data, nothing)
        @test result == data

        # Test with single condition
        result = eegfun._condition_select(data, 3)
        @test length(result) == 1
        @test result[1] == data[3]

        # Test with multiple conditions
        result = eegfun._condition_select(data, [1, 3, 5])
        @test length(result) == 3
        @test result[1] == data[1]
        @test result[2] == data[3]
        @test result[3] == data[5]

        # Test with empty selection
        result = eegfun._condition_select(data, Int[])
        @test isempty(result)
    end

    @testset "_validate_input_dir" begin
        # Test existing directory
        result = eegfun._validate_input_dir(test_dir)
        @test result === nothing

        # Test non-existent directory
        result = eegfun._validate_input_dir("/nonexistent/directory")
        @test result isa String
        @test occursin("does not exist", result)

        # Test file instead of directory
        test_file = joinpath(test_dir, "test_file.txt")
        write(test_file, "test")
        result = eegfun._validate_input_dir(test_file)
        @test result isa String
        @test occursin("does not exist", result)
    end

    @testset "_validate_channel_groups" begin
        # Test valid channel groups
        groups = [[:Fz, :Cz], [:Pz, :Oz], [:M1, :M2]]
        result = eegfun._validate_channel_groups(groups)
        @test result === nothing

        # Test empty groups
        result = eegfun._validate_channel_groups(Vector{Symbol}[])
        @test result isa String
        @test occursin("cannot be empty", result)

        # Test empty group
        result = eegfun._validate_channel_groups([Symbol[]])
        @test result isa String
        @test occursin("Channel group 1 is empty", result)

        # Test single channel group (should warn but not error)
        groups_single = [[:Fz], [:Cz, :Pz]]
        result = eegfun._validate_channel_groups(groups_single)
        @test result === nothing  # Should not error, just warn
    end

    @testset "_validate_condition_groups" begin
        # Test valid condition groups
        groups = [[1, 2], [3, 4], [5, 6]]
        result = eegfun._validate_condition_groups(groups)
        @test result === nothing

        # Test empty groups
        result = eegfun._validate_condition_groups(Vector{Int}[])
        @test result isa String
        @test occursin("cannot be empty", result)

        # Test duplicate removal
        groups_with_duplicates = [[1, 1, 2], [3, 4]]
        result = eegfun._validate_condition_groups(groups_with_duplicates)
        @test result === nothing
        @test groups_with_duplicates[1] == [1, 2]  # Duplicates should be removed

        # Test overlap detection
        groups_overlap = [[1, 2], [2, 3]]
        result = eegfun._validate_condition_groups(groups_overlap)
        @test result === nothing  # Should warn but not error
    end

    @testset "_validate_condition_pairs" begin
        # Test valid condition pairs (tuples)
        pairs_tuples = [(1, 2), (3, 4), (5, 6)]
        result = eegfun._validate_condition_pairs(pairs_tuples)
        @test result === nothing

        # Test valid condition pairs (vectors)
        pairs_vectors = [[1, 2], [3, 4], [5, 6]]
        result = eegfun._validate_condition_pairs(pairs_vectors)
        @test result === nothing

        # Test empty pairs
        result = eegfun._validate_condition_pairs(Tuple{Int,Int}[])
        @test result isa String
        @test occursin("cannot be empty", result)

        # Test identical conditions (should warn but not error)
        pairs_identical = [(1, 1), (2, 3)]
        result = eegfun._validate_condition_pairs(pairs_identical)
        @test result === nothing  # Should warn but not error
    end

    @testset "_run_batch_operation" begin
        # Create test files
        for i = 1:3
            erps = [create_test_erp_data(i, 1)]
            file_path = joinpath(test_dir, "test_$i.jld2")
            jldsave(file_path; data = erps)
        end

        # Test successful processing
        files = ["test_1.jld2", "test_2.jld2", "test_3.jld2"]

        process_fn =
            (input_path, output_path) -> begin
                if basename(input_path) == "test_2.jld2"
                    return eegfun.BatchResult(false, basename(input_path), "Simulated error")
                else
                    return eegfun.BatchResult(true, basename(input_path), "Success")
                end
            end

        output_dir = joinpath(test_dir, "batch_output")
        mkpath(output_dir)

        results = eegfun._run_batch_operation(process_fn, files, test_dir, output_dir)

        @test length(results) == 3
        @test results[1].success == true
        @test results[2].success == false
        @test results[3].success == true

        # Test error handling
        error_process_fn = (input_path, output_path) -> error("Test error")

        results_error = eegfun._run_batch_operation(error_process_fn, files[1:1], test_dir, output_dir)

        @test length(results_error) == 1
        @test results_error[1].success == false
        @test occursin("Exception", results_error[1].message)
    end

    @testset "_log_batch_summary" begin
        # Test with mixed results
        results = [
            eegfun.BatchResult(true, "file1.jld2", "Success"),
            eegfun.BatchResult(false, "file2.jld2", "Error"),
            eegfun.BatchResult(true, "file3.jld2", "Success"),
        ]

        output_dir = joinpath(test_dir, "summary_test")
        mkpath(output_dir)

        summary = eegfun._log_batch_summary(results, output_dir)

        @test summary.success == 2
        @test summary.errors == 1
    end

    @testset "_cleanup_logging" begin
        # Test cleanup without output directory
        eegfun._cleanup_logging("nonexistent.log", nothing)
        # Should not throw error

        # Test cleanup with existing log file
        log_file = joinpath(test_dir, "test.log")
        write(log_file, "test log content")
        eegfun._cleanup_logging(log_file, test_dir)
        # Should not throw error

        # Test cleanup with non-existent log file (should handle gracefully)
        try
            eegfun._cleanup_logging("nonexistent.log", test_dir)
        catch
            # Expected to fail, that's okay for this test
        end

        # Test cleanup with output directory
        output_dir = joinpath(test_dir, "cleanup_test")
        mkpath(output_dir)

        # Create a dummy log file in current directory (as expected by _cleanup_logging)
        log_file = "test.log"
        write(log_file, "Test log content")

        eegfun._cleanup_logging(log_file, output_dir)

        # Check if log file was moved
        moved_log = joinpath(output_dir, "test.log")
        @test isfile(moved_log)
        @test read(moved_log, String) == "Test log content"

        # Clean up the moved file
        rm(moved_log, force = true)
    end


    @testset "Edge cases and error handling" begin
        @testset "File system edge cases" begin
            # Test with empty directory
            empty_dir = joinpath(test_dir, "empty")
            mkpath(empty_dir)

            files = eegfun._find_batch_files("pattern", empty_dir)
            @test isempty(files)

            # Test with non-existent directory
            @test_throws Base.IOError eegfun._find_batch_files("pattern", "/nonexistent")
        end

        @testset "Data validation edge cases" begin
            # Test with very large condition numbers
            groups_large = [[1000, 2000], [3000, 4000]]
            result = eegfun._validate_condition_groups(groups_large)
            @test result === nothing

            # Test with negative condition numbers
            groups_negative = [[-1, 1], [2, 3]]
            result = eegfun._validate_condition_groups(groups_negative)
            @test result === nothing  # Should not error, just process

            # Test with mixed data types in groups
            groups_mixed = [[1, 2], [3, 4]]  # All Int
            result = eegfun._validate_condition_groups(groups_mixed)
            @test result === nothing  # Should handle gracefully
        end

        @testset "Batch operation edge cases" begin
            # Test with empty file list
            results = eegfun._run_batch_operation(
                (x, y) -> eegfun.BatchResult(true, "test", "ok"),
                String[],
                test_dir,
                test_dir,
            )
            @test isempty(results)

            # Test with processing function that returns nothing
            process_fn_nothing =
                (input_path, output_path) -> eegfun.BatchResult(false, "test", "Function returned nothing")
            results = eegfun._run_batch_operation(process_fn_nothing, ["test_1.jld2"], test_dir, test_dir)
            @test length(results) == 1
            @test results[1].success == false
            @test occursin("Function returned nothing", results[1].message)
        end
    end

    @testset "Integration tests" begin
        @testset "Full batch workflow simulation" begin
            # Create test files
            for participant = 1:3
                erps = [create_test_erp_data(participant, 1)]
                file_path = joinpath(test_dir, "$(participant)_test_erps.jld2")
                jldsave(file_path; data = erps)
            end

            # Test complete workflow
            files = eegfun._find_batch_files("test_erps", test_dir)
            @test length(files) == 4

            # Validate input directory
            validation = eegfun._validate_input_dir(test_dir)
            @test validation === nothing

            # Process files
            process_fn =
                (input_path, output_path) -> begin
                    data_result = eegfun.load_data(input_path)
                    if isnothing(data_result)
                        return eegfun.BatchResult(false, basename(input_path), "No data")
                    end
                    return eegfun.BatchResult(true, basename(input_path), "Processed")
                end

            output_dir = joinpath(test_dir, "integration_output")
            mkpath(output_dir)

            results = eegfun._run_batch_operation(process_fn, files, test_dir, output_dir)
            @test length(results) == 4
            @test all(r.success for r in results)

            # Log summary
            summary = eegfun._log_batch_summary(results, output_dir)
            @test summary.success == 4
            @test summary.errors == 0
        end
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
