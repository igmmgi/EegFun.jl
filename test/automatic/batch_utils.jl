"""
Test suite for src/analysis/batch_utils.jl
"""

using Test
using JLD2
using DataFrames
using Logging

@testset "Batch Utils" begin
    # Create temporary test directory
    test_dir = mktempdir()

    # Save initial global logger state to restore later
    initial_logger = global_logger()

    # Track log files created in current directory for cleanup
    created_log_files = String[]

    @testset "BatchResult struct" begin
        # Test BatchResult creation and access
        result = eegfun.BatchResult(true, "test_file.jld2", "Success message")

        @test result.success == true
        @test result.filename == "test_file.jld2"
        @test result.message == "Success message"

        # Test immutable nature
        @test_throws ErrorException result.success = false

        # Test with false success
        result_fail = eegfun.BatchResult(false, "error.jld2", "Error message")
        @test result_fail.success == false
        @test result_fail.filename == "error.jld2"
        @test result_fail.message == "Error message"
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

        # Test with single Int
        config_single = eegfun.BatchConfig("test", "/input", "/output", 1, 1)
        @test config_single.participants == 1
        @test config_single.conditions == 1
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
        files_filtered = eegfun._find_batch_files("erps_cleaned", test_dir, eegfun.participants([1, 3, 5]))
        @test length(files_filtered) == 3
        @test "1_erps_cleaned.jld2" in files_filtered
        @test "3_erps_cleaned.jld2" in files_filtered
        @test "5_erps_cleaned.jld2" in files_filtered

        # Test single participant
        files_single = eegfun._find_batch_files("erps_cleaned", test_dir, eegfun.participants(2))
        @test length(files_single) == 1
        @test "2_erps_cleaned.jld2" in files_single

        # Test no participants (should return all)
        files_all = eegfun._find_batch_files("erps_cleaned", test_dir, nothing)
        @test length(files_all) == 5

        # Test non-matching pattern
        files_none = eegfun._find_batch_files("nonexistent", test_dir)
        @test isempty(files_none)

        # Test with empty directory
        empty_dir = joinpath(test_dir, "empty")
        mkpath(empty_dir)
        files_empty = eegfun._find_batch_files("pattern", empty_dir)
        @test isempty(files_empty)
    end

    @testset "_condition_select" begin
        # Create test data
        data = [create_test_erp_data(1, i) for i = 1:5]

        # Test with nothing (should return original)
        result = eegfun._condition_select(data, nothing)
        @test result == data

        # Test with single condition (Int)
        result = eegfun._condition_select(data, 3)
        @test length(result) == 1
        @test result[1] == data[3]

        # Test with multiple conditions (Vector{Int})
        result = eegfun._condition_select(data, [1, 3, 5])
        @test length(result) == 3
        @test result[1] == data[1]
        @test result[2] == data[3]
        @test result[3] == data[5]

        # Test with empty selection
        result = eegfun._condition_select(data, Int[])
        @test isempty(result)

        # Test with empty data
        result = eegfun._condition_select([], [1, 2])
        @test isempty(result)

        # Test with empty data and nothing
        result = eegfun._condition_select([], nothing)
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

        # Test multiple single-channel groups
        groups_multi_single = [[:Fz], [:Cz], [:Pz]]
        result = eegfun._validate_channel_groups(groups_multi_single)
        @test result === nothing  # Should warn but not error

        # Test large groups
        groups_large = [[:Fz, :Cz, :Pz, :Oz, :M1, :M2]]
        result = eegfun._validate_channel_groups(groups_large)
        @test result === nothing
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
        groups_copy = deepcopy(groups_with_duplicates)
        result = eegfun._validate_condition_groups(groups_copy)
        @test result === nothing
        @test groups_copy[1] == [1, 2]  # Duplicates should be removed

        # Test multiple duplicates
        groups_multi_dup = [[1, 1, 1, 2, 2], [3, 4]]
        groups_copy2 = deepcopy(groups_multi_dup)
        result = eegfun._validate_condition_groups(groups_copy2)
        @test result === nothing
        @test groups_copy2[1] == [1, 2]  # All duplicates removed

        # Test overlap detection (should warn but not error)
        groups_overlap = [[1, 2], [2, 3]]
        groups_copy3 = deepcopy(groups_overlap)
        result = eegfun._validate_condition_groups(groups_copy3)
        @test result === nothing  # Should warn but not error

        # Test multiple overlaps
        groups_multi_overlap = [[1, 2], [2, 3], [3, 4]]
        groups_copy4 = deepcopy(groups_multi_overlap)
        result = eegfun._validate_condition_groups(groups_copy4)
        @test result === nothing

        # Test single condition groups
        groups_single = [[1], [2], [3]]
        result = eegfun._validate_condition_groups(groups_single)
        @test result === nothing
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

        result = eegfun._validate_condition_pairs(Vector{Int}[])
        @test result isa String
        @test occursin("cannot be empty", result)

        # Test identical conditions (should warn but not error)
        pairs_identical = [(1, 1), (2, 3)]
        result = eegfun._validate_condition_pairs(pairs_identical)
        @test result === nothing  # Should warn but not error

        # Test multiple identical pairs
        pairs_multi_identical = [(1, 1), (2, 2), (3, 3)]
        result = eegfun._validate_condition_pairs(pairs_multi_identical)
        @test result === nothing

        # Test mixed tuples and vectors (should error at type level, but test what we can)
        # Note: This would be a type error, so we can't easily test it
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
        @test results[1].filename == "test_1.jld2"
        @test results[2].filename == "test_2.jld2"
        @test results[3].filename == "test_3.jld2"

        # Test error handling
        error_process_fn = (input_path, output_path) -> error("Test error")

        results_error = eegfun._run_batch_operation(error_process_fn, files[1:1], test_dir, output_dir)

        @test length(results_error) == 1
        @test results_error[1].success == false
        @test occursin("Exception", results_error[1].message)
        @test occursin("Test error", results_error[1].message)

        # Test with empty file list
        results_empty = eegfun._run_batch_operation(process_fn, String[], test_dir, output_dir)
        @test isempty(results_empty)

        # Test with custom operation name
        results_custom = eegfun._run_batch_operation(
            process_fn,
            files[1:1],
            test_dir,
            output_dir;
            operation_name = "Custom Operation",
        )
        @test length(results_custom) == 1
    end

    @testset "_log_batch_summary" begin
        # Test with mixed results
        results = Vector{eegfun.BatchResult}([
            eegfun.BatchResult(true, "file1.jld2", "Success"),
            eegfun.BatchResult(false, "file2.jld2", "Error"),
            eegfun.BatchResult(true, "file3.jld2", "Success"),
        ])

        output_dir = joinpath(test_dir, "summary_test")
        mkpath(output_dir)

        summary = eegfun._log_batch_summary(results, output_dir)

        @test summary.success == 2
        @test summary.errors == 1

        # Test with all successes
        results_all_success = Vector{eegfun.BatchResult}([
            eegfun.BatchResult(true, "file1.jld2", "Success"),
            eegfun.BatchResult(true, "file2.jld2", "Success"),
        ])
        summary_all = eegfun._log_batch_summary(results_all_success, output_dir)
        @test summary_all.success == 2
        @test summary_all.errors == 0

        # Test with all errors
        results_all_error = Vector{eegfun.BatchResult}([
            eegfun.BatchResult(false, "file1.jld2", "Error 1"),
            eegfun.BatchResult(false, "file2.jld2", "Error 2"),
        ])
        summary_all_error = eegfun._log_batch_summary(results_all_error, output_dir)
        @test summary_all_error.success == 0
        @test summary_all_error.errors == 2

        # Test with empty results
        summary_empty = eegfun._log_batch_summary(Vector{eegfun.BatchResult}(), output_dir)
        @test summary_empty.success == 0
        @test summary_empty.errors == 0
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

        # Test cleanup with output directory
        output_dir = joinpath(test_dir, "cleanup_test")
        mkpath(output_dir)

        # Create a dummy log file in current directory (as expected by _cleanup_logging)
        log_file = "test_cleanup.log"
        write(log_file, "Test log content")
        push!(created_log_files, log_file)

        eegfun._cleanup_logging(log_file, output_dir)

        # Check if log file was moved
        moved_log = joinpath(output_dir, "test_cleanup.log")
        @test isfile(moved_log)
        @test read(moved_log, String) == "Test log content"

        # Clean up the moved file
        rm(moved_log, force = true)

        # Remove from tracking since it was moved
        filter!(f -> f != log_file, created_log_files)

        # Test cleanup when log file and destination are the same
        log_file_same = joinpath(output_dir, "same.log")
        write(log_file_same, "Same file")
        eegfun._cleanup_logging(log_file_same, output_dir)
        @test isfile(log_file_same)  # Should still exist (not moved to itself)
    end


    @testset "@log_call macro" begin
        # Test @log_call with function name only (current simplified form)
        function test_func1(x, y; opt1 = 1, opt2 = "test")
            eegfun.@log_call "test_func1"
            return x + y + opt1
        end

        # Call the function (should log)
        result = test_func1(1, 2; opt1 = 3, opt2 = "test")
        @test result == 6

        # Test @log_call with different function name
        function test_func2(x, y, z; opt1 = 1)
            eegfun.@log_call "test_func2"
            return x + y + z + opt1
        end

        result = test_func2(1, 2, 3; opt1 = 4)
        @test result == 10

        # Test @log_call in function with single arg
        function test_func3(x; opt1 = 1)
            eegfun.@log_call "test_func3"
            return x + opt1
        end

        result = test_func3(5; opt1 = 2)
        @test result == 7

        # Test @log_call in function with no positional args
        function test_func4(; opt1 = 1, opt2 = "test")
            eegfun.@log_call "test_func4"
            return opt1
        end

        result = test_func4(opt1 = 3, opt2 = "test")
        @test result == 3

    end

    @testset "Edge cases and error handling" begin
        @testset "File system edge cases" begin
            # Test with non-existent directory
            @test_throws Base.IOError eegfun._find_batch_files("pattern", "/nonexistent")

            # Test with file instead of directory
            test_file = joinpath(test_dir, "test_file.txt")
            write(test_file, "test")
            @test_throws Base.IOError eegfun._find_batch_files("pattern", test_file)
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

            # Test with zero
            groups_zero = [[0, 1], [2, 3]]
            result = eegfun._validate_condition_groups(groups_zero)
            @test result === nothing
        end

        @testset "Batch operation edge cases" begin
            # Test with empty file list (use NullLogger to avoid stream issues)
            results = with_logger(NullLogger()) do
                eegfun._run_batch_operation((x, y) -> eegfun.BatchResult(true, "test", "ok"), String[], test_dir, test_dir)
            end
            @test isempty(results)

            # Test with very long filenames
            long_filename = "a" ^ 200 * ".jld2"
            erps = [create_test_erp_data(1, 1)]
            long_file_path = joinpath(test_dir, long_filename)
            jldsave(long_file_path; data = erps)

            files_long = [long_filename]
            results_long = with_logger(NullLogger()) do
                eegfun._run_batch_operation(
                    (x, y) -> eegfun.BatchResult(true, basename(x), "Success"),
                    files_long,
                    test_dir,
                    test_dir,
                )
            end
            @test length(results_long) == 1
            @test results_long[1].success == true
        end

        @testset "Condition selection edge cases" begin
            data = [create_test_erp_data(1, i) for i = 1:5]

            # Test with out-of-bounds indices (throws BoundsError for indices > length)
            @test_throws BoundsError eegfun._condition_select(data, [10])

            # Test with negative indices (throws BoundsError from array indexing)
            @test_throws BoundsError eegfun._condition_select(data, [-1])

            # Test with zero index (throws BoundsError from array indexing)
            @test_throws BoundsError eegfun._condition_select(data, [0])
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
            @test length(files) == 3

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

            # Use NullLogger to avoid stream initialization issues during tests
            results = with_logger(NullLogger()) do
                eegfun._run_batch_operation(process_fn, files, test_dir, output_dir)
            end
            @test length(results) == 3
            @test all(r.success for r in results)

            # Log summary (also with NullLogger)
            summary = with_logger(NullLogger()) do
                eegfun._log_batch_summary(results, output_dir)
            end
            @test summary.success == 3
            @test summary.errors == 0
        end
    end

    # Cleanup: Restore global logger state and remove any leftover log files
    # Helper to safely execute cleanup operations
    safe_cleanup(f) =
        try
            ;
            f();
        catch
            ;
        end

    # Close any open global logging and restore initial logger
    safe_cleanup(() -> eegfun.close_global_logging())
    global_logger(initial_logger)

    # Clean up log files created in current directory
    for log_file in created_log_files
        safe_cleanup(() -> isfile(log_file) && rm(log_file, force = true))
    end

    # Remove test directory
    safe_cleanup(() -> rm(test_dir, recursive = true))
end
