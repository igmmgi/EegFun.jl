using Test
using DataFrames
using EegFun
using JLD2
using Statistics

@testset "Batch Combine Conditions" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2, 3]
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:4)  # 4 conditions, 10 epochs each
                filename = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
                jldsave(filename; data = epochs)
                @test isfile(filename)
            end
        end

        @testset "Basic condition combining" begin
            output_dir = joinpath(test_dir, "combined_basic")

            # Combine conditions 1,2 into group 1 and 3,4 into group 2
            result = EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            @test result !== nothing
            @test result.success == 3
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that combined files exist
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "3_epochs_cleaned.jld2"))

            # Load and verify combined data
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 2  # 2 groups

            # Check epoch counts: each group should have 6 epochs (3 from each original condition)
            @test length(combined_epochs[1].data) == 20  # Group 1: conditions 1+2
            @test length(combined_epochs[2].data) == 20  # Group 2: conditions 3+4

            # Verify combined condition labels
            # Group 1 is a combination of original conditions 1 and 2
            @test combined_epochs[1].condition == 1  # group_idx becomes condition number
            # Group 2 is a combination of original conditions 3 and 4
            @test combined_epochs[2].condition == 2  # group_idx becomes condition number
            # Note: condition column is no longer in DataFrames, it's in struct fields
        end

        @testset "Combine specific participants" begin
            output_dir = joinpath(test_dir, "combined_participant")

            result = EegFun.condition_combine(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(2),
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "3_epochs_cleaned.jld2"))
        end

        @testset "Combine multiple participants" begin
            output_dir = joinpath(test_dir, "combined_multi_participants")

            result = EegFun.condition_combine(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 3]),
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "3_epochs_cleaned.jld2"))
        end

        @testset "Single condition groups" begin
            output_dir = joinpath(test_dir, "combined_single")

            # Each condition becomes its own group (no actual combining)
            result = EegFun.condition_combine("epochs_cleaned", [[1], [2], [3], [4]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            # Load and verify
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 4  # 4 groups

            # Each group should have 3 epochs (original count)
            for i = 1:4
                @test length(combined_epochs[i].data) == 10
                @test combined_epochs[i].condition == i
            end
        end

        @testset "Overlapping condition groups" begin
            output_dir = joinpath(test_dir, "combined_overlapping")

            # Overlapping groups: [1,2], [2,3], [3,4]
            result = EegFun.condition_combine("epochs_cleaned", [[1, 2], [2, 3], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 3  # 3 groups

            # Check epoch counts
            @test length(combined_epochs[1].data) == 20  # Group 1: conditions 1+2
            @test length(combined_epochs[2].data) == 20  # Group 2: conditions 2+3
            @test length(combined_epochs[3].data) == 20  # Group 3: conditions 3+4
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception EegFun.condition_combine("epochs_cleaned", [[1, 2]], input_dir = "/nonexistent/path")

            # Invalid pattern (doesn't contain 'epochs')
            @test_throws Exception EegFun.condition_combine("erps_cleaned", [[1, 2]], input_dir = test_dir)

            # Invalid condition groups
            @test_throws Exception EegFun.condition_combine(
                "epochs_cleaned",
                [],  # Empty groups
                input_dir = test_dir,
            )

            @test_throws Exception EegFun.condition_combine(
                "epochs_cleaned",
                [[1, 2], []],  # Empty group
                input_dir = test_dir,
            )
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "combined_invalid_condition")

            # Request condition 5 when only 4 exist
            result = EegFun.condition_combine(
                "epochs_cleaned",
                [[1, 5]],  # Condition 5 doesn't exist
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            result = EegFun.condition_combine("epochs_cleaned", [[1, 2]], input_dir = empty_dir)

            @test result === nothing  # No files to process
        end

        @testset "Data integrity - concatenation preserves all epochs" begin
            output_dir = joinpath(test_dir, "combined_integrity")

            # Get original epoch counts
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "data")
            original_counts = [length(cond.data) for cond in original_epochs]

            # Combine conditions 1 and 2
            EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            # Load combined data
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")

            # Verify total epoch count is preserved
            @test length(combined_epochs[1].data) == original_counts[1] + original_counts[2]
            @test length(combined_epochs[2].data) == original_counts[3] + original_counts[4]

            # Verify all original epochs are present (check by total count)
            # Each group should have 6 epochs total (3 from each original condition)
            @test length(combined_epochs[1].data) == 20  # 3 from condition 1 + 3 from condition 2
        end

        @testset "Metadata preservation" begin
            output_dir = joinpath(test_dir, "combined_metadata")

            EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")

            # Verify metadata columns exist
            @test hasproperty(combined_epochs[1].data[1], :time)
            # condition and condition_name are now in struct, not DataFrame
            @test hasproperty(combined_epochs[1], :condition)
            @test hasproperty(combined_epochs[1], :condition_name)
            @test hasproperty(combined_epochs[1].data[1], :epoch)
            # Verify channel data is preserved
            @test hasproperty(combined_epochs[1].data[1], :Ch1)
            @test hasproperty(combined_epochs[1].data[1], :Ch2)
        end

        @testset "Layout and sample rate preservation" begin
            output_dir = joinpath(test_dir, "combined_layout")

            # Get original metadata
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "data")
            original_layout = original_epochs[1].layout
            original_fs = original_epochs[1].sample_rate

            EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")

            # Verify layout and sample rate are preserved
            @test combined_epochs[1].layout.data == original_layout.data
            @test combined_epochs[1].sample_rate == original_fs
            @test combined_epochs[1] isa EegFun.EpochData
        end

        @testset "Different epoch counts per condition" begin
            # TODO: adapr create_test_epoch_data to have different number of epochs per condition
            # Create data with different number of epochs per condition
            dat = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 2)

            # Save and process
            var_dir = joinpath(test_dir, "var_epochs")
            mkpath(var_dir)
            jldsave(joinpath(var_dir, "1_epochs_var.jld2"); data = dat)

            output_dir = joinpath(test_dir, "combined_var")
            result = EegFun.condition_combine("epochs_var", [[1, 2]], input_dir = var_dir, output_dir = output_dir)

            @test result.success == 1

            # Load and verify epoch counts
            combined_epochs = load(joinpath(output_dir, "1_epochs_var.jld2"), "data")
            @test length(combined_epochs) == 1
            @test length(combined_epochs[1].data) == 4  # 2 + 5 epochs
        end

        @testset "Empty condition groups" begin
            output_dir = joinpath(test_dir, "combined_empty_groups")

            # Test with empty condition groups (should fail validation)
            @test_throws Exception EegFun.condition_combine("epochs_cleaned", [], input_dir = test_dir, output_dir = output_dir)
        end

        @testset "Single condition per group" begin
            output_dir = joinpath(test_dir, "combined_single_per_group")

            # Each condition in its own group (no actual combining)
            result = EegFun.condition_combine("epochs_cleaned", [[1], [2], [3], [4]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 4

            # Each group should have the same number of epochs as original
            for i = 1:4
                @test length(combined_epochs[i].data) == 10  # Original epoch count
            end
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "combined_with_log")

            result = EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            # Check log file exists
            log_file = joinpath(output_dir, "condition_combine.log")
            @test isfile(log_file)

            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "condition_combine")
            @test contains(log_contents, "epochs_cleaned")
            @test contains(log_contents, "Found")
        end

        @testset "Output directory naming" begin
            # Test default output directory naming
            result = EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir)

            # Should create directory with pattern and groups in name
            expected_dir = joinpath(test_dir, "combined_epochs_cleaned_1-2_3-4")
            @test isdir(expected_dir)
        end

        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "combined_return_check")

            result = EegFun.condition_combine("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            # Check result structure
            @test hasfield(typeof(result), :success)
            @test hasfield(typeof(result), :errors)
            @test result.success isa Integer
            @test result.errors isa Integer
            @test result.success >= 0
            @test result.errors >= 0
            @test result.success + result.errors == 3  # Total files processed
        end

        @testset "Duplicate conditions within groups" begin
            output_dir = joinpath(test_dir, "combined_duplicates")

            # Test with duplicate conditions within groups: [[1, 1, 2], [3, 3, 4]]
            result = EegFun.condition_combine("epochs_cleaned", [[1, 1, 2], [3, 3, 4]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            # Load and verify - duplicates should be removed
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 2

            # Each group should have 6 epochs (3 from each unique condition)
            @test length(combined_epochs[1].data) == 20  # Group 1: conditions 1+2 (duplicates removed)
            @test length(combined_epochs[2].data) == 20  # Group 2: conditions 3+4 (duplicates removed)
        end

        @testset "Empty groups after duplicate removal" begin
            output_dir = joinpath(test_dir, "combined_empty_after_duplicates")

            # Test with group that becomes single condition after removing duplicates: [[1, 1, 1]]
            result = EegFun.condition_combine("epochs_cleaned", [[1, 1, 1], [2, 3]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            # Should create 2 groups (1 becomes [1] after removing duplicates)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 2

            # First group should have 3 epochs (from condition 1), second group should have 6 epochs (from conditions 2+3)
            @test length(combined_epochs[1].data) == 10  # Group 1: condition 1 only
            @test length(combined_epochs[2].data) == 20  # Group 2: conditions 2+3
        end

        @testset "Negative condition numbers" begin
            output_dir = joinpath(test_dir, "combined_negative")

            # Test with negative condition numbers (should fail for all files)
            result = EegFun.condition_combine("epochs_cleaned", [[1, -1], [2, 3]], input_dir = test_dir, output_dir = output_dir)

            # Should fail for all files because negative conditions don't exist
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Zero condition numbers" begin
            output_dir = joinpath(test_dir, "combined_zero")

            # Test with zero condition numbers (should fail for all files)
            result = EegFun.condition_combine("epochs_cleaned", [[0, 1], [2, 3]], input_dir = test_dir, output_dir = output_dir)

            # Should fail for all files because zero conditions don't exist
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Very large condition numbers" begin
            output_dir = joinpath(test_dir, "combined_large")

            # Test with very large condition numbers that don't exist
            result = EegFun.condition_combine("epochs_cleaned", [[1000, 1001]], input_dir = test_dir, output_dir = output_dir)

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Single condition in multiple groups" begin
            output_dir = joinpath(test_dir, "combined_single_multiple")

            # Test with same condition in multiple groups: [[1], [1], [2]]
            result = EegFun.condition_combine("epochs_cleaned", [[1], [1], [2]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            # Should create 3 groups (even though 2 are identical)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 3

            # Each group should have 3 epochs
            for i = 1:3
                @test length(combined_epochs[i].data) == 10
            end
        end

        @testset "All conditions in single group" begin
            output_dir = joinpath(test_dir, "combined_all_single")

            # Test with all conditions in one group: [[1, 2, 3, 4]]
            result = EegFun.condition_combine("epochs_cleaned", [[1, 2, 3, 4]], input_dir = test_dir, output_dir = output_dir)

            @test result.success == 3

            # Should create 1 group with all epochs
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 1
            @test length(combined_epochs[1].data) == 40  # 3 epochs from each of 4 conditions
        end

        @testset "Non-integer condition groups" begin
            output_dir = joinpath(test_dir, "combined_non_integer")

            # Test with non-integer condition groups (should fail at type level)
            @test_throws MethodError EegFun.condition_combine(
                "epochs_cleaned",
                [["1", "2"], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )
        end

        @testset "Nested empty groups" begin
            output_dir = joinpath(test_dir, "combined_nested_empty")

            # Test with nested empty groups: [[], [1, 2]]
            @test_throws Exception EegFun.condition_combine("epochs_cleaned", [[], [1, 2]], input_dir = test_dir, output_dir = output_dir)
        end

        @testset "Very large number of groups" begin
            output_dir = joinpath(test_dir, "combined_many_groups")

            # Test with many groups: [[1], [2], [3], [4], [1], [2], [3], [4]]
            result = EegFun.condition_combine(
                "epochs_cleaned",
                [[1], [2], [3], [4], [1], [2], [3], [4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Should create 8 groups (including duplicates)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "data")
            @test length(combined_epochs) == 8
        end

        @testset "File corruption handling" begin
            corrupt_dir = joinpath(test_dir, "corrupt_files")
            mkpath(corrupt_dir)

            # Create a corrupted JLD2 file
            corrupt_file = joinpath(corrupt_dir, "1_epochs_corrupt.jld2")
            write(corrupt_file, "This is not a valid JLD2 file")

            output_dir = joinpath(test_dir, "combined_corrupt")
            result = EegFun.condition_combine("epochs_corrupt", [[1, 2]], input_dir = corrupt_dir, output_dir = output_dir)

            # Should fail for the corrupted file
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Missing epochs variable" begin
            missing_var_dir = joinpath(test_dir, "missing_var")
            mkpath(missing_var_dir)

            # Create file with invalid data type (String instead of Vector{EpochData}) to test error handling
            jldsave(joinpath(missing_var_dir, "1_epochs_missing.jld2"); data = "invalid_data")

            output_dir = joinpath(test_dir, "combined_missing_var")
            result = EegFun.condition_combine("epochs_missing", [[1, 2]], input_dir = missing_var_dir, output_dir = output_dir)

            # Should fail for the file with invalid data type
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Empty epochs data" begin
            empty_epochs_dir = joinpath(test_dir, "empty_epochs")
            mkpath(empty_epochs_dir)

            # Create file with empty epochs array
            jldsave(joinpath(empty_epochs_dir, "1_epochs_empty.jld2"); data = EegFun.EpochData[])

            output_dir = joinpath(test_dir, "combined_empty_epochs")
            result = EegFun.condition_combine("epochs_empty", [[1, 2]], input_dir = empty_epochs_dir, output_dir = output_dir)

            # Should fail because no conditions to combine
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Many channels" begin
            # Test with more channels (5)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)

            dat = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 3, n_channels = 1000)

            jldsave(joinpath(many_ch_dir, "1_epochs_many.jld2"); data = dat)

            # Combine
            output_dir = joinpath(test_dir, "combined_many_ch")
            result = EegFun.condition_combine("epochs_many", [[1, 2]], input_dir = many_ch_dir, output_dir = output_dir)

            @test result.success == 1

            combined_epochs = load(joinpath(output_dir, "1_epochs_many.jld2"), "data")

            # Verify all channels are present
            # Verify epoch count
            @test length(combined_epochs[1].data) == 6  # 3 from each condition
        end

    finally
        # Cleanup
        rm(test_dir, recursive = true, force = true)
    end
end
