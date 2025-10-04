using Test
using DataFrames
using eegfun
using JLD2
using Statistics

@testset "Batch Combine Conditions" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Helper to create test EpochData with multiple conditions
        function create_test_epoch_data(n_conditions::Int = 3, n_epochs_per_condition::Int = 4)
            epochs = eegfun.EpochData[]
            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            for cond = 1:n_conditions
                # Create multiple epochs per condition
                dfs = DataFrame[]
                for ep = 1:n_epochs_per_condition
                    # Add some trial-to-trial variability and condition-specific signal
                    ch1 = sin.(2π .* (5 + cond) .* t) .+ 0.1 .* randn(n_samples)
                    ch2 = cos.(2π .* (5 + cond) .* t) .+ 0.1 .* randn(n_samples)

                    df = DataFrame(
                        time = collect(t),
                        sample = 1:n_samples,
                        condition = fill(cond, n_samples),
                        condition_name = fill("condition_$cond", n_samples),
                        epoch = fill(ep, n_samples),
                        Fz = ch1,
                        Cz = ch2,
                    )
                    push!(dfs, df)
                end

                layout =
                    eegfun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

                push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
            end

            return epochs
        end

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2, 3]
                epochs = create_test_epoch_data(4, 3)  # 4 conditions, 3 epochs each
                filename = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
                save(filename, "epochs", epochs)
                @test isfile(filename)
            end
        end

        @testset "Basic condition combining" begin
            output_dir = joinpath(test_dir, "combined_basic")

            # Combine conditions 1,2 into group 1 and 3,4 into group 2
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result !== nothing
            @test result.success == 3
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that combined files exist
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "3_epochs_cleaned.jld2"))

            # Load and verify combined data
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 2  # 2 groups

            # Check epoch counts: each group should have 6 epochs (3 from each original condition)
            @test length(combined_epochs[1].data) == 6  # Group 1: conditions 1+2
            @test length(combined_epochs[2].data) == 6  # Group 2: conditions 3+4

            # Verify condition labels are preserved from original data
            # Group 1 contains conditions 1 and 2, so should have both labels
            group1_conditions = unique(vcat([df.condition for df in combined_epochs[1].data]...))
            @test 1 in group1_conditions && 2 in group1_conditions

            # Group 2 contains conditions 3 and 4, so should have both labels  
            group2_conditions = unique(vcat([df.condition for df in combined_epochs[2].data]...))
            @test 3 in group2_conditions && 4 in group2_conditions
        end

        @testset "Combine specific participants" begin
            output_dir = joinpath(test_dir, "combined_participant")

            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 2,
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "3_epochs_cleaned.jld2"))
        end

        @testset "Combine multiple participants" begin
            output_dir = joinpath(test_dir, "combined_multi_participants")

            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = [1, 3],
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
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1], [2], [3], [4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Load and verify
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 4  # 4 groups

            # Each group should have 3 epochs (original count)
            for i = 1:4
                @test length(combined_epochs[i].data) == 3
                @test all(combined_epochs[i].data[1].condition .== i)
            end
        end

        @testset "Overlapping condition groups" begin
            output_dir = joinpath(test_dir, "combined_overlapping")

            # Overlapping groups: [1,2], [2,3], [3,4]
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [2, 3], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 3  # 3 groups

            # Check epoch counts
            @test length(combined_epochs[1].data) == 6  # Group 1: conditions 1+2
            @test length(combined_epochs[2].data) == 6  # Group 2: conditions 2+3
            @test length(combined_epochs[3].data) == 6  # Group 3: conditions 3+4
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2]],
                input_dir = "/nonexistent/path",
            )

            # Invalid pattern (doesn't contain 'epochs')
            @test_throws Exception eegfun.combine_conditions("erps_cleaned", [[1, 2]], input_dir = test_dir)

            # Invalid condition groups
            @test_throws Exception eegfun.combine_conditions(
                "epochs_cleaned",
                [],  # Empty groups
                input_dir = test_dir,
            )

            @test_throws Exception eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], []],  # Empty group
                input_dir = test_dir,
            )
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "combined_invalid_condition")

            # Request condition 5 when only 4 exist
            result = eegfun.combine_conditions(
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

            result = eegfun.combine_conditions("epochs_cleaned", [[1, 2]], input_dir = empty_dir)

            @test result === nothing  # No files to process
        end

        @testset "Data integrity - concatenation preserves all epochs" begin
            output_dir = joinpath(test_dir, "combined_integrity")

            # Get original epoch counts
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            original_counts = [length(cond.data) for cond in original_epochs]

            # Combine conditions 1 and 2
            eegfun.combine_conditions("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            # Load combined data
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")

            # Verify total epoch count is preserved
            @test length(combined_epochs[1].data) == original_counts[1] + original_counts[2]
            @test length(combined_epochs[2].data) == original_counts[3] + original_counts[4]

            # Verify all original epochs are present (check by total count)
            # Each group should have 6 epochs total (3 from each original condition)
            @test length(combined_epochs[1].data) == 6  # 3 from condition 1 + 3 from condition 2
        end

        @testset "Metadata preservation" begin
            output_dir = joinpath(test_dir, "combined_metadata")

            eegfun.combine_conditions("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")

            # Verify metadata columns exist
            @test hasproperty(combined_epochs[1].data[1], :time)
            @test hasproperty(combined_epochs[1].data[1], :condition)
            @test hasproperty(combined_epochs[1].data[1], :condition_name)
            @test hasproperty(combined_epochs[1].data[1], :epoch)
            @test hasproperty(combined_epochs[1].data[1], :sample)

            # Verify channel data is preserved
            @test hasproperty(combined_epochs[1].data[1], :Fz)
            @test hasproperty(combined_epochs[1].data[1], :Cz)
        end

        @testset "Layout and sample rate preservation" begin
            output_dir = joinpath(test_dir, "combined_layout")

            # Get original metadata
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            original_layout = original_epochs[1].layout
            original_fs = original_epochs[1].sample_rate

            eegfun.combine_conditions("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir, output_dir = output_dir)

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")

            # Verify layout and sample rate are preserved
            @test combined_epochs[1].layout.data == original_layout.data
            @test combined_epochs[1].sample_rate == original_fs
            @test combined_epochs[1] isa eegfun.EpochData
        end

        @testset "Different epoch counts per condition" begin
            # Create data with different number of epochs per condition
            epochs_var = eegfun.EpochData[]
            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Condition 1: 2 epochs
            dfs1 = DataFrame[]
            for ep = 1:2
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples),
                )
                push!(dfs1, df)
            end

            # Condition 2: 5 epochs
            dfs2 = DataFrame[]
            for ep = 1:5
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(2, n_samples),
                    condition_name = fill("condition_2", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples),
                )
                push!(dfs2, df)
            end

            layout = eegfun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)

            push!(epochs_var, eegfun.EpochData(dfs1, layout, fs, eegfun.AnalysisInfo()))
            push!(epochs_var, eegfun.EpochData(dfs2, layout, fs, eegfun.AnalysisInfo()))

            # Save and process
            var_dir = joinpath(test_dir, "var_epochs")
            mkpath(var_dir)
            save(joinpath(var_dir, "1_epochs_var.jld2"), "epochs", epochs_var)

            output_dir = joinpath(test_dir, "combined_var")
            result = eegfun.combine_conditions("epochs_var", [[1, 2]], input_dir = var_dir, output_dir = output_dir)

            @test result.success == 1

            # Load and verify epoch counts
            combined_epochs = load(joinpath(output_dir, "1_epochs_var.jld2"), "epochs")
            @test length(combined_epochs) == 1
            @test length(combined_epochs[1].data) == 7  # 2 + 5 epochs
        end

        @testset "Empty condition groups" begin
            output_dir = joinpath(test_dir, "combined_empty_groups")

            # Test with empty condition groups (should fail validation)
            @test_throws Exception eegfun.combine_conditions(
                "epochs_cleaned",
                [],
                input_dir = test_dir,
                output_dir = output_dir,
            )
        end

        @testset "Single condition per group" begin
            output_dir = joinpath(test_dir, "combined_single_per_group")

            # Each condition in its own group (no actual combining)
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1], [2], [3], [4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 4

            # Each group should have the same number of epochs as original
            for i = 1:4
                @test length(combined_epochs[i].data) == 3  # Original epoch count
            end
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "combined_with_log")

            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Check log file exists
            log_file = joinpath(output_dir, "combine_conditions.log")
            @test isfile(log_file)

            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch condition combining started")
            @test contains(log_contents, "combine_conditions")
            @test contains(log_contents, "epochs_cleaned")
        end

        @testset "Output directory naming" begin
            # Test default output directory naming
            result = eegfun.combine_conditions("epochs_cleaned", [[1, 2], [3, 4]], input_dir = test_dir)

            # Should create directory with pattern and groups in name
            expected_dir = joinpath(test_dir, "combined_epochs_cleaned_1-2_3-4")
            @test isdir(expected_dir)
        end

        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "combined_return_check")

            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

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
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 1, 2], [3, 3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Load and verify - duplicates should be removed
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 2

            # Each group should have 6 epochs (3 from each unique condition)
            @test length(combined_epochs[1].data) == 6  # Group 1: conditions 1+2 (duplicates removed)
            @test length(combined_epochs[2].data) == 6  # Group 2: conditions 3+4 (duplicates removed)
        end

        @testset "Empty groups after duplicate removal" begin
            output_dir = joinpath(test_dir, "combined_empty_after_duplicates")

            # Test with group that becomes single condition after removing duplicates: [[1, 1, 1]]
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 1, 1], [2, 3]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Should create 2 groups (1 becomes [1] after removing duplicates)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 2

            # First group should have 3 epochs (from condition 1), second group should have 6 epochs (from conditions 2+3)
            @test length(combined_epochs[1].data) == 3  # Group 1: condition 1 only
            @test length(combined_epochs[2].data) == 6  # Group 2: conditions 2+3
        end

        @testset "Negative condition numbers" begin
            output_dir = joinpath(test_dir, "combined_negative")

            # Test with negative condition numbers (should fail for all files)
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, -1], [2, 3]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should fail for all files because negative conditions don't exist
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Zero condition numbers" begin
            output_dir = joinpath(test_dir, "combined_zero")

            # Test with zero condition numbers (should fail for all files)
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[0, 1], [2, 3]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should fail for all files because zero conditions don't exist
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Very large condition numbers" begin
            output_dir = joinpath(test_dir, "combined_large")

            # Test with very large condition numbers that don't exist
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1000, 1001]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 3
        end

        @testset "Single condition in multiple groups" begin
            output_dir = joinpath(test_dir, "combined_single_multiple")

            # Test with same condition in multiple groups: [[1], [1], [2]]
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1], [1], [2]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Should create 3 groups (even though 2 are identical)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 3

            # Each group should have 3 epochs
            for i = 1:3
                @test length(combined_epochs[i].data) == 3
            end
        end

        @testset "All conditions in single group" begin
            output_dir = joinpath(test_dir, "combined_all_single")

            # Test with all conditions in one group: [[1, 2, 3, 4]]
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1, 2, 3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Should create 1 group with all epochs
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 1
            @test length(combined_epochs[1].data) == 12  # 3 epochs from each of 4 conditions
        end

        @testset "Non-integer condition groups" begin
            output_dir = joinpath(test_dir, "combined_non_integer")

            # Test with non-integer condition groups (should fail at type level)
            @test_throws MethodError eegfun.combine_conditions(
                "epochs_cleaned",
                [["1", "2"], [3, 4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )
        end

        @testset "Nested empty groups" begin
            output_dir = joinpath(test_dir, "combined_nested_empty")

            # Test with nested empty groups: [[], [1, 2]]
            @test_throws Exception eegfun.combine_conditions(
                "epochs_cleaned",
                [[], [1, 2]],
                input_dir = test_dir,
                output_dir = output_dir,
            )
        end

        @testset "Very large number of groups" begin
            output_dir = joinpath(test_dir, "combined_many_groups")

            # Test with many groups: [[1], [2], [3], [4], [1], [2], [3], [4]]
            result = eegfun.combine_conditions(
                "epochs_cleaned",
                [[1], [2], [3], [4], [1], [2], [3], [4]],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 3

            # Should create 8 groups (including duplicates)
            combined_epochs = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "epochs")
            @test length(combined_epochs) == 8
        end

        @testset "File corruption handling" begin
            corrupt_dir = joinpath(test_dir, "corrupt_files")
            mkpath(corrupt_dir)

            # Create a corrupted JLD2 file
            corrupt_file = joinpath(corrupt_dir, "1_epochs_corrupt.jld2")
            write(corrupt_file, "This is not a valid JLD2 file")

            output_dir = joinpath(test_dir, "combined_corrupt")
            result =
                eegfun.combine_conditions("epochs_corrupt", [[1, 2]], input_dir = corrupt_dir, output_dir = output_dir)

            # Should fail for the corrupted file
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Missing epochs variable" begin
            missing_var_dir = joinpath(test_dir, "missing_var")
            mkpath(missing_var_dir)

            # Create file with wrong variable name
            save(joinpath(missing_var_dir, "1_epochs_missing.jld2"), "wrong_var", create_test_epoch_data(2, 3))

            output_dir = joinpath(test_dir, "combined_missing_var")
            result = eegfun.combine_conditions(
                "epochs_missing",
                [[1, 2]],
                input_dir = missing_var_dir,
                output_dir = output_dir,
            )

            # Should fail for the file with missing epochs variable
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Empty epochs data" begin
            empty_epochs_dir = joinpath(test_dir, "empty_epochs")
            mkpath(empty_epochs_dir)

            # Create file with empty epochs array
            save(joinpath(empty_epochs_dir, "1_epochs_empty.jld2"), "epochs", eegfun.EpochData[])

            output_dir = joinpath(test_dir, "combined_empty_epochs")
            result = eegfun.combine_conditions(
                "epochs_empty",
                [[1, 2]],
                input_dir = empty_epochs_dir,
                output_dir = output_dir,
            )

            # Should fail because no conditions to combine
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Inconsistent sample rates across conditions" begin
            inconsistent_dir = joinpath(test_dir, "inconsistent_fs")
            mkpath(inconsistent_dir)

            # Create epochs with different sample rates
            fs1 = 256
            fs2 = 512
            n_samples = 101
            t1 = range(-0.2, 0.2, length = n_samples)
            t2 = range(-0.1, 0.1, length = n_samples)

            # Condition 1: 256 Hz
            dfs1 = DataFrame[]
            for ep = 1:2
                df = DataFrame(
                    time = collect(t1),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = sin.(2π .* 5 .* t1),
                )
                push!(dfs1, df)
            end

            # Condition 2: 512 Hz
            dfs2 = DataFrame[]
            for ep = 1:2
                df = DataFrame(
                    time = collect(t2),
                    sample = 1:n_samples,
                    condition = fill(2, n_samples),
                    condition_name = fill("condition_2", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = cos.(2π .* 5 .* t2),
                )
                push!(dfs2, df)
            end

            layout = eegfun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)

            epochs = [
                eegfun.EpochData(dfs1, layout, fs1, eegfun.AnalysisInfo()),
                eegfun.EpochData(dfs2, layout, fs2, eegfun.AnalysisInfo()),
            ]

            save(joinpath(inconsistent_dir, "1_epochs_inconsistent.jld2"), "epochs", epochs)

            output_dir = joinpath(test_dir, "combined_inconsistent")
            result = eegfun.combine_conditions(
                "epochs_inconsistent",
                [[1, 2]],
                input_dir = inconsistent_dir,
                output_dir = output_dir,
            )

            # Should succeed but use first condition's sample rate
            @test result.success == 1

            combined_epochs = load(joinpath(output_dir, "1_epochs_inconsistent.jld2"), "epochs")
            @test combined_epochs[1].sample_rate == fs1  # Should use first condition's sample rate
        end

        @testset "Inconsistent layouts across conditions" begin
            inconsistent_layout_dir = joinpath(test_dir, "inconsistent_layout")
            mkpath(inconsistent_layout_dir)

            # Create epochs with different layouts
            layout1 = eegfun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)
            layout2 = eegfun.Layout(DataFrame(label = [:Cz], inc = [0.0], azi = [0.0]), nothing, nothing)

            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Condition 1: Fz channel
            dfs1 = DataFrame[]
            for ep = 1:2
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = sin.(2π .* 5 .* t),
                )
                push!(dfs1, df)
            end

            # Condition 2: Cz channel
            dfs2 = DataFrame[]
            for ep = 1:2
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(2, n_samples),
                    condition_name = fill("condition_2", n_samples),
                    epoch = fill(ep, n_samples),
                    Cz = cos.(2π .* 5 .* t),
                )
                push!(dfs2, df)
            end

            epochs = [
                eegfun.EpochData(dfs1, layout1, fs, eegfun.AnalysisInfo()),
                eegfun.EpochData(dfs2, layout2, fs, eegfun.AnalysisInfo()),
            ]

            save(joinpath(inconsistent_layout_dir, "1_epochs_inconsistent_layout.jld2"), "epochs", epochs)

            output_dir = joinpath(test_dir, "combined_inconsistent_layout")
            result = eegfun.combine_conditions(
                "epochs_inconsistent_layout",
                [[1, 2]],
                input_dir = inconsistent_layout_dir,
                output_dir = output_dir,
            )

            # Should succeed but use first condition's layout
            @test result.success == 1

            combined_epochs = load(joinpath(output_dir, "1_epochs_inconsistent_layout.jld2"), "epochs")
            @test combined_epochs[1].layout.data == layout1.data  # Should use first condition's layout
        end

        @testset "Memory usage with large datasets" begin
            large_dir = joinpath(test_dir, "large_dataset")
            mkpath(large_dir)

            # Create a larger dataset (more epochs per condition)
            epochs = create_test_epoch_data(2, 50)  # 2 conditions, 50 epochs each
            save(joinpath(large_dir, "1_epochs_large.jld2"), "epochs", epochs)

            output_dir = joinpath(test_dir, "combined_large")
            result = eegfun.combine_conditions("epochs_large", [[1, 2]], input_dir = large_dir, output_dir = output_dir)

            @test result.success == 1

            # Verify the combined data
            combined_epochs = load(joinpath(output_dir, "1_epochs_large.jld2"), "epochs")
            @test length(combined_epochs) == 1
            @test length(combined_epochs[1].data) == 100  # 50 from each condition
        end

        @testset "Many channels" begin
            # Test with more channels (5)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)

            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Create 5 channels
            channel_names = Symbol.("Ch" .* string.(1:5))

            epochs = eegfun.EpochData[]
            for cond = 1:2
                dfs = DataFrame[]
                for ep = 1:3
                    df = DataFrame(
                        time = collect(t),
                        sample = 1:n_samples,
                        condition = fill(cond, n_samples),
                        condition_name = fill("condition_$cond", n_samples),
                        epoch = fill(ep, n_samples),
                    )

                    # Add channel data
                    for (i, ch) in enumerate(channel_names)
                        df[!, ch] = sin.(2π .* i .* t) .+ 0.1 .* randn(n_samples)
                    end

                    push!(dfs, df)
                end

                layout_df = DataFrame(label = channel_names, inc = zeros(5), azi = zeros(5))
                layout = eegfun.Layout(layout_df, nothing, nothing)

                push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
            end

            save(joinpath(many_ch_dir, "1_epochs_many.jld2"), "epochs", epochs)

            # Combine
            output_dir = joinpath(test_dir, "combined_many_ch")
            result =
                eegfun.combine_conditions("epochs_many", [[1, 2]], input_dir = many_ch_dir, output_dir = output_dir)

            @test result.success == 1

            combined_epochs = load(joinpath(output_dir, "1_epochs_many.jld2"), "epochs")

            # Verify all channels are present
            for ch in channel_names
                @test hasproperty(combined_epochs[1].data[1], ch)
            end

            # Verify epoch count
            @test length(combined_epochs[1].data) == 6  # 3 from each condition
        end

    finally
        # Cleanup
        rm(test_dir, recursive = true, force = true)
    end
end
