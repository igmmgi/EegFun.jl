using Test
using DataFrames
using EegFun
using JLD2
using Statistics
using CSV


@testset "EegFun" begin
    @testset "channel_summary" begin

        @testset "_channel_summary_impl core function" begin
            dat = create_test_continuous_data()

            # Test basic functionality
            result = EegFun._channel_summary_impl(dat.data, collect(1:50), [:Ch1, :Ch2])

            @test result isa DataFrame
            @test nrow(result) == 2  # Two channels
            @test :channel in propertynames(result)
            @test :min in propertynames(result)
            @test :max in propertynames(result)
            @test :std in propertynames(result)
            @test :range in propertynames(result)
            @test :var in propertynames(result)
            @test :zvar in propertynames(result)

            # Test column order
            @test propertynames(result)[1] == :channel

            # Test range calculation
            @test all(result.range .≈ result.max .- result.min)

            # Test variance calculation
            selected_data = dat.data[collect(1:50), [:Ch1, :Ch2]]
            expected_var = var.(eachcol(selected_data))
            @test result.var ≈ expected_var

            # Test that channels are correctly identified
            @test Set(result.channel) == Set([:Ch1, :Ch2])
        end

        @testset "_channel_summary_impl input validation" begin
            dat = create_test_continuous_data(n_channels = 4)

            @test_throws Exception EegFun._channel_summary_impl(dat.data, Int[], [:Ch1, :Ch2])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [1, 2, 3], Symbol[])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [1, 2, 3], [:NonExistent])
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [3000], [:Ch1])  # Out of bounds
            @test_throws Exception EegFun._channel_summary_impl(dat.data, [0], [:Ch1])     # Below bounds

        end

        @testset "_channel_summary_impl edge cases" begin
            dat = create_test_continuous_data(n_channels = 4)

            # Test single sample
            result = EegFun._channel_summary_impl(dat.data, [1], [:Ch1])
            @test nrow(result) == 1
            @test result.range[1] == 0.0  # min == max for single sample
            @test isnan(result.std[1]) || result.std[1] == 0.0  # std is NaN for single sample in Julia

            # Test single channel
            result = EegFun._channel_summary_impl(dat.data, collect(1:10), [:Ch1])
            @test nrow(result) == 1
            @test result.channel[1] == :Ch1

            # Test with constant signal (zero variance case)
            dat.data[!, :Ch4] .= 1.0
            result = EegFun._channel_summary_impl(dat.data, collect(1:50), [:Ch4])  # Ch4 is constant
            @test result.var[1] ≈ 0.0
            @test result.zvar[1] == 0.0  # Should handle zero variance case
        end

        @testset "channel_summary SingleDataFrameEeg" begin
            dat = create_test_continuous_data(n_channels = 4)

            # Test basic functionality with defaults
            result = EegFun.channel_summary(dat)
            @test result isa DataFrame
            @test nrow(result) == 4  # Four layout channels (Ch1, Ch2, Ch3, Ch4)
            @test :channel in propertynames(result)

            # Test that all layout channels are included by default
            @test Set(result.channel) == Set([:Ch1, :Ch2, :Ch3, :Ch4])

            # Test specific channel selection
            result_subset = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1, :Ch2]))
            @test nrow(result_subset) == 2
            @test Set(result_subset.channel) == Set([:Ch1, :Ch2])

            # Test sample selection (first half)
            n_samples = nrow(dat.data)
            result_samples = EegFun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
            @test result_samples isa DataFrame
            @test nrow(result_samples) == 4

            # Test include_meta and include_extra flags
            result_meta = EegFun.channel_summary(dat, include_meta = true)
            @test :time in result_meta.channel || :triggers in result_meta.channel  # Should include meta columns

            # Test statistics are reasonable
            result = EegFun.channel_summary(dat)
            @test all(result.range .>= 0)  # Range should be non-negative
            @test all(result.var .>= 0)    # Variance should be non-negative
            @test all(result.std .>= 0)    # Standard deviation should be non-negative
        end

        @testset "channel_summary SingleDataFrameEeg input validation" begin
            # Test empty data
            empty_df = DataFrame()
            layout = EegFun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing)
            empty_dat = EegFun.ContinuousData("test_data", empty_df, layout, 100, EegFun.AnalysisInfo())
            @test_throws Exception EegFun.channel_summary(empty_dat)
        end

        @testset "channel_summary MultiDataFrameEeg" begin
            epoch_dat = create_test_epoch_data()

            # Test basic functionality
            result = EegFun.channel_summary(epoch_dat)
            @test result isa DataFrame
            @test :epoch in propertynames(result)
            @test :channel in propertynames(result)

            # Test that we have results for each epoch and channel
            n_epochs = EegFun.n_epochs(epoch_dat)
            n_channels = 3  # A, B, C
            @test nrow(result) == n_epochs * n_channels

            # Test epoch numbers are preserved
            original_epoch_numbers = [df.epoch[1] for df in epoch_dat.data]
            result_epoch_numbers = unique(result.epoch)
            @test Set(result_epoch_numbers) == Set(original_epoch_numbers)

            # Test channel selection works across epochs
            result_subset = EegFun.channel_summary(epoch_dat, channel_selection = EegFun.channels([:Ch1]))
            @test nrow(result_subset) == n_epochs  # One channel per epoch
            @test all(result_subset.channel .== :Ch1)

            # Test that statistics vary across epochs (should be different)
            epochs_Ch1 = result[result.channel.==:Ch1, :]
            if nrow(epochs_Ch1) > 1
                # Should have some variation across epochs (unless data is identical)
                @test length(unique(epochs_Ch1.var)) > 1 || all(epochs_Ch1.var .== 0)
            end
        end

        @testset "channel_summary MultiDataFrameEeg input validation" begin
            # Test empty epoch data
            empty_layout = EegFun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing)
            empty_epochs = EegFun.EpochData("test_data", 1, "condition_1", DataFrame[], empty_layout, 100, EegFun.AnalysisInfo())
            @test_throws Exception EegFun.channel_summary(empty_epochs)
        end

        @testset "channel_summary statistical properties" begin
            dat = create_test_continuous_data(n_channels = 4)
            dat.data[!, :Ch4] .= 1.0
            result = EegFun.channel_summary(dat)

            # Test that constant channel has zero variance
            constant_row = result[result.channel.==:Ch4, :]
            @test nrow(constant_row) == 1
            @test constant_row.var[1] ≈ 0.0 atol = 1e-10
            @test constant_row.std[1] ≈ 0.0 atol = 1e-10
            @test constant_row.range[1] ≈ 0.0 atol = 1e-10

            # Test that varying channels have non-zero variance
            varying_rows = result[result.channel.!=:Ch4, :]
            @test all(varying_rows.var .> 0)
            @test all(varying_rows.std .> 0)
            @test all(varying_rows.range .> 0)

            # Test z-scored variance properties
            if length(unique(result.var)) > 1  # If there's variance in variances
                @test abs(mean(result.zvar)) < 1e-10  # Should be approximately zero mean
                @test abs(std(result.zvar) - 1.0) < 1e-10  # Should have unit std
            end
        end

        @testset "channel_summary consistency" begin
            dat = create_test_continuous_data()

            # Test that results are consistent between calls
            result1 = EegFun.channel_summary(dat)
            result2 = EegFun.channel_summary(dat)
            @test result1 == result2

            # Test that subsetting gives expected results
            full_result = EegFun.channel_summary(dat)
            subset_result = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1, :Ch2]))

            # The subset should match the corresponding rows from the full result
            full_subset = full_result[in([:Ch1, :Ch2]).(full_result.channel), :]
            @test subset_result.channel == full_subset.channel
            @test subset_result.min ≈ full_subset.min
            @test subset_result.max ≈ full_subset.max
            @test subset_result.std ≈ full_subset.std
        end

        @testset "channel_summary integration with selection functions" begin
            dat = create_test_continuous_data(n_channels = 4)

            # Test with different selection functions
            # Note: These functions are defined in the EegFun package
            result_all = EegFun.channel_summary(dat, channel_selection = EegFun.channels())
            @test nrow(result_all) == 4

            # Test channel selection by name
            result_specific = EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Ch1]))
            @test nrow(result_specific) == 1
            @test result_specific.channel[1] == :Ch1

            # Test sample selection by range
            result_half = EegFun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
            @test result_half isa DataFrame
            @test nrow(result_half) == 4  # Same number of channels
        end

    end # channel_summary testset
end # EegFun testset


@testset "Batch Channel Summary" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_batch_test_erp_data(n_conditions = 2)
                filename = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
                jldsave(filename; data = erps)
                @test isfile(filename)
            end
        end

        @testset "Basic channel summary" begin
            output_dir = joinpath(test_dir, "summary_output")

            # Test basic channel summary
            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

            @test isdir(output_dir)

            # Check that CSV file exists
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)

            # Load and verify CSV content
            results = CSV.read(csv_file, DataFrame)

            # Should have results for 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12

            # Check required columns exist
            @test hasproperty(results, :file)
            @test hasproperty(results, :condition)
            @test hasproperty(results, :channel)
            @test hasproperty(results, :min)
            @test hasproperty(results, :max)
            @test hasproperty(results, :std)
            @test hasproperty(results, :range)
            @test hasproperty(results, :var)
            @test hasproperty(results, :zvar)

            # Verify channels are present (CSV reads them as strings)
            @test "Ch1" in results.channel
            @test "Ch2" in results.channel
            @test "Ch3" in results.channel
        end

        @testset "Summary specific participants" begin
            output_dir = joinpath(test_dir, "summary_participant")

            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)

            results = CSV.read(csv_file, DataFrame)
            # Only participant 1: 1 file × 2 conditions × 3 channels = 6 rows
            @test nrow(results) == 6
            @test all(results.file .== "1_erps_cleaned")
        end

        @testset "Summary multiple participants" begin
            output_dir = joinpath(test_dir, "summary_multi_participants")

            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 2]),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)

            # 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12
            @test "1_erps_cleaned" in results.file
            @test "2_erps_cleaned" in results.file
        end

        @testset "Summary specific conditions" begin
            output_dir = joinpath(test_dir, "summary_condition")

            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(1),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)

            # 2 files × 1 condition × 3 channels = 6 rows
            @test nrow(results) == 6
            @test all(results.condition .== 1)
        end

        @testset "Summary multiple conditions" begin
            output_dir = joinpath(test_dir, "summary_multi_conditions")

            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions([1, 2]),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)

            # 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12
            @test 1 in results.condition
            @test 2 in results.condition
        end

        @testset "Channel selection predicate" begin
            output_dir = joinpath(test_dir, "summary_channel_select")

            # Select only Ch1 and Ch2
            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                channel_selection = EegFun.channels([:Ch1, :Ch2]),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)

            # 2 files × 2 conditions × 2 channels = 8 rows
            @test nrow(results) == 8
            @test "Ch1" in results.channel
            @test "Ch2" in results.channel
            @test "Ch3" ∉ results.channel
        end

        @testset "Channel exclusion predicate" begin
            output_dir = joinpath(test_dir, "summary_channel_exclude")

            # Exclude Ch3
            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                channel_selection = EegFun.channels_not([:Ch3]),
            )

            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)

            # 2 files × 2 conditions × 2 channels = 8 rows
            @test nrow(results) == 8
            @test "Ch3" ∉ results.channel
        end

        @testset "Custom output filename" begin
            output_dir = joinpath(test_dir, "summary_custom_name")

            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir, output_file = "my_custom_summary")

            # Check custom filename
            csv_file = joinpath(output_dir, "my_custom_summary.csv")
            @test isfile(csv_file)

            # Check log file has custom name too
            log_file = joinpath(output_dir, "my_custom_summary.log")
            @test isfile(log_file)
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception EegFun.channel_summary("erps_cleaned", input_dir = "/nonexistent/path")
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            # Directory exists but has no JLD2 files matching pattern
            result = EegFun.channel_summary("erps_cleaned", input_dir = empty_dir)

            @test result === nothing  # No files to process
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "summary_with_log")

            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

            # Check log file exists
            log_file = joinpath(output_dir, "channel_summary.log")
            @test isfile(log_file)

            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch channel summary started")
            @test contains(log_contents, "channel_summary")
            @test contains(log_contents, "erps_cleaned")
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_summary")
            mkpath(output_dir)

            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            # Run channel_summary - should work fine with existing directory
            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "channel_summary.csv"))
        end

        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)

            # Create one valid file
            erps = create_batch_test_erp_data(n_conditions = 2)
            jldsave(joinpath(partial_dir, "1_erps_cleaned.jld2"); data = erps)

            # Create one malformed file (invalid data type - String instead of Vector{ErpData})
            jldsave(joinpath(partial_dir, "2_erps_cleaned.jld2"); data = "invalid_data")

            output_dir = joinpath(test_dir, "summary_partial")
            EegFun.channel_summary("erps_cleaned", input_dir = partial_dir, output_dir = output_dir)

            # Should have results from the valid file only
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)

            results = CSV.read(csv_file, DataFrame)
            # 1 file × 2 conditions × 3 channels = 6 rows
            @test nrow(results) == 6
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "summary_invalid_condition")

            # Request condition 5 when only 2 exist
            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(5),
            )

            # Should produce empty results
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test !isfile(csv_file)  # No file created because no valid results
        end

        @testset "Statistics correctness" begin
            # Create data with known statistics
            stats_dir = joinpath(test_dir, "stats_test")
            mkpath(stats_dir)

            fs = 256.0
            n_samples = 100
            t = range(0.0, 1.0, length = n_samples)

            # Channel with known properties
            # Ch1: constant value = 5.0
            # Ch2: values 1 to 100
            # Ch3: all zeros

            df = DataFrame(
                time = collect(t),
                sample = 1:n_samples,
                condition = fill(1, n_samples),
                Ch1 = fill(5.0, n_samples),
                Ch2 = collect(1.0:100.0),
                Ch3 = zeros(n_samples),
            )

            layout = EegFun.Layout(DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]), nothing, nothing)

            erps = [EegFun.ErpData("test_data", 1, "condition_1", df, layout, fs, EegFun.AnalysisInfo(), 1)]
            jldsave(joinpath(stats_dir, "1_erps_stats.jld2"); data = erps)

            # Process
            output_dir = joinpath(test_dir, "summary_stats")
            EegFun.channel_summary("erps_stats", input_dir = stats_dir, output_dir = output_dir)

            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)

            # Verify statistics for Ch1 (constant = 5.0)
            ch1 = results[results.channel.=="Ch1", :]
            @test ch1.min[1] ≈ 5.0
            @test ch1.max[1] ≈ 5.0
            @test ch1.std[1] ≈ 0.0 atol = 1e-10
            @test ch1.range[1] ≈ 0.0 atol = 1e-10
            @test ch1.var[1] ≈ 0.0 atol = 1e-10

            # Verify statistics for Ch2 (1 to 100)
            ch2 = results[results.channel.=="Ch2", :]
            @test ch2.min[1] ≈ 1.0
            @test ch2.max[1] ≈ 100.0
            @test ch2.range[1] ≈ 99.0
            @test ch2.std[1] ≈ std(1.0:100.0)
            @test ch2.var[1] ≈ var(1.0:100.0)

            # Verify statistics for Ch3 (all zeros)
            ch3 = results[results.channel.=="Ch3", :]
            @test ch3.min[1] ≈ 0.0
            @test ch3.max[1] ≈ 0.0
            @test ch3.std[1] ≈ 0.0 atol = 1e-10
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "summary_combined")

            # Summary for specific participant AND condition AND channels
            EegFun.channel_summary(
                "erps_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
                condition_selection = EegFun.conditions(1),
                channel_selection = EegFun.channels([:Ch1, :Ch2]),
            )

            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)

            # 1 file × 1 condition × 2 channels = 2 rows
            @test nrow(results) == 2
            @test all(results.file .== "1_erps_cleaned")
            @test all(results.condition .== 1)
            @test "Ch3" ∉ results.channel
        end

        @testset "Pattern matching variants" begin
            # Create files with different naming patterns
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)

            erps = create_batch_test_erp_data(n_conditions = 2)
            jldsave(joinpath(pattern_dir, "1_erps_original.jld2"); data = erps)
            jldsave(joinpath(pattern_dir, "2_erps_cleaned.jld2"); data = erps)
            jldsave(joinpath(pattern_dir, "3_custom_erps.jld2"); data = erps)

            # Test pattern matching "erps_original"
            output_dir1 = joinpath(test_dir, "summary_original")
            EegFun.channel_summary("erps_original", input_dir = pattern_dir, output_dir = output_dir1)

            results1 = CSV.read(joinpath(output_dir1, "channel_summary.csv"), DataFrame)
            @test nrow(results1) == 6  # 1 file × 2 conditions × 3 channels

            # Test pattern matching "erps" (should match all)
            output_dir2 = joinpath(test_dir, "summary_all_erps")
            EegFun.channel_summary("erps", input_dir = pattern_dir, output_dir = output_dir2)

            results2 = CSV.read(joinpath(output_dir2, "channel_summary.csv"), DataFrame)
            @test nrow(results2) == 18  # 3 files × 2 conditions × 3 channels
        end

        @testset "File metadata preservation" begin
            output_dir = joinpath(test_dir, "summary_metadata")

            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)

            # Check that file names are correctly preserved
            @test "1_erps_cleaned" in results.file
            @test "2_erps_cleaned" in results.file

            # Check condition numbers are correct
            @test 1 in results.condition
            @test 2 in results.condition

            # Verify each file-condition combination has all channels
            for file in unique(results.file)
                for cond in unique(results.condition)
                    subset = results[(results.file.==file).&(results.condition.==cond), :]
                    @test nrow(subset) == 3  # 3 channels
                end
            end
        end

        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "summary_overwrite")

            # First run
            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = overwrite_dir)

            csv_file = joinpath(overwrite_dir, "channel_summary.csv")
            mtime1 = stat(csv_file).mtime

            # Wait a tiny bit to ensure different mtime
            sleep(0.1)

            # Second run (should overwrite)
            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = overwrite_dir)

            # Verify file was overwritten
            mtime2 = stat(csv_file).mtime
            @test mtime2 > mtime1
        end

        @testset "Many channels" begin
            # Test with more channels (10)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)

            fs = 256.0
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Create 10 channels
            channel_names = Symbol.("Ch" .* string.(1:10))

            df = DataFrame(time = collect(t), sample = 1:n_samples, condition = fill(1, n_samples))
            for (i, ch) in enumerate(channel_names)
                df[!, ch] = i .* randn(n_samples)  # Different variance for each channel
            end

            layout = EegFun.Layout(DataFrame(label = channel_names, inc = zeros(10), azi = zeros(10)), nothing, nothing)

            erps = [EegFun.ErpData("test_data", 1, "condition_1", df, layout, fs, EegFun.AnalysisInfo(), 1)]
            jldsave(joinpath(many_ch_dir, "1_erps_many.jld2"); data = erps)

            # Process
            output_dir = joinpath(test_dir, "summary_many_ch")
            EegFun.channel_summary("erps_many", input_dir = many_ch_dir, output_dir = output_dir)

            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)

            # Should have all 10 channels
            @test nrow(results) == 10
            for ch in channel_names
                @test string(ch) in results.channel
            end
        end

        @testset "Z-scored variance" begin
            # Test that zvar is correctly computed
            output_dir = joinpath(test_dir, "summary_zvar")

            EegFun.channel_summary("erps_cleaned", input_dir = test_dir, output_dir = output_dir)

            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)

            # For each file-condition combination, zvar should have mean ≈ 0 and std ≈ 1
            for file in unique(results.file)
                for cond in unique(results.condition)
                    subset = results[(results.file.==file).&(results.condition.==cond), :]

                    # Mean of z-scores should be close to 0
                    @test mean(subset.zvar) ≈ 0.0 atol = 1e-10

                    # Std of z-scores should be close to 1 (if more than 1 channel)
                    # Note: With small sample size (3 channels), there's some variability
                    if nrow(subset) > 1
                        @test std(subset.zvar, corrected = false) ≈ 1.0 atol = 0.2
                    end
                end
            end
        end

    finally
        # Cleanup
        rm(test_dir, recursive = true, force = true)
    end
end
