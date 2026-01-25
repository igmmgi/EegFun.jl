using Test
using DataFrames
using EegFun
using JLD2
using Statistics

@testset "channel_average" begin

    dat = create_test_data(n = 500)

    # 1) Append averaged columns only (Symbols input)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data)
    @test all(dat.data.Ch1_Ch2 .== (dat.data.Ch1 .+ dat.data.Ch2) ./ 2)

    # 2) Reduce to only averages and create averaged layout
    dat = create_test_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)
    @test all(propertynames(dat.data) .== [:time, :sample, :triggers, :Ch1_Ch2])
    @test size(dat.layout.data, 1) == 1
    @test :inc ∈ propertynames(dat.layout.data)
    @test :azi ∈ propertynames(dat.layout.data)

    # 3) Append averaged channels to layout
    dat = create_test_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    @test :Ch2_Ch3 ∈ propertynames(dat.data)
    @test any(dat.layout.data.label .== :Ch2_Ch3)

    # 4) Mixed Symbol input
    dat = create_test_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch3])])
    @test :Ch1_Ch3 ∈ propertynames(dat.data)
    @test :Ch1_Ch2 ∉ propertynames(dat.data)

    # 5) Auto-label :avg for all channels
    dat = create_test_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2, :Ch3])]; reduce = true)
    @test :avg ∈ propertynames(dat.data)

    # 6) Custom output_labels applied and length mismatch errors
    dat = create_test_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; output_labels = [:output_label])
    @test :output_label ∈ propertynames(dat.data)
    @test_throws Any EegFun.channel_average!(
        dat,
        channel_selections = [EegFun.channels([:A, :B]), EegFun.channels([:A, :C])];
        output_labels = [:x],
    )

    # 7) Duplicate labels in layout (should accumulate)
    # First add B_C, then add again to verify rows accumulate
    dat = create_test_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    n1 = sum(dat.layout.data.label .== :Ch2_Ch3)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch2, :Ch3])])
    n2 = sum(dat.layout.data.label .== :Ch2_Ch3)
    @test n1 == 1 && n2 == 2

    # 9) EpochData append and reduce
    # Create simple EpochData (2 epochs) from dat
    dat = create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])])

    @test :Ch1_Ch2 ∈ propertynames(dat.data[1]) && :Ch1_Ch2 ∈ propertynames(dat.data[2])
    @test :Ch1_Ch2 ∈ propertynames(dat.data[3]) && :Ch1_Ch2 ∈ propertynames(dat.data[4])

    dat = create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)

    # Expect meta columns (leading) + A_B only
    cols = propertynames(dat.data[1])
    @test cols[end] == :Ch1_Ch2
    @test :Ch1 ∉ cols && :Ch2 ∉ cols && :Ch3 ∉ cols

    # 10) ErpData reduce path
    dat = create_test_epoch_data(n = 500)
    dat = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])]; reduce = true)
    # condition and condition_name are now in struct, not DataFrame
    @test all(propertynames(dat.data[1]) .== [:time, :sample, :epoch, :Ch1_Ch2])
    @test all(propertynames(dat.data[end]) .== [:time, :sample, :epoch, :Ch1_Ch2])
    @test hasproperty(dat, :condition)
    @test hasproperty(dat, :condition_name)

    # 11) Test default behavior (average all channels)
    dat = create_test_epoch_data(n = 500)
    EegFun.channel_average!(dat)  # Should use default channels() from second function
    @test :avg ∈ propertynames(dat.data[1])
    @test all(dat.data[1].avg .== (dat.data[1].Ch1 .+ dat.data[1].Ch2 .+ dat.data[1].Ch3) ./ 3)
    @test all(dat.data[end].avg .== (dat.data[end].Ch1 .+ dat.data[end].Ch2 .+ dat.data[end].Ch3) ./ 3)

    # 12) Test single channel selection with custom label
    dat = create_test_epoch_data(n = 500)
    EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:custom])
    @test :custom ∈ propertynames(dat.data[1])
    @test :custom ∈ propertynames(dat.data[end])
    @test all(dat.data[1].custom .== (dat.data[1].Ch1 .+ dat.data[1].Ch2) ./ 2)
    @test all(dat.data[end].custom .== (dat.data[end].Ch1 .+ dat.data[end].Ch2) ./ 2)

end

@testset "Batch Average Channels" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Create test data files

        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_batch_test_erp_data(2; n_channels = 7)
                filename = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
                jldsave(filename; data = erps)
                @test isfile(filename)
            end
        end

        @testset "Basic channel averaging" begin
            output_dir = joinpath(test_dir, "combined_output")

            # Average two channel groups
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that averaged files exist
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))

            # Load and verify averaged data
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test length(erps) == 2  # 2 conditions
            @test erps[1] isa EegFun.ErpData

            # Should have original channels + 2 averaged channels
            @test hasproperty(erps[1].data, :Ch1)
            @test hasproperty(erps[1].data, :Ch2)
            @test hasproperty(erps[1].data, :Ch3)
            @test hasproperty(erps[1].data, :Ch4)
            @test hasproperty(erps[1].data, :avg_1)  # Default label
            @test hasproperty(erps[1].data, :avg_2)  # Default label

            # Verify averaging: avg_1 should be mean of Ch1 and Ch2
            @test erps[1].data.avg_1 ≈ (erps[1].data.Ch1 .+ erps[1].data.Ch2) ./ 2
        end

        @testset "Custom output labels" begin
            output_dir = joinpath(test_dir, "combined_custom_labels")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])],
                output_labels = [:Group1, :Group2],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :Group1)
            @test hasproperty(erps[1].data, :Group2)
            @test !hasproperty(erps[1].data, :avg_1)
        end

        @testset "Reduce mode" begin
            output_dir = joinpath(test_dir, "combined_reduce")

            # With reduce=true, only averaged channels should remain
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])],
                output_labels = [:Group1, :Group2],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Should have metadata + only averaged channels
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :Group1)
            @test hasproperty(erps[1].data, :Group2)

            # Original channels should be removed
            @test !hasproperty(erps[1].data, :Ch1)
            @test !hasproperty(erps[1].data, :Ch2)
            @test !hasproperty(erps[1].data, :Ch3)
            @test !hasproperty(erps[1].data, :Ch4)
            @test !hasproperty(erps[1].data, :Ch5)
        end

        @testset "Average specific participants" begin
            output_dir = joinpath(test_dir, "combined_participant")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Average multiple participants" begin
            output_dir = joinpath(test_dir, "combined_multi_participants")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 2]),
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Average specific conditions" begin
            output_dir = joinpath(test_dir, "combined_condition")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(1),
            )

            @test result.success == 2

            # Load and verify only one condition
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test length(erps) == 1
            @test erps[1].data[1, :condition] == 1
        end

        @testset "Average multiple conditions" begin
            output_dir = joinpath(test_dir, "combined_multi_conditions")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions([1, 2]),
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test length(erps) == 2
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = "/nonexistent/path",
            )

            # Mismatched output labels
            @test_throws Exception EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])],
                output_labels = [:Group1],  # Only 1 label for 2 groups
                input_dir = test_dir,
            )
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            result = EegFun.channel_average("erps_cleaned", [EegFun.channels([:Ch1, :Ch2])], input_dir = empty_dir)

            @test result === nothing
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "combined_with_log")

            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            log_file = joinpath(output_dir, "channel_average.log")
            @test isfile(log_file)

            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch channel averaging started")
            @test contains(log_contents, "channel_average")
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_combine")
            mkpath(output_dir)

            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test isfile(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
        end

        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)

            # Create one valid file
            erps = create_batch_test_erp_data(2; n_channels = 7)
            jldsave(joinpath(partial_dir, "1_erps_cleaned.jld2"); data = erps)

            # Create one malformed file (invalid data type - String instead of Vector{ErpData})
            jldsave(joinpath(partial_dir, "2_erps_cleaned.jld2"); data = "invalid_data")

            output_dir = joinpath(test_dir, "combined_partial")
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = partial_dir,
                output_dir = output_dir,
            )

            @test result.success == 1
            @test result.errors == 1
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "combined_return_check")

            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test hasfield(typeof(result), :success)
            @test hasfield(typeof(result), :errors)
            @test result.success isa Integer
            @test result.errors isa Integer
            @test result.success >= 0
            @test result.errors >= 0
            @test result.success + result.errors == 2
        end

        @testset "Averaging math correctness" begin
            # Verify that averaging is mathematically correct
            math_dir = joinpath(test_dir, "math_test")
            mkpath(math_dir)

            fs = 256.0
            n_samples = 11
            t = range(0.0, 1.0, length = n_samples)

            # Create channels with known values
            ch1_vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
            ch2_vals = [11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]

            df = DataFrame(
                time = collect(t),
                sample = 1:n_samples,
                condition = fill(1, n_samples),
                Ch1 = ch1_vals,
                Ch2 = ch2_vals,
            )

            layout =
                EegFun.Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

            erps = [EegFun.ErpData("test_data", 1, "condition_1", df, layout, fs, EegFun.AnalysisInfo(), 1)]
            jldsave(joinpath(math_dir, "1_erps_math.jld2"); data = erps)

            # Average
            output_dir = joinpath(test_dir, "combined_math")
            EegFun.channel_average(
                "erps_math",
                [EegFun.channels([:Ch1, :Ch2])],
                output_labels = [:avg],
                input_dir = math_dir,
                output_dir = output_dir,
            )

            erps_combined = load(joinpath(output_dir, "1_erps_math.jld2"), "data")

            # Expected average: [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]
            expected_avg = fill(6.0, 11)
            @test erps_combined[1].data.avg ≈ expected_avg
        end

        @testset "Layout preservation" begin
            output_dir = joinpath(test_dir, "combined_layout")

            # Get original layout info
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "data")
            original_n_channels = nrow(original_erps[1].layout.data)

            # Average without reduce (should append to layout)
            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Layout should have original channels + 1 averaged
            @test nrow(erps[1].layout.data) == original_n_channels + 1

            # Verify averaged channel is in layout
            @test :avg_1 in erps[1].layout.data.label
        end

        @testset "Layout in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_layout_reduce")

            # Average with reduce (layout should only have averaged channels)
            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4])],
                output_labels = [:Group1, :Group2],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Layout should only have 2 channels
            @test nrow(erps[1].layout.data) == 2
            @test :Group1 in erps[1].layout.data.label
            @test :Group2 in erps[1].layout.data.label
            @test :Ch1 ∉ erps[1].layout.data.label
        end

        @testset "Multiple channel groups" begin
            output_dir = joinpath(test_dir, "combined_multiple_groups")

            # Average 3 different groups
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2]), EegFun.channels([:Ch3, :Ch4]), EegFun.channels([:Ch5, :Ch6, :Ch7])],
                output_labels = [:Group1, :Group2, :Group3],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :Group1)
            @test hasproperty(erps[1].data, :Group2)
            @test hasproperty(erps[1].data, :Group3)

            expected_group3 = (erps[1].data.Ch5 .+ erps[1].data.Ch6 .+ erps[1].data.Ch7) ./ 3
            @test erps[1].data.Group3 ≈ expected_group3
        end

        @testset "Single channel group" begin
            output_dir = joinpath(test_dir, "combined_single_group")

            # Average just one group
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                output_labels = [:Group1],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :Group1)
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "combined_filters")

            # Average with participant AND condition filters
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participant_selection = EegFun.participants(1),
                condition_selection = EegFun.conditions(1),
            )

            @test result.success == 1

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test length(erps) == 1
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Pattern matching variants" begin
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)

            erps = create_batch_test_erp_data(2; n_channels = 7)
            jldsave(joinpath(pattern_dir, "1_erps_original.jld2"); data = erps)
            jldsave(joinpath(pattern_dir, "2_erps_cleaned.jld2"); data = erps)
            jldsave(joinpath(pattern_dir, "3_custom_erps.jld2"); data = erps)

            # Test pattern matching "erps_original"
            output_dir1 = joinpath(test_dir, "combined_original")
            result1 = EegFun.channel_average(
                "erps_original",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = pattern_dir,
                output_dir = output_dir1,
            )
            @test result1.success == 1

            # Test pattern matching "erps" (should match all)
            output_dir2 = joinpath(test_dir, "combined_all_erps")
            result2 = EegFun.channel_average(
                "erps",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = pattern_dir,
                output_dir = output_dir2,
            )
            @test result2.success == 3
        end

        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "combined_overwrite")

            # First run
            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            file1 = joinpath(overwrite_dir, "1_erps_cleaned.jld2")
            mtime1 = stat(file1).mtime

            sleep(0.1)

            # Second run (should overwrite)
            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            mtime2 = stat(file1).mtime
            @test mtime2 > mtime1
        end

        @testset "Channel exclusion predicate" begin
            output_dir = joinpath(test_dir, "combined_exclusion")

            # Average all channels except Ch5
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels_not([:Ch5])],
                output_labels = [:not_Ch5],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :not_Ch5)
            @test !hasproperty(erps[1].data, :Ch5)
            @test !hasproperty(erps[1].data, :Ch1)  # All original channels removed in reduce mode
        end

        @testset "EpochData support" begin
            # Create EpochData test files
            epochs_dir = joinpath(test_dir, "epochs_test")
            mkpath(epochs_dir)

            epochs = create_test_epoch_data()
            jldsave(joinpath(epochs_dir, "1_epochs.jld2"); data = [epochs])

            # Average channels in epoch data
            output_dir = joinpath(test_dir, "combined_epochs")
            result = EegFun.channel_average(
                "epochs",
                [EegFun.channels([:Ch1, :Ch2])],
                output_labels = [:Group1],
                input_dir = epochs_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            epochs_combined = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test epochs_combined[1] isa EegFun.EpochData
            @test hasproperty(epochs_combined[1].data[1], :Group1)
            @test hasproperty(epochs_combined[1].data[1], :Ch1)  # Original channels preserved
        end

        @testset "Invalid channel selection" begin
            # Request non-existent channels
            invalid_dir = joinpath(test_dir, "invalid_channels")
            mkpath(invalid_dir)

            erps = create_batch_test_erp_data(1; n_channels = 7)
            jldsave(joinpath(invalid_dir, "1_erps.jld2"); data = erps)

            output_dir = joinpath(test_dir, "combined_invalid")

            # This should fail because channels don't exist
            result = EegFun.channel_average(
                "erps",
                [EegFun.channels([:NonExistent1, :NonExistent2])],
                input_dir = invalid_dir,
                output_dir = output_dir,
            )

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Single channel in group" begin
            output_dir = joinpath(test_dir, "combined_single_channel")

            # Average just one channel (edge case)
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch5])],
                output_labels = [:Ch5_avg],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :Ch5_avg)

            # Single channel average should equal the original
            @test erps[1].data.Ch5_avg ≈ erps[1].data.Ch5
        end

        @testset "All channels combined" begin
            output_dir = joinpath(test_dir, "combined_all_channels")

            # Average ALL channels
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels()],  # All channels
                output_labels = [:global_avg],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Should only have metadata + global_avg
            @test hasproperty(erps[1].data, :global_avg)
            @test !hasproperty(erps[1].data, :Ch1)

            # Verify it's the average of all original channels
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "data")
            all_ch_data = hcat(
                original_erps[1].data.Ch1,
                original_erps[1].data.Ch2,
                original_erps[1].data.Ch3,
                original_erps[1].data.Ch4,
                original_erps[1].data.Ch5,
                original_erps[1].data.Ch6,
                original_erps[1].data.Ch7,
            )
            expected_avg = vec(mean(all_ch_data, dims = 2))
            @test erps[1].data.global_avg ≈ expected_avg
        end

        @testset "Overlapping channel groups" begin
            output_dir = joinpath(test_dir, "combined_overlapping")

            # Same channel in multiple groups
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2, :Ch5]), EegFun.channels([:Ch5, :Ch6, :Ch7])],
                output_labels = [:Group1, :Group2],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")
            @test hasproperty(erps[1].data, :Group1)
            @test hasproperty(erps[1].data, :Group2)

            # Verify both averages include Ch5
            expected_group1 = (erps[1].data.Ch1 .+ erps[1].data.Ch2 .+ erps[1].data.Ch5) ./ 3
            expected_group2 = (erps[1].data.Ch5 .+ erps[1].data.Ch6 .+ erps[1].data.Ch7) ./ 3
            @test erps[1].data.Group1 ≈ expected_group1
            @test erps[1].data.Group2 ≈ expected_group2
        end

        @testset "Sample rate and analysis info preservation" begin
            output_dir = joinpath(test_dir, "combined_metadata_preserve")

            # Get original metadata
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "data")
            original_fs = original_erps[1].sample_rate

            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Verify sample rate preserved
            @test erps[1].sample_rate == original_fs

            # Verify it's still ErpData
            @test erps[1] isa EegFun.ErpData
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "combined_invalid_condition")

            # Request condition 5 when only 2 exist
            # With predicate-based selection, this results in empty selection but successful processing
            result = EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions(5),
            )

            # Files are processed successfully but with empty condition selection
            @test result.success == 2
            @test result.errors == 0
        end

        @testset "Many channels" begin
            # Test with 20 channels
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)

            fs = 256.0
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Create 20 channels
            channel_names = Symbol.("Ch" .* string.(1:20))

            df = DataFrame(time = collect(t), sample = 1:n_samples, condition = fill(1, n_samples))
            for ch in channel_names
                df[!, ch] = randn(n_samples)
            end

            layout = EegFun.Layout(DataFrame(label = channel_names, inc = zeros(20), azi = zeros(20)), nothing, nothing)

            erps = [EegFun.ErpData("test_data", 1, "condition_1", df, layout, fs, EegFun.AnalysisInfo(), 1)]
            jldsave(joinpath(many_ch_dir, "1_erps_many.jld2"); data = erps)

            # Average into 4 groups
            output_dir = joinpath(test_dir, "combined_many_ch")
            result = EegFun.channel_average(
                "erps_many",
                [
                    EegFun.channels(channel_names[1:5]),
                    EegFun.channels(channel_names[6:10]),
                    EegFun.channels(channel_names[11:15]),
                    EegFun.channels(channel_names[16:20]),
                ],
                output_labels = [:Group1, :Group2, :Group3, :Group4],
                input_dir = many_ch_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            erps_combined = load(joinpath(output_dir, "1_erps_many.jld2"), "data")
            @test hasproperty(erps_combined[1].data, :Group1)
            @test hasproperty(erps_combined[1].data, :Group2)
            @test hasproperty(erps_combined[1].data, :Group3)
            @test hasproperty(erps_combined[1].data, :Group4)
        end

        @testset "Metadata columns preserved" begin
            output_dir = joinpath(test_dir, "combined_metadata_cols")

            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Verify metadata columns still exist
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :sample)
            # condition is now in struct, not DataFrame
            @test hasproperty(erps[1], :condition)
        end

        @testset "Metadata columns in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_metadata_reduce")

            EegFun.channel_average(
                "erps_cleaned",
                [EegFun.channels([:Ch1, :Ch2])],
                output_labels = [:Group1],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "data")

            # Metadata should be preserved in reduce mode
            @test hasproperty(erps[1].data, :time)
            # condition is now in struct, not DataFrame
            @test hasproperty(erps[1], :condition)

            # Only averaged channel and metadata
            @test hasproperty(erps[1].data, :Group1)
            @test !hasproperty(erps[1].data, :Ch1)
        end

        @testset "Empty result handling" begin
            empty_dir = joinpath(test_dir, "empty_result")
            mkpath(empty_dir)

            # Create file but select channels that result in empty selection
            erps = create_batch_test_erp_data(1; n_channels = 7)
            jldsave(joinpath(empty_dir, "1_erps.jld2"); data = erps)

            # Use a predicate that selects nothing
            output_dir = joinpath(test_dir, "combined_empty")
            result = EegFun.channel_average(
                "erps",
                [ch -> Symbol[]],  # Returns empty
                input_dir = empty_dir,
                output_dir = output_dir,
            )

            # Should fail because no channels selected
            @test result.success == 0
            @test result.errors == 1
        end

    finally
        # Cleanup
        rm(test_dir, recursive = true, force = true)
    end
end
