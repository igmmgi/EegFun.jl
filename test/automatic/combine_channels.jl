using Test
using DataFrames
using eegfun
using JLD2
using Statistics

@testset "Batch Combine Channels" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Helper to create test ErpData with multiple channels
        function create_test_erp_data(n_conditions::Int = 2)
            erps = eegfun.ErpData[]
            fs = 256.0
            n_samples = 513
            t = range(-1.0, 1.0, length = n_samples)

            for cond = 1:n_conditions
                # Create multiple channels with known relationships
                # Fp1, Fp2: frontal pair
                # PO7, PO8: parietal-occipital pair
                # Fz, Cz, Pz: midline electrodes
                fp1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                fp2 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                po7 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                po8 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                fz = 0.5 .* (fp1 .+ fp2)
                cz = randn(n_samples)
                pz = 0.5 .* (po7 .+ po8)

                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(cond, n_samples),
                    Fp1 = fp1,
                    Fp2 = fp2,
                    PO7 = po7,
                    PO8 = po8,
                    Fz = fz,
                    Cz = cz,
                    Pz = pz,
                )

                layout = eegfun.Layout(
                    DataFrame(
                        label = [:Fp1, :Fp2, :PO7, :PO8, :Fz, :Cz, :Pz],
                        inc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        azi = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ),
                    nothing,
                    nothing,
                )

                push!(erps, eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 10))
            end

            return erps
        end

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_test_erp_data(2)
                filename = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
                save(filename, "erps", erps)
                @test isfile(filename)
            end
        end

        @testset "Basic channel combining" begin
            output_dir = joinpath(test_dir, "combined_output")

            # Combine two channel groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that combined files exist
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))

            # Load and verify combined data
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 2  # 2 conditions
            @test erps[1] isa eegfun.ErpData

            # Should have original channels + 2 combined channels
            @test hasproperty(erps[1].data, :Fp1)
            @test hasproperty(erps[1].data, :Fp2)
            @test hasproperty(erps[1].data, :PO7)
            @test hasproperty(erps[1].data, :PO8)
            @test hasproperty(erps[1].data, :combined_1)  # Default label
            @test hasproperty(erps[1].data, :combined_2)  # Default label

            # Verify averaging: combined_1 should be mean of Fp1 and Fp2
            @test erps[1].data.combined_1 ≈ (erps[1].data.Fp1 .+ erps[1].data.Fp2) ./ 2
        end

        @testset "Custom output labels" begin
            output_dir = joinpath(test_dir, "combined_custom_labels")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)
            @test !hasproperty(erps[1].data, :combined_1)
        end

        @testset "Reduce mode" begin
            output_dir = joinpath(test_dir, "combined_reduce")

            # With reduce=true, only combined channels should remain
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Should have metadata + only combined channels
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)

            # Original channels should be removed
            @test !hasproperty(erps[1].data, :Fp1)
            @test !hasproperty(erps[1].data, :Fp2)
            @test !hasproperty(erps[1].data, :PO7)
            @test !hasproperty(erps[1].data, :PO8)
            @test !hasproperty(erps[1].data, :Fz)
        end

        @testset "Combine specific participants" begin
            output_dir = joinpath(test_dir, "combined_participant")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 1,
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Combine multiple participants" begin
            output_dir = joinpath(test_dir, "combined_multi_participants")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = [1, 2],
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Combine specific conditions" begin
            output_dir = joinpath(test_dir, "combined_condition")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = 1,
            )

            @test result.success == 2

            # Load and verify only one condition
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test erps[1].data[1, :condition] == 1
        end

        @testset "Combine multiple conditions" begin
            output_dir = joinpath(test_dir, "combined_multi_conditions")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = [1, 2],
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 2
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = "/nonexistent/path",
            )

            # Mismatched output labels
            @test_throws Exception eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal],  # Only 1 label for 2 groups
                input_dir = test_dir,
            )
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            result = eegfun.combine_channels("erps_cleaned", [eegfun.channels([:Fp1, :Fp2])], input_dir = empty_dir)

            @test result === nothing
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "combined_with_log")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            log_file = joinpath(output_dir, "combine_channels.log")
            @test isfile(log_file)

            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch channel combining started")
            @test contains(log_contents, "combine_channels")
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_combine")
            mkpath(output_dir)

            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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
            erps = create_test_erp_data(2)
            save(joinpath(partial_dir, "1_erps_cleaned.jld2"), "erps", erps)

            # Create one malformed file
            save(joinpath(partial_dir, "2_erps_cleaned.jld2"), "invalid_var", erps)

            output_dir = joinpath(test_dir, "combined_partial")
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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

        @testset "Combining math correctness" begin
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
                eegfun.Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(math_dir, "1_erps_math.jld2"), "erps", erps)

            # Combine
            output_dir = joinpath(test_dir, "combined_math")
            eegfun.combine_channels(
                "erps_math",
                [eegfun.channels([:Ch1, :Ch2])],
                output_labels = [:avg],
                input_dir = math_dir,
                output_dir = output_dir,
            )

            erps_combined = load(joinpath(output_dir, "1_erps_math.jld2"), "erps")

            # Expected average: [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]
            expected_avg = fill(6.0, 11)
            @test erps_combined[1].data.avg ≈ expected_avg
        end

        @testset "Layout preservation" begin
            output_dir = joinpath(test_dir, "combined_layout")

            # Get original layout info
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            original_n_channels = nrow(original_erps[1].layout.data)

            # Combine without reduce (should append to layout)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Layout should have original channels + 1 combined
            @test nrow(erps[1].layout.data) == original_n_channels + 1

            # Verify combined channel is in layout
            @test :combined_1 in erps[1].layout.data.label
        end

        @testset "Layout in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_layout_reduce")

            # Combine with reduce (layout should only have combined channels)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Layout should only have 2 channels
            @test nrow(erps[1].layout.data) == 2
            @test :frontal in erps[1].layout.data.label
            @test :parietal in erps[1].layout.data.label
            @test :Fp1 ∉ erps[1].layout.data.label
        end

        @testset "Multiple channel groups" begin
            output_dir = joinpath(test_dir, "combined_multiple_groups")

            # Combine 3 different groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8]), eegfun.channels([:Fz, :Cz, :Pz])],
                output_labels = [:frontal, :parietal, :midline],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)
            @test hasproperty(erps[1].data, :midline)

            # Verify midline is average of Fz, Cz, Pz
            expected_midline = (erps[1].data.Fz .+ erps[1].data.Cz .+ erps[1].data.Pz) ./ 3
            @test erps[1].data.midline ≈ expected_midline
        end

        @testset "Single channel group" begin
            output_dir = joinpath(test_dir, "combined_single_group")

            # Combine just one group
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "combined_filters")

            # Combine with participant AND condition filters
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 1,
                conditions = 1,
            )

            @test result.success == 1

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Pattern matching variants" begin
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)

            erps = create_test_erp_data(2)
            save(joinpath(pattern_dir, "1_erps_original.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "2_erps_cleaned.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "3_custom_erps.jld2"), "erps", erps)

            # Test pattern matching "erps_original"
            output_dir1 = joinpath(test_dir, "combined_original")
            result1 = eegfun.combine_channels(
                "erps_original",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = pattern_dir,
                output_dir = output_dir1,
            )
            @test result1.success == 1

            # Test pattern matching "erps" (should match all)
            output_dir2 = joinpath(test_dir, "combined_all_erps")
            result2 = eegfun.combine_channels(
                "erps",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = pattern_dir,
                output_dir = output_dir2,
            )
            @test result2.success == 3
        end

        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "combined_overwrite")

            # First run
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            file1 = joinpath(overwrite_dir, "1_erps_cleaned.jld2")
            mtime1 = stat(file1).mtime

            sleep(0.1)

            # Second run (should overwrite)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            mtime2 = stat(file1).mtime
            @test mtime2 > mtime1
        end

        @testset "Channel exclusion predicate" begin
            output_dir = joinpath(test_dir, "combined_exclusion")

            # Combine all channels except Fz
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels_not([:Fz])],
                output_labels = [:not_fz],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :not_fz)
            @test !hasproperty(erps[1].data, :Fz)
            @test !hasproperty(erps[1].data, :Fp1)  # All original channels removed in reduce mode
        end

        @testset "EpochData support" begin
            # Create EpochData test files
            epochs_dir = joinpath(test_dir, "epochs_test")
            mkpath(epochs_dir)

            function create_test_epoch_data()
                epochs = eegfun.EpochData[]
                fs = 256
                n_samples = 101
                t = range(-0.2, 0.2, length = n_samples)

                for cond = 1:2
                    dfs = DataFrame[]
                    for ep = 1:3
                        fp1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                        fp2 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)

                        df = DataFrame(
                            time = collect(t),
                            sample = 1:n_samples,
                            condition = fill(cond, n_samples),
                            condition_name = fill("condition_$cond", n_samples),
                            epoch = fill(ep, n_samples),
                            Fp1 = fp1,
                            Fp2 = fp2,
                        )
                        push!(dfs, df)
                    end

                    layout = eegfun.Layout(
                        DataFrame(label = [:Fp1, :Fp2], inc = [0.0, 0.0], azi = [0.0, 0.0]),
                        nothing,
                        nothing,
                    )

                    push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
                end

                return epochs
            end

            epochs = create_test_epoch_data()
            save(joinpath(epochs_dir, "1_epochs.jld2"), "epochs", epochs)

            # Combine channels in epoch data
            output_dir = joinpath(test_dir, "combined_epochs")
            result = eegfun.combine_channels(
                "epochs",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = epochs_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            epochs_combined = load(joinpath(output_dir, "1_epochs.jld2"), "epochs")
            @test epochs_combined[1] isa eegfun.EpochData
            @test hasproperty(epochs_combined[1].data[1], :frontal)
            @test hasproperty(epochs_combined[1].data[1], :Fp1)  # Original channels preserved
        end

        @testset "Invalid channel selection" begin
            # Request non-existent channels
            invalid_dir = joinpath(test_dir, "invalid_channels")
            mkpath(invalid_dir)

            erps = create_test_erp_data(1)
            save(joinpath(invalid_dir, "1_erps.jld2"), "erps", erps)

            output_dir = joinpath(test_dir, "combined_invalid")

            # This should fail because channels don't exist
            result = eegfun.combine_channels(
                "erps",
                [eegfun.channels([:NonExistent1, :NonExistent2])],
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
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fz])],
                output_labels = [:fz_avg],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :fz_avg)

            # Single channel average should equal the original
            @test erps[1].data.fz_avg ≈ erps[1].data.Fz
        end

        @testset "All channels combined" begin
            output_dir = joinpath(test_dir, "combined_all_channels")

            # Combine ALL channels
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels()],  # All channels
                output_labels = [:global_avg],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Should only have metadata + global_avg
            @test hasproperty(erps[1].data, :global_avg)
            @test !hasproperty(erps[1].data, :Fp1)

            # Verify it's the average of all original channels
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            all_ch_data = hcat(
                original_erps[1].data.Fp1,
                original_erps[1].data.Fp2,
                original_erps[1].data.PO7,
                original_erps[1].data.PO8,
                original_erps[1].data.Fz,
                original_erps[1].data.Cz,
                original_erps[1].data.Pz,
            )
            expected_avg = vec(mean(all_ch_data, dims = 2))
            @test erps[1].data.global_avg ≈ expected_avg
        end

        @testset "Overlapping channel groups" begin
            output_dir = joinpath(test_dir, "combined_overlapping")

            # Same channel in multiple groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2, :Fz]), eegfun.channels([:Fz, :Cz, :Pz])],
                output_labels = [:frontal_with_fz, :midline],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal_with_fz)
            @test hasproperty(erps[1].data, :midline)

            # Verify both averages include Fz
            expected_frontal = (erps[1].data.Fp1 .+ erps[1].data.Fp2 .+ erps[1].data.Fz) ./ 3
            expected_midline = (erps[1].data.Fz .+ erps[1].data.Cz .+ erps[1].data.Pz) ./ 3
            @test erps[1].data.frontal_with_fz ≈ expected_frontal
            @test erps[1].data.midline ≈ expected_midline
        end

        @testset "Sample rate and analysis info preservation" begin
            output_dir = joinpath(test_dir, "combined_metadata_preserve")

            # Get original metadata
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            original_fs = original_erps[1].sample_rate

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Verify sample rate preserved
            @test erps[1].sample_rate == original_fs

            # Verify it's still ErpData
            @test erps[1] isa eegfun.ErpData
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "combined_invalid_condition")

            # Request condition 5 when only 2 exist
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = 5,
            )

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 2
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

            layout = eegfun.Layout(DataFrame(label = channel_names, inc = zeros(20), azi = zeros(20)), nothing, nothing)

            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(many_ch_dir, "1_erps_many.jld2"), "erps", erps)

            # Combine into 4 groups
            output_dir = joinpath(test_dir, "combined_many_ch")
            result = eegfun.combine_channels(
                "erps_many",
                [
                    eegfun.channels(channel_names[1:5]),
                    eegfun.channels(channel_names[6:10]),
                    eegfun.channels(channel_names[11:15]),
                    eegfun.channels(channel_names[16:20]),
                ],
                output_labels = [:group1, :group2, :group3, :group4],
                input_dir = many_ch_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            erps_combined = load(joinpath(output_dir, "1_erps_many.jld2"), "erps")
            @test hasproperty(erps_combined[1].data, :group1)
            @test hasproperty(erps_combined[1].data, :group2)
            @test hasproperty(erps_combined[1].data, :group3)
            @test hasproperty(erps_combined[1].data, :group4)
        end

        @testset "Metadata columns preserved" begin
            output_dir = joinpath(test_dir, "combined_metadata_cols")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Verify metadata columns still exist
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :sample)
            @test hasproperty(erps[1].data, :condition)
        end

        @testset "Metadata columns in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_metadata_reduce")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Metadata should be preserved in reduce mode
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :condition)

            # Only combined channel and metadata
            @test hasproperty(erps[1].data, :frontal)
            @test !hasproperty(erps[1].data, :Fp1)
        end

        @testset "Empty result handling" begin
            empty_dir = joinpath(test_dir, "empty_result")
            mkpath(empty_dir)

            # Create file but select channels that result in empty selection
            erps = create_test_erp_data(1)
            save(joinpath(empty_dir, "1_erps.jld2"), "erps", erps)

            # Use a predicate that selects nothing
            output_dir = joinpath(test_dir, "combined_empty")
            result = eegfun.combine_channels(
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


@testset "Batch Combine Channels" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Helper to create test ErpData with multiple channels
        function create_test_erp_data(n_conditions::Int = 2)
            erps = eegfun.ErpData[]
            fs = 256.0
            n_samples = 513
            t = range(-1.0, 1.0, length = n_samples)

            for cond = 1:n_conditions
                # Create multiple channels with known relationships
                # Fp1, Fp2: frontal pair
                # PO7, PO8: parietal-occipital pair
                # Fz, Cz, Pz: midline electrodes
                fp1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                fp2 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                po7 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                po8 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                fz = 0.5 .* (fp1 .+ fp2)
                cz = randn(n_samples)
                pz = 0.5 .* (po7 .+ po8)

                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(cond, n_samples),
                    Fp1 = fp1,
                    Fp2 = fp2,
                    PO7 = po7,
                    PO8 = po8,
                    Fz = fz,
                    Cz = cz,
                    Pz = pz,
                )

                layout = eegfun.Layout(
                    DataFrame(
                        label = [:Fp1, :Fp2, :PO7, :PO8, :Fz, :Cz, :Pz],
                        inc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        azi = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ),
                    nothing,
                    nothing,
                )

                push!(erps, eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 10))
            end

            return erps
        end

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_test_erp_data(2)
                filename = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
                save(filename, "erps", erps)
                @test isfile(filename)
            end
        end

        @testset "Basic channel combining" begin
            output_dir = joinpath(test_dir, "combined_output")

            # Combine two channel groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that combined files exist
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))

            # Load and verify combined data
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 2  # 2 conditions
            @test erps[1] isa eegfun.ErpData

            # Should have original channels + 2 combined channels
            @test hasproperty(erps[1].data, :Fp1)
            @test hasproperty(erps[1].data, :Fp2)
            @test hasproperty(erps[1].data, :PO7)
            @test hasproperty(erps[1].data, :PO8)
            @test hasproperty(erps[1].data, :combined_1)  # Default label
            @test hasproperty(erps[1].data, :combined_2)  # Default label

            # Verify averaging: combined_1 should be mean of Fp1 and Fp2
            @test erps[1].data.combined_1 ≈ (erps[1].data.Fp1 .+ erps[1].data.Fp2) ./ 2
        end

        @testset "Custom output labels" begin
            output_dir = joinpath(test_dir, "combined_custom_labels")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)
            @test !hasproperty(erps[1].data, :combined_1)
        end

        @testset "Reduce mode" begin
            output_dir = joinpath(test_dir, "combined_reduce")

            # With reduce=true, only combined channels should remain
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Should have metadata + only combined channels
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)

            # Original channels should be removed
            @test !hasproperty(erps[1].data, :Fp1)
            @test !hasproperty(erps[1].data, :Fp2)
            @test !hasproperty(erps[1].data, :PO7)
            @test !hasproperty(erps[1].data, :PO8)
            @test !hasproperty(erps[1].data, :Fz)
        end

        @testset "Combine specific participants" begin
            output_dir = joinpath(test_dir, "combined_participant")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 1,
            )

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Combine multiple participants" begin
            output_dir = joinpath(test_dir, "combined_multi_participants")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = [1, 2],
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Combine specific conditions" begin
            output_dir = joinpath(test_dir, "combined_condition")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = 1,
            )

            @test result.success == 2

            # Load and verify only one condition
            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test erps[1].data[1, :condition] == 1
        end

        @testset "Combine multiple conditions" begin
            output_dir = joinpath(test_dir, "combined_multi_conditions")

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = [1, 2],
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 2
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = "/nonexistent/path",
            )

            # Mismatched output labels
            @test_throws Exception eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal],  # Only 1 label for 2 groups
                input_dir = test_dir,
            )
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            result = eegfun.combine_channels("erps_cleaned", [eegfun.channels([:Fp1, :Fp2])], input_dir = empty_dir)

            @test result === nothing
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "combined_with_log")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            log_file = joinpath(output_dir, "combine_channels.log")
            @test isfile(log_file)

            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch channel combining started")
            @test contains(log_contents, "combine_channels")
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_combine")
            mkpath(output_dir)

            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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
            erps = create_test_erp_data(2)
            save(joinpath(partial_dir, "1_erps_cleaned.jld2"), "erps", erps)

            # Create one malformed file
            save(joinpath(partial_dir, "2_erps_cleaned.jld2"), "invalid_var", erps)

            output_dir = joinpath(test_dir, "combined_partial")
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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

            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
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

        @testset "Combining math correctness" begin
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
                eegfun.Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(math_dir, "1_erps_math.jld2"), "erps", erps)

            # Combine
            output_dir = joinpath(test_dir, "combined_math")
            eegfun.combine_channels(
                "erps_math",
                [eegfun.channels([:Ch1, :Ch2])],
                output_labels = [:avg],
                input_dir = math_dir,
                output_dir = output_dir,
            )

            erps_combined = load(joinpath(output_dir, "1_erps_math.jld2"), "erps")

            # Expected average: [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0]
            expected_avg = fill(6.0, 11)
            @test erps_combined[1].data.avg ≈ expected_avg
        end

        @testset "Layout preservation" begin
            output_dir = joinpath(test_dir, "combined_layout")

            # Get original layout info
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            original_n_channels = nrow(original_erps[1].layout.data)

            # Combine without reduce (should append to layout)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Layout should have original channels + 1 combined
            @test nrow(erps[1].layout.data) == original_n_channels + 1

            # Verify combined channel is in layout
            @test :combined_1 in erps[1].layout.data.label
        end

        @testset "Layout in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_layout_reduce")

            # Combine with reduce (layout should only have combined channels)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8])],
                output_labels = [:frontal, :parietal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Layout should only have 2 channels
            @test nrow(erps[1].layout.data) == 2
            @test :frontal in erps[1].layout.data.label
            @test :parietal in erps[1].layout.data.label
            @test :Fp1 ∉ erps[1].layout.data.label
        end

        @testset "Multiple channel groups" begin
            output_dir = joinpath(test_dir, "combined_multiple_groups")

            # Combine 3 different groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2]), eegfun.channels([:PO7, :PO8]), eegfun.channels([:Fz, :Cz, :Pz])],
                output_labels = [:frontal, :parietal, :midline],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
            @test hasproperty(erps[1].data, :parietal)
            @test hasproperty(erps[1].data, :midline)

            # Verify midline is average of Fz, Cz, Pz
            expected_midline = (erps[1].data.Fz .+ erps[1].data.Cz .+ erps[1].data.Pz) ./ 3
            @test erps[1].data.midline ≈ expected_midline
        end

        @testset "Single channel group" begin
            output_dir = joinpath(test_dir, "combined_single_group")

            # Combine just one group
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal)
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "combined_filters")

            # Combine with participant AND condition filters
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 1,
                conditions = 1,
            )

            @test result.success == 1

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test !isfile(joinpath(output_dir, "2_erps_cleaned.jld2"))
        end

        @testset "Pattern matching variants" begin
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)

            erps = create_test_erp_data(2)
            save(joinpath(pattern_dir, "1_erps_original.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "2_erps_cleaned.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "3_custom_erps.jld2"), "erps", erps)

            # Test pattern matching "erps_original"
            output_dir1 = joinpath(test_dir, "combined_original")
            result1 = eegfun.combine_channels(
                "erps_original",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = pattern_dir,
                output_dir = output_dir1,
            )
            @test result1.success == 1

            # Test pattern matching "erps" (should match all)
            output_dir2 = joinpath(test_dir, "combined_all_erps")
            result2 = eegfun.combine_channels(
                "erps",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = pattern_dir,
                output_dir = output_dir2,
            )
            @test result2.success == 3
        end

        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "combined_overwrite")

            # First run
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            file1 = joinpath(overwrite_dir, "1_erps_cleaned.jld2")
            mtime1 = stat(file1).mtime

            sleep(0.1)

            # Second run (should overwrite)
            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = overwrite_dir,
            )

            mtime2 = stat(file1).mtime
            @test mtime2 > mtime1
        end

        @testset "Channel exclusion predicate" begin
            output_dir = joinpath(test_dir, "combined_exclusion")

            # Combine all channels except Fz
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels_not([:Fz])],
                output_labels = [:not_fz],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :not_fz)
            @test !hasproperty(erps[1].data, :Fz)
            @test !hasproperty(erps[1].data, :Fp1)  # All original channels removed in reduce mode
        end

        @testset "EpochData support" begin
            # Create EpochData test files
            epochs_dir = joinpath(test_dir, "epochs_test")
            mkpath(epochs_dir)

            function create_test_epoch_data()
                epochs = eegfun.EpochData[]
                fs = 256
                n_samples = 101
                t = range(-0.2, 0.2, length = n_samples)

                for cond = 1:2
                    dfs = DataFrame[]
                    for ep = 1:3
                        fp1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                        fp2 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)

                        df = DataFrame(
                            time = collect(t),
                            sample = 1:n_samples,
                            condition = fill(cond, n_samples),
                            condition_name = fill("condition_$cond", n_samples),
                            epoch = fill(ep, n_samples),
                            Fp1 = fp1,
                            Fp2 = fp2,
                        )
                        push!(dfs, df)
                    end

                    layout = eegfun.Layout(
                        DataFrame(label = [:Fp1, :Fp2], inc = [0.0, 0.0], azi = [0.0, 0.0]),
                        nothing,
                        nothing,
                    )

                    push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
                end

                return epochs
            end

            epochs = create_test_epoch_data()
            save(joinpath(epochs_dir, "1_epochs.jld2"), "epochs", epochs)

            # Combine channels in epoch data
            output_dir = joinpath(test_dir, "combined_epochs")
            result = eegfun.combine_channels(
                "epochs",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = epochs_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            epochs_combined = load(joinpath(output_dir, "1_epochs.jld2"), "epochs")
            @test epochs_combined[1] isa eegfun.EpochData
            @test hasproperty(epochs_combined[1].data[1], :frontal)
            @test hasproperty(epochs_combined[1].data[1], :Fp1)  # Original channels preserved
        end

        @testset "Invalid channel selection" begin
            # Request non-existent channels
            invalid_dir = joinpath(test_dir, "invalid_channels")
            mkpath(invalid_dir)

            erps = create_test_erp_data(1)
            save(joinpath(invalid_dir, "1_erps.jld2"), "erps", erps)

            output_dir = joinpath(test_dir, "combined_invalid")

            # This should fail because channels don't exist
            result = eegfun.combine_channels(
                "erps",
                [eegfun.channels([:NonExistent1, :NonExistent2])],
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
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fz])],
                output_labels = [:fz_avg],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :fz_avg)

            # Single channel average should equal the original
            @test erps[1].data.fz_avg ≈ erps[1].data.Fz
        end

        @testset "All channels combined" begin
            output_dir = joinpath(test_dir, "combined_all_channels")

            # Combine ALL channels
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels()],  # All channels
                output_labels = [:global_avg],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Should only have metadata + global_avg
            @test hasproperty(erps[1].data, :global_avg)
            @test !hasproperty(erps[1].data, :Fp1)

            # Verify it's the average of all original channels
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            all_ch_data = hcat(
                original_erps[1].data.Fp1,
                original_erps[1].data.Fp2,
                original_erps[1].data.PO7,
                original_erps[1].data.PO8,
                original_erps[1].data.Fz,
                original_erps[1].data.Cz,
                original_erps[1].data.Pz,
            )
            expected_avg = vec(mean(all_ch_data, dims = 2))
            @test erps[1].data.global_avg ≈ expected_avg
        end

        @testset "Overlapping channel groups" begin
            output_dir = joinpath(test_dir, "combined_overlapping")

            # Same channel in multiple groups
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2, :Fz]), eegfun.channels([:Fz, :Cz, :Pz])],
                output_labels = [:frontal_with_fz, :midline],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result.success == 2

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")
            @test hasproperty(erps[1].data, :frontal_with_fz)
            @test hasproperty(erps[1].data, :midline)

            # Verify both averages include Fz
            expected_frontal = (erps[1].data.Fp1 .+ erps[1].data.Fp2 .+ erps[1].data.Fz) ./ 3
            expected_midline = (erps[1].data.Fz .+ erps[1].data.Cz .+ erps[1].data.Pz) ./ 3
            @test erps[1].data.frontal_with_fz ≈ expected_frontal
            @test erps[1].data.midline ≈ expected_midline
        end

        @testset "Sample rate and analysis info preservation" begin
            output_dir = joinpath(test_dir, "combined_metadata_preserve")

            # Get original metadata
            original_erps = load(joinpath(test_dir, "1_erps_cleaned.jld2"), "erps")
            original_fs = original_erps[1].sample_rate

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Verify sample rate preserved
            @test erps[1].sample_rate == original_fs

            # Verify it's still ErpData
            @test erps[1] isa eegfun.ErpData
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "combined_invalid_condition")

            # Request condition 5 when only 2 exist
            result = eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = 5,
            )

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 2
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

            layout = eegfun.Layout(DataFrame(label = channel_names, inc = zeros(20), azi = zeros(20)), nothing, nothing)

            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(many_ch_dir, "1_erps_many.jld2"), "erps", erps)

            # Combine into 4 groups
            output_dir = joinpath(test_dir, "combined_many_ch")
            result = eegfun.combine_channels(
                "erps_many",
                [
                    eegfun.channels(channel_names[1:5]),
                    eegfun.channels(channel_names[6:10]),
                    eegfun.channels(channel_names[11:15]),
                    eegfun.channels(channel_names[16:20]),
                ],
                output_labels = [:group1, :group2, :group3, :group4],
                input_dir = many_ch_dir,
                output_dir = output_dir,
            )

            @test result.success == 1

            erps_combined = load(joinpath(output_dir, "1_erps_many.jld2"), "erps")
            @test hasproperty(erps_combined[1].data, :group1)
            @test hasproperty(erps_combined[1].data, :group2)
            @test hasproperty(erps_combined[1].data, :group3)
            @test hasproperty(erps_combined[1].data, :group4)
        end

        @testset "Metadata columns preserved" begin
            output_dir = joinpath(test_dir, "combined_metadata_cols")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = false,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Verify metadata columns still exist
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :sample)
            @test hasproperty(erps[1].data, :condition)
        end

        @testset "Metadata columns in reduce mode" begin
            output_dir = joinpath(test_dir, "combined_metadata_reduce")

            eegfun.combine_channels(
                "erps_cleaned",
                [eegfun.channels([:Fp1, :Fp2])],
                output_labels = [:frontal],
                input_dir = test_dir,
                output_dir = output_dir,
                reduce = true,
            )

            erps = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "erps")

            # Metadata should be preserved in reduce mode
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :condition)

            # Only combined channel and metadata
            @test hasproperty(erps[1].data, :frontal)
            @test !hasproperty(erps[1].data, :Fp1)
        end

        @testset "Empty result handling" begin
            empty_dir = joinpath(test_dir, "empty_result")
            mkpath(empty_dir)

            # Create file but select channels that result in empty selection
            erps = create_test_erp_data(1)
            save(joinpath(empty_dir, "1_erps.jld2"), "erps", erps)

            # Use a predicate that selects nothing
            output_dir = joinpath(test_dir, "combined_empty")
            result = eegfun.combine_channels(
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
