"""
Test suite for src/analysis/erp_measurements.jl
"""

using Test
using JLD2
using DataFrames
using CSV

@testset "Batch ERP Measurements" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic ERP measurements" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements")

        # Test mean amplitude measurement
        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 6  # 3 participants × 2 conditions
        @test ncol(result) >= 5  # participant, condition, condition_name, Ch1, Ch2, Ch3

        # Verify output file was created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test "erp_measurements.csv" in output_files

        # Load and verify CSV output
        csv_file = joinpath(output_dir, "erp_measurements.csv")
        csv_data = CSV.read(csv_file, DataFrame)
        @test nrow(csv_data) == 6
        @test "participant" in names(csv_data)
        @test "condition" in names(csv_data)
        @test "Ch1" in names(csv_data)
    end

    @testset "Different analysis types" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_types")

        for analysis_type in ["max_peak", "min_peak", "max_peak_lat", "min_peak_lat"]
            result = eegfun.erp_measurements(
                "erps_cleaned",
                analysis_type,
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result isa DataFrame
            @test nrow(result) == 6

            # Verify measurements are reasonable
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result, ch)
                    values = result[!, ch]
                    @test all(isfinite.(values))

                    if analysis_type == "max_peak"
                        @test all(values .> 0)  # Should be positive
                    elseif analysis_type == "min_peak"
                        @test all(isfinite.(values))  # Should be finite (may be positive or negative)
                    elseif analysis_type in ["max_peak_lat", "min_peak_lat"]
                        @test all(0.1 .<= values .<= 0.2)  # Should be in analysis window
                    end
                end
            end
        end
    end

    @testset "Area and integral measurements" begin
        # Create test ERP files
        for participant = 1:2
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]
            file_path = joinpath(test_dir, "$(participant)_erps_area.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_area")

        for analysis_type in ["rectified_area", "integral", "positive_area", "negative_area"]
            result = eegfun.erp_measurements(
                "erps_area",
                analysis_type,
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result isa DataFrame
            @test nrow(result) == 4  # 2 participants × 2 conditions

            # Verify measurements are finite
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result, ch)
                    values = result[!, ch]
                    @test all(isfinite.(values))

                    if analysis_type == "rectified_area"
                        @test all(values .>= 0)  # Should be non-negative
                    elseif analysis_type == "positive_area"
                        @test all(values .>= 0)  # Should be non-negative
                    elseif analysis_type == "negative_area"
                        @test all(values .>= 0)  # Should be non-negative (absolute value)
                    end
                end
            end
        end
    end

    @testset "Fractional latency measurements" begin
        # Create test ERP files
        for participant = 1:2
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]
            file_path = joinpath(test_dir, "$(participant)_erps_fractional.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_fractional")

        for analysis_type in ["fractional_area_lat", "fractional_peak_lat"]
            result = eegfun.erp_measurements(
                "erps_fractional",
                analysis_type,
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result isa DataFrame
            @test nrow(result) == 4  # 2 participants × 2 conditions

            # Verify measurements are finite and in time range
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result, ch)
                    values = result[!, ch]
                    @test all(isfinite.(values))
                    @test all(0.1 .<= values .<= 0.2)  # Should be in analysis window
                end
            end
        end
    end

    @testset "Measurement kwargs" begin
        # Create test ERP files
        for participant = 1:2
            erps = [create_test_erp_data(participant, 1)]
            file_path = joinpath(test_dir, "$(participant)_erps_kwargs.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_kwargs")

        # Test local_window
        result1 = eegfun.erp_measurements(
            "erps_kwargs",
            "max_peak",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            local_window = 5,
        )
        @test result1 isa DataFrame

        # Test fractional_area_fraction
        result2 = eegfun.erp_measurements(
            "erps_kwargs",
            "fractional_area_lat",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_area_fraction = 0.3,
        )
        @test result2 isa DataFrame

        # Test fractional_peak_fraction and fractional_peak_direction
        result3 = eegfun.erp_measurements(
            "erps_kwargs",
            "fractional_peak_lat",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_peak_fraction = 0.7,
            fractional_peak_direction = :offset,
        )
        @test result3 isa DataFrame
    end

    @testset "Kwargs validation" begin
        # Create test ERP files
        for participant = 1:2
            erps = [create_test_erp_data(participant, 1)]
            file_path = joinpath(test_dir, "$(participant)_erps_validate.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_validate")

        # Test invalid local_window
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "max_peak",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            local_window = 0,
        )

        # Test invalid fractional_area_fraction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_area_lat",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_area_fraction = 1.5,
        )

        # Test invalid fractional_peak_fraction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_peak_lat",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_peak_fraction = -0.1,
        )

        # Test invalid fractional_peak_direction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_peak_lat",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_peak_direction = :invalid,
        )
    end

    @testset "Epoch data processing" begin
        # Create test epoch files
        for participant = 1:2
            epochs = create_test_epoch_data(conditions = 2, n_channels = 3)  # This returns Vector{EpochData}

            file_path = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
            jldsave(file_path; data = epochs)
        end

        output_dir = joinpath(test_dir, "measurements_epochs")

        result = eegfun.erp_measurements(
            "epochs_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 40  # 2 participants × 2 conditions × 10 epochs each

        # Verify epoch column is present
        @test "epoch" in names(result)
    end

    @testset "Participant and condition filtering" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_filtered")

        # Test participant filtering
        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            participant_selection = eegfun.participants([1, 2]),
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 4  # 2 participants × 2 conditions
        @test all(result.participant .∈ [[1, 2]])

        # Test condition filtering
        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            condition_selection = eegfun.conditions([1]),
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 3  # 3 participants × 1 condition
        @test all(result.condition .== 1)
    end

    @testset "Channel selection" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_channels")

        # Test specific channel selection
        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            channel_selection = eegfun.channels([:Ch1, :Ch2]),
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 6
        @test "Ch1" in names(result)
        @test "Ch2" in names(result)
        @test "Ch3" ∉ names(result)  # Ch3 should be excluded
    end

    @testset "Baseline correction" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            jldsave(file_path; data = erps)
        end

        output_dir = joinpath(test_dir, "measurements_baseline")

        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            baseline_window = eegfun.samples((-0.2, 0.0)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa DataFrame
        @test nrow(result) == 6

        # Verify measurements are reasonable (should be different from non-baseline corrected)
        result_no_baseline = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        # Values should be different due to baseline correction
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(result, ch) && hasproperty(result_no_baseline, ch)
                @test !all(result[!, ch] .== result_no_baseline[!, ch])
            end
        end
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = "/nonexistent/dir",
            )
        end

        @testset "Invalid analysis type" begin
            @test_throws Exception eegfun.erp_measurements(
                "erps_cleaned",
                "invalid_type",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
            )
        end

        @testset "Invalid analysis window" begin
            @test_throws Exception eegfun.erp_measurements("erps_cleaned", "mean_amp", analysis_window = eegfun.samples((0.2, 0.1)), input_dir = test_dir)
        end

        @testset "Analysis window outside data range" begin
            output_dir = joinpath(test_dir, "measurements_outside")

            result = eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amp",
                analysis_window = eegfun.samples((2.1, 3.0)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should return nothing when no data in window
            @test result === nothing
        end
    end

    @testset "Edge cases" begin
        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "measurements_none")

            result = eegfun.erp_measurements(
                "nonexistent_pattern",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result === nothing
        end

        @testset "Empty data files" begin
            # Create file with empty data
            empty_file = joinpath(test_dir, "empty_erps_cleaned.jld2")
            jldsave(empty_file; data = eegfun.ErpData[])

            output_dir = joinpath(test_dir, "measurements_empty")

            result = eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                participant_selection = eegfun.participants(999),  # Non-existent participant
                output_dir = output_dir,
            )

            @test result === nothing
        end

        @testset "Files with no recognized data variable" begin
            # Create file with unrecognized variable
            unrecognized_file = joinpath(test_dir, "unrecognized_erps_cleaned.jld2")
            jldsave(unrecognized_file; other_data = "test")

            output_dir = joinpath(test_dir, "measurements_unrecognized")

            result = eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                participant_selection = eegfun.participants(999),  # Non-existent participant
                output_dir = output_dir,
            )

            @test result === nothing
        end

        @testset "Empty baseline window" begin
            # Create test ERP files
            for participant = 1:2
                erps = [create_test_erp_data(participant, 1)]
                file_path = joinpath(test_dir, "$(participant)_erps_baseline_edge.jld2")
                jldsave(file_path; data = erps)
            end

            output_dir = joinpath(test_dir, "measurements_baseline_edge")

            # Baseline window with no samples (outside data range)
            result = eegfun.erp_measurements(
                "erps_baseline_edge",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                baseline_window = eegfun.samples((10.0, 11.0)),  # No samples in this range
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should still work (baseline skipped if no samples)
            @test result isa DataFrame
        end

        @testset "Single sample analysis window" begin
            # Create test ERP files
            for participant = 1:2
                erps = [create_test_erp_data(participant, 1)]
                file_path = joinpath(test_dir, "$(participant)_erps_single.jld2")
                jldsave(file_path; data = erps)
            end

            output_dir = joinpath(test_dir, "measurements_single")

            # Very narrow window (might result in single sample)
            result = eegfun.erp_measurements(
                "erps_single",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.1001)),  # Very narrow window
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should handle gracefully
            if !isnothing(result)
                @test result isa DataFrame
            end
        end

        @testset "Empty channel selection" begin
            # Create test ERP files
            for participant = 1:2
                erps = [create_test_erp_data(participant, 1)]
                file_path = joinpath(test_dir, "$(participant)_erps_nochannels.jld2")
                jldsave(file_path; data = erps)
            end

            output_dir = joinpath(test_dir, "measurements_nochannels")

            # Select no channels
            result = eegfun.erp_measurements(
                "erps_nochannels",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                channel_selection = eegfun.channels(Symbol[]),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should return nothing
            @test result === nothing
        end
    end

    @testset "Output handling" begin
        @testset "Custom output directory and file" begin
            custom_dir = joinpath(test_dir, "custom_output")
            custom_file = "custom_measurements"

            result = eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amp",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = custom_dir,
                output_file = custom_file,
            )

            @test isdir(custom_dir)
            @test isfile(joinpath(custom_dir, "$(custom_file).csv"))
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.erp_measurements("erps_cleaned", "mean_amp", analysis_window = eegfun.samples((0.1, 0.2)), input_dir = test_dir)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "measurements_mean_amp")
            @test isdir(expected_dir)
        end
    end

    @testset "Data integrity and measurements" begin
        output_dir = joinpath(test_dir, "measurements_integrity")

        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            participant_selection = eegfun.participants(1),
            output_dir = output_dir,
        )

        # Verify participant ID extraction
        @test all(result.participant .== 1)

        # Verify condition numbers
        @test all(result.condition .∈ [[1, 2]])

        # Verify measurements are finite
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(result, ch)
                @test all(isfinite.(result[!, ch]))
            end
        end

        # Verify column ordering (metadata first, then channels)
        metadata_cols = [:participant, :condition, :condition_name]
        channel_cols = [:Ch1, :Ch2, :Ch3]

        for (i, col) in enumerate(metadata_cols)
            if hasproperty(result, col)
                @test names(result)[i] == string(col)
            end
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "measurements_logging")

        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amp",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        # Check that log file was created
        log_file = joinpath(output_dir, "erp_measurements.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("ERP measurement analysis started", log_content)
        @test occursin("Found 5 JLD2 files", log_content)
        @test occursin("Analysis complete", log_content)
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
