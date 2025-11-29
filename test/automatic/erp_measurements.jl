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
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 6  # 3 participants × 2 conditions
        @test ncol(result.data) >= 5  # participant, condition, condition_name, Ch1, Ch2, Ch3

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

        for analysis_type in ["max_peak_amplitude", "min_peak_amplitude", "max_peak_latency", "min_peak_latency", "peak_to_peak_amplitude", "peak_to_peak_latency"]
            result = eegfun.erp_measurements(
                "erps_cleaned",
                analysis_type,
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result isa eegfun.ErpMeasurementsResult
            @test nrow(result.data) == 6

            # Verify measurements are reasonable
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result.data, ch)
                    values = result.data[!, ch]
                    @test all(isfinite.(values))

                    if analysis_type == "max_peak_amplitude"
                        @test all(values .> 0)  # Should be positive
                    elseif analysis_type == "min_peak_amplitude"
                        @test all(isfinite.(values))  # Should be finite (may be positive or negative)
                    elseif analysis_type in ["max_peak_latency", "min_peak_latency"]
                        @test all(0.1 .<= values .<= 0.2)  # Should be in analysis window
                    elseif analysis_type == "peak_to_peak_latency"
                        @test all(values .>= 0)  # Should be non-negative (time difference)
                        @test all(values .<= 0.1)  # Should be <= analysis window duration (0.2 - 0.1 = 0.1)
                    elseif analysis_type == "peak_to_peak_amplitude"
                        @test all(values .> 0)  # Should be positive (difference between max and min)
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

            @test result isa eegfun.ErpMeasurementsResult
            @test nrow(result.data) == 4  # 2 participants × 2 conditions

            # Verify measurements are finite
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result.data, ch)
                    values = result.data[!, ch]
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

        for analysis_type in ["fractional_area_latency", "fractional_peak_latency"]
            result = eegfun.erp_measurements(
                "erps_fractional",
                analysis_type,
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = output_dir,
            )

            @test result isa eegfun.ErpMeasurementsResult
            @test nrow(result.data) == 4  # 2 participants × 2 conditions

            # Verify measurements are finite and in time range
            for ch in [:Ch1, :Ch2, :Ch3]
                if hasproperty(result.data, ch)
                    values = result.data[!, ch]
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
            "max_peak_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            local_window = 5,
        )
        @test result1 isa eegfun.ErpMeasurementsResult

        # Test fractional_area_fraction
        result2 = eegfun.erp_measurements(
            "erps_kwargs",
            "fractional_area_latency",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_area_fraction = 0.3,
        )
        @test result2 isa eegfun.ErpMeasurementsResult

        # Test fractional_peak_fraction and fractional_peak_direction
        result3 = eegfun.erp_measurements(
            "erps_kwargs",
            "fractional_peak_latency",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_peak_fraction = 0.7,
            fractional_peak_direction = :offset,
        )
        @test result3 isa eegfun.ErpMeasurementsResult
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
            "max_peak_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            local_window = 0,
        )

        # Test invalid fractional_area_fraction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_area_latency",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_area_fraction = 1.5,
        )

        # Test invalid fractional_peak_fraction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_peak_latency",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
            fractional_peak_fraction = -0.1,
        )

        # Test invalid fractional_peak_direction
        @test_throws Exception eegfun.erp_measurements(
            "erps_validate",
            "fractional_peak_latency",
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
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 40  # 2 participants × 2 conditions × 10 epochs each

        # Verify epoch column is present
        @test "epoch" in names(result.data)
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
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            participant_selection = eegfun.participants([1, 2]),
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 4  # 2 participants × 2 conditions
        @test all(result.data.participant .∈ [[1, 2]])

        # Test condition filtering
        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            condition_selection = eegfun.conditions([1]),
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 3  # 3 participants × 1 condition
        @test all(result.data.condition .== 1)
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
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            channel_selection = eegfun.channels([:Ch1, :Ch2]),
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 6
        @test "Ch1" in names(result.data)
        @test "Ch2" in names(result.data)
        @test "Ch3" ∉ names(result.data)  # Ch3 should be excluded
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
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            baseline_window = eegfun.samples((-0.2, 0.0)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        @test result isa eegfun.ErpMeasurementsResult
        @test nrow(result.data) == 6

        # Verify measurements are reasonable (should be different from non-baseline corrected)
        result_no_baseline = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            output_dir = output_dir,
        )

        # Values should be different due to baseline correction
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(result.data, ch) && hasproperty(result_no_baseline.data, ch)
                @test !all(result.data[!, ch] .== result_no_baseline.data[!, ch])
            end
        end
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amplitude",
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
            @test_throws Exception eegfun.erp_measurements("erps_cleaned", "mean_amplitude", analysis_window = eegfun.samples((0.2, 0.1)), input_dir = test_dir)
        end

        @testset "Analysis window outside data range" begin
            output_dir = joinpath(test_dir, "measurements_outside")

            result = eegfun.erp_measurements(
                "erps_cleaned",
                "mean_amplitude",
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
                "mean_amplitude",
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
                "mean_amplitude",
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
                "mean_amplitude",
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
                "mean_amplitude",
                analysis_window = eegfun.samples((0.1, 0.2)),
                baseline_window = eegfun.samples((10.0, 11.0)),  # No samples in this range
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should still work (baseline skipped if no samples)
            @test result isa eegfun.ErpMeasurementsResult
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
                "mean_amplitude",
                analysis_window = eegfun.samples((0.1, 0.1001)),  # Very narrow window
                input_dir = test_dir,
                output_dir = output_dir,
            )

            # Should handle gracefully
            if !isnothing(result)
                @test result isa eegfun.ErpMeasurementsResult
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
                "mean_amplitude",
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
                "mean_amplitude",
                analysis_window = eegfun.samples((0.1, 0.2)),
                input_dir = test_dir,
                output_dir = custom_dir,
                output_file = custom_file,
            )

            @test isdir(custom_dir)
            @test isfile(joinpath(custom_dir, "$(custom_file).csv"))
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.erp_measurements("erps_cleaned", "mean_amplitude", analysis_window = eegfun.samples((0.1, 0.2)), input_dir = test_dir)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "measurements_mean_amplitude")
            @test isdir(expected_dir)
        end
    end

    @testset "Data integrity and measurements" begin
        output_dir = joinpath(test_dir, "measurements_integrity")

        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amplitude",
            analysis_window = eegfun.samples((0.1, 0.2)),
            input_dir = test_dir,
            participant_selection = eegfun.participants(1),
            output_dir = output_dir,
        )

        # Verify participant ID extraction
        @test all(result.data.participant .== 1)

        # Verify condition numbers
        @test all(result.data.condition .∈ [[1, 2]])

        # Verify measurements are finite
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(result.data, ch)
                @test all(isfinite.(result.data[!, ch]))
            end
        end

        # Verify column ordering (metadata first, then channels)
        metadata_cols = [:participant, :condition, :condition_name]
        channel_cols = [:Ch1, :Ch2, :Ch3]

        for (i, col) in enumerate(metadata_cols)
            if hasproperty(result.data, col)
                @test names(result.data)[i] == string(col)
            end
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "measurements_logging")

        result = eegfun.erp_measurements(
            "erps_cleaned",
            "mean_amplitude",
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

    @testset "Known value validation" begin
        # Create a temporary directory for known-value tests
        known_test_dir = mktempdir()
        
        # Helper function to create ERP with known signal
        function create_known_signal_erp(participant::Int, condition::Int, signal_func::Function; 
                                         fs::Int = 1000, t_start::Float64 = 0.0, t_end::Float64 = 1.0)
            time = collect(range(t_start, t_end, length = Int((t_end - t_start) * fs) + 1))
            df = DataFrame(
                time = time,
                sample = 1:length(time),
                condition = fill(condition, length(time)),
                condition_name = fill("condition_$condition", length(time)),
                participant = fill(participant, length(time)),
            )
            
            # Create signal using the provided function
            df[!, :Ch1] = signal_func.(time)
            
            layout = eegfun.Layout(
                DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]),
                nothing,
                nothing,
            )
            
            analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)
            return eegfun.ErpData("test", condition, "condition_$condition", df, layout, fs, analysis_info, 1)
        end
        
        @testset "Mean amplitude" begin
            # Constant signal: mean should equal the constant
            constant_value = 5.0
            erp = create_known_signal_erp(1, 1, t -> constant_value, t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_mean_test.jld2")
            jldsave(file_path; data = [erp])
            
            result = eegfun.erp_measurements(
                "mean_test",
                "mean_amplitude",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            @test nrow(result.data) == 1
            @test isapprox(result.data[1, :Ch1], constant_value, rtol = 1e-6)
        end
        
        @testset "Peak amplitude and latency" begin
            # Signal with known peak: y = -(t - 0.5)^2 + 1
            # Peak at t=0.5 with value 1.0
            peak_time = 0.5
            peak_value = 1.0
            erp = create_known_signal_erp(1, 1, t -> -(t - peak_time)^2 + peak_value, t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_peak_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test max_peak_amplitude
            result = eegfun.erp_measurements(
                "peak_test",
                "max_peak_amplitude",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_peak"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            @test nrow(result.data) == 1
            # Should be close to peak value (allowing for sampling/discretization)
            @test isapprox(result.data[1, :Ch1], peak_value, rtol = 0.01)
            
            # Test max_peak_latency
            result_lat = eegfun.erp_measurements(
                "peak_test",
                "max_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_peak_lat"),
            )
            
            @test result_lat isa eegfun.ErpMeasurementsResult
            @test nrow(result_lat.data) == 1
            # Latency should be close to peak_time
            @test isapprox(result_lat.data[1, :Ch1], peak_time, rtol = 0.01)
            
            # Test min_peak_amplitude with inverted signal
            min_peak_time = 0.5
            min_peak_value = -1.0
            erp_min = create_known_signal_erp(1, 1, t -> (t - min_peak_time)^2 - 1.0, t_start = 0.0, t_end = 1.0)
            
            file_path_min = joinpath(known_test_dir, "1_min_peak_test.jld2")
            jldsave(file_path_min; data = [erp_min])
            
            result_min = eegfun.erp_measurements(
                "min_peak_test",
                "min_peak_amplitude",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_min"),
            )
            
            @test result_min isa eegfun.ErpMeasurementsResult
            @test isapprox(result_min.data[1, :Ch1], min_peak_value, rtol = 0.01)
            
            # Test min_peak_latency
            result_min_lat = eegfun.erp_measurements(
                "min_peak_test",
                "min_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_min_lat"),
            )
            
            @test result_min_lat isa eegfun.ErpMeasurementsResult
            @test isapprox(result_min_lat.data[1, :Ch1], min_peak_time, rtol = 0.01)
        end
        
        @testset "Peak-to-peak measurements" begin
            # Signal with known max and min peaks
            # y = -(t - 0.5)^2 + 1 for t < 0.5, y = (t - 0.5)^2 - 1 for t >= 0.5
            # Max peak at t=0.5 with value 1.0
            # Min peak at t=0.0 and t=1.0 with value 0.0 (but let's make it clearer)
            # Actually, let's use: y = sin(2πt) from 0 to 1, which has max at 0.25 and min at 0.75
            max_peak_time = 0.25
            min_peak_time = 0.75
            max_peak_value = 1.0
            min_peak_value = -1.0
            
            erp = create_known_signal_erp(1, 1, t -> sin(2π * t), t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_peak_to_peak_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test peak_to_peak_amplitude
            result_amp = eegfun.erp_measurements(
                "peak_to_peak_test",
                "peak_to_peak_amplitude",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_peak_to_peak_amp"),
            )
            
            @test result_amp isa eegfun.ErpMeasurementsResult
            @test nrow(result_amp.data) == 1
            # Peak-to-peak should be max - min = 1.0 - (-1.0) = 2.0
            @test isapprox(result_amp.data[1, :Ch1], max_peak_value - min_peak_value, rtol = 0.05)
            
            # Test peak_to_peak_latency
            result_lat = eegfun.erp_measurements(
                "peak_to_peak_test",
                "peak_to_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_peak_to_peak_lat"),
            )
            
            @test result_lat isa eegfun.ErpMeasurementsResult
            @test nrow(result_lat.data) == 1
            # Peak-to-peak latency should be |max_time - min_time| = |0.25 - 0.75| = 0.5
            expected_lat_diff = abs(max_peak_time - min_peak_time)
            @test isapprox(result_lat.data[1, :Ch1], expected_lat_diff, rtol = 0.05)
        end
        
        @testset "Area measurements" begin
            # Constant signal: area = value * duration
            constant_value = 3.0
            t_start = 0.0
            t_end = 1.0
            expected_area = constant_value * (t_end - t_start)
            
            erp = create_known_signal_erp(1, 1, t -> constant_value, t_start = t_start, t_end = t_end)
            
            file_path = joinpath(known_test_dir, "1_area_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test rectified_area (should be same as integral for positive constant)
            result = eegfun.erp_measurements(
                "area_test",
                "rectified_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_area"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            @test isapprox(result.data[1, :Ch1], expected_area, rtol = 0.01)
            
            # Test integral (same for constant)
            result_int = eegfun.erp_measurements(
                "area_test",
                "integral",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_integral"),
            )
            
            @test isapprox(result_int.data[1, :Ch1], expected_area, rtol = 0.01)
            
            # Test positive_area (same for positive constant)
            result_pos = eegfun.erp_measurements(
                "area_test",
                "positive_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_pos"),
            )
            
            @test isapprox(result_pos.data[1, :Ch1], expected_area, rtol = 0.01)
            
            # Test negative_area (should be 0 for positive constant)
            result_neg = eegfun.erp_measurements(
                "area_test",
                "negative_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_neg"),
            )
            
            @test isapprox(result_neg.data[1, :Ch1], 0.0, atol = 0.01)
            
            # Test with negative constant
            neg_value = -2.0
            expected_neg_area = abs(neg_value) * (t_end - t_start)
            erp_neg = create_known_signal_erp(1, 1, t -> neg_value, t_start = t_start, t_end = t_end)
            
            file_path_neg = joinpath(known_test_dir, "1_area_neg_test.jld2")
            jldsave(file_path_neg; data = [erp_neg])
            
            result_neg2 = eegfun.erp_measurements(
                "area_neg_test",
                "negative_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_neg2"),
            )
            
            @test isapprox(result_neg2.data[1, :Ch1], expected_neg_area, rtol = 0.01)
            
            # Test rectified_area with negative constant (should be same as negative_area)
            result_rect_neg = eegfun.erp_measurements(
                "area_neg_test",
                "rectified_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_rect_neg"),
            )
            
            @test isapprox(result_rect_neg.data[1, :Ch1], expected_neg_area, rtol = 0.01)
            
            # Test integral with negative constant (should be negative)
            result_int_neg = eegfun.erp_measurements(
                "area_neg_test",
                "integral",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_int_neg"),
            )
            
            expected_integral_neg = neg_value * (t_end - t_start)
            @test isapprox(result_int_neg.data[1, :Ch1], expected_integral_neg, rtol = 0.01)
            
            # Test positive_area with negative constant (should be 0)
            result_pos_neg = eegfun.erp_measurements(
                "area_neg_test",
                "positive_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_pos_neg"),
            )
            
            @test isapprox(result_pos_neg.data[1, :Ch1], 0.0, atol = 0.01)
        end
        
        @testset "Mixed positive/negative area measurements" begin
            # Signal that crosses zero: y = t - 0.5, from 0 to 1
            # Positive from 0.5 to 1.0, negative from 0.0 to 0.5
            # Total area: ∫₀¹ (t - 0.5) dt = [t²/2 - 0.5t]₀¹ = 0.5 - 0.5 - 0 = 0
            # Positive area: ∫₀.₅¹ (t - 0.5) dt = [t²/2 - 0.5t]₀.₅¹ = 0.5 - 0.5 - (0.125 - 0.25) = 0.125
            # Negative area: ∫₀⁰.⁵ (t - 0.5) dt = [t²/2 - 0.5t]₀⁰.⁵ = 0.125 - 0.25 - 0 = -0.125 (abs = 0.125)
            # Rectified area: 0.125 + 0.125 = 0.25
            t_start = 0.0
            t_end = 1.0
            zero_crossing = 0.5
            
            erp = create_known_signal_erp(1, 1, t -> t - zero_crossing, t_start = t_start, t_end = t_end)
            
            file_path = joinpath(known_test_dir, "1_mixed_area_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Expected values
            # Signal: y = t - 0.5 from 0 to 1
            # Integral: ∫₀¹ (t - 0.5) dt = [t²/2 - 0.5t]₀¹ = 0.5 - 0.5 - 0 = 0
            # Positive area: ∫₀.₅¹ (t - 0.5) dt = [t²/2 - 0.5t]₀.₅¹ = (0.5 - 0.5) - (0.125 - 0.25) = 0.125
            # Negative area: ∫₀⁰.⁵ |t - 0.5| dt = ∫₀⁰.⁵ (0.5 - t) dt = [0.5t - t²/2]₀⁰.⁵ = 0.25 - 0.125 = 0.125
            # Rectified area: 0.125 + 0.125 = 0.25
            expected_integral = 0.0
            expected_positive_area = 0.5 * (t_end - zero_crossing)^2  # = 0.5 * 0.5² = 0.125
            expected_negative_area = 0.5 * (zero_crossing - t_start)^2  # = 0.5 * 0.5² = 0.125
            expected_rectified_area = expected_positive_area + expected_negative_area  # = 0.25
            
            # Test integral
            result_int = eegfun.erp_measurements(
                "mixed_area_test",
                "integral",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_mixed_int"),
            )
            
            @test isapprox(result_int.data[1, :Ch1], expected_integral, atol = 0.01)
            
            # Test positive_area
            result_pos = eegfun.erp_measurements(
                "mixed_area_test",
                "positive_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_mixed_pos"),
            )
            
            @test isapprox(result_pos.data[1, :Ch1], expected_positive_area, rtol = 0.01)
            
            # Test negative_area
            result_neg = eegfun.erp_measurements(
                "mixed_area_test",
                "negative_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_mixed_neg"),
            )
            
            @test isapprox(result_neg.data[1, :Ch1], expected_negative_area, rtol = 0.01)
            
            # Test rectified_area
            result_rect = eegfun.erp_measurements(
                "mixed_area_test",
                "rectified_area",
                analysis_window = eegfun.samples((t_start, t_end)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_mixed_rect"),
            )
            
            @test isapprox(result_rect.data[1, :Ch1], expected_rectified_area, rtol = 0.01)
        end
        
        @testset "Robust peak detection" begin
            # Signal with multiple local peaks but one clear global peak
            # y = -2*(t - 0.5)^2 + 1 + 0.1*sin(20πt)
            # Main peak at t=0.5, but with small oscillations
            peak_time = 0.5
            peak_value = 1.0
            
            erp = create_known_signal_erp(1, 1, 
                t -> -2*(t - peak_time)^2 + peak_value + 0.1*sin(20π*t),
                t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_robust_peak_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test with small local_window (should find robust peak)
            result = eegfun.erp_measurements(
                "robust_peak_test",
                "max_peak_amplitude",
                analysis_window = eegfun.samples((0.0, 1.0)),
                local_window = 5,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_robust"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            # Should find the main peak, not the small oscillations
            @test result.data[1, :Ch1] > 0.9  # Should be close to 1.0
            
            # Test latency
            result_lat = eegfun.erp_measurements(
                "robust_peak_test",
                "max_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                local_window = 5,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_robust_lat"),
            )
            
            @test isapprox(result_lat.data[1, :Ch1], peak_time, rtol = 0.05)
        end
        
        @testset "Fractional area latency edge cases" begin
            # Test with different fractions
            t_start = 0.0
            t_end = 1.0
            
            erp = create_known_signal_erp(1, 1, t -> t, t_start = t_start, t_end = t_end)
            
            file_path = joinpath(known_test_dir, "1_frac_area_edge_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test fraction = 0.0 (should return start time)
            result_0 = eegfun.erp_measurements(
                "frac_area_edge_test",
                "fractional_area_latency",
                analysis_window = eegfun.samples((t_start, t_end)),
                fractional_area_fraction = 0.0,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_0"),
            )
            
            @test isapprox(result_0.data[1, :Ch1], t_start, rtol = 0.01)
            
            # Test fraction = 1.0 (should return end time)
            result_1 = eegfun.erp_measurements(
                "frac_area_edge_test",
                "fractional_area_latency",
                analysis_window = eegfun.samples((t_start, t_end)),
                fractional_area_fraction = 1.0,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_1"),
            )
            
            @test isapprox(result_1.data[1, :Ch1], t_end, rtol = 0.01)
            
            # Test fraction = 0.25
            # For linear signal y = t: cumulative area = t²/2
            # Total area = 0.5
            # For fraction 0.25: cumulative area = 0.5 * 0.25 = 0.125
            # So t²/2 = 0.125, which means t² = 0.25, so t = 0.5
            fraction = 0.25
            expected_latency = sqrt(2 * 0.5 * fraction)  # t²/2 = total_area * fraction, so t = √(2 * total_area * fraction)
            # For fraction 0.25: t = √(2 * 0.5 * 0.25) = √0.25 = 0.5
            result_25 = eegfun.erp_measurements(
                "frac_area_edge_test",
                "fractional_area_latency",
                analysis_window = eegfun.samples((t_start, t_end)),
                fractional_area_fraction = fraction,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_25"),
            )
            
            @test isapprox(result_25.data[1, :Ch1], expected_latency, rtol = 0.05)
        end
        
        @testset "Fractional peak latency edge cases" begin
            # Triangular signal with known peak
            peak_time = 0.5
            peak_value = 1.0
            
            erp = create_known_signal_erp(1, 1, 
                t -> t < peak_time ? (peak_value / peak_time) * t : peak_value - (peak_value / (1.0 - peak_time)) * (t - peak_time),
                t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_frac_peak_edge_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test fraction = 0.0 (should return start/end depending on direction)
            result_onset_0 = eegfun.erp_measurements(
                "frac_peak_edge_test",
                "fractional_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                fractional_peak_fraction = 0.0,
                fractional_peak_direction = :onset,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_peak_onset_0"),
            )
            
            @test isapprox(result_onset_0.data[1, :Ch1], 0.0, rtol = 0.01)
            
            result_offset_0 = eegfun.erp_measurements(
                "frac_peak_edge_test",
                "fractional_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                fractional_peak_fraction = 0.0,
                fractional_peak_direction = :offset,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_peak_offset_0"),
            )
            
            @test isapprox(result_offset_0.data[1, :Ch1], 1.0, rtol = 0.01)
            
            # Test fraction = 1.0 (should return peak time)
            result_onset_1 = eegfun.erp_measurements(
                "frac_peak_edge_test",
                "fractional_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                fractional_peak_fraction = 1.0,
                fractional_peak_direction = :onset,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_peak_onset_1"),
            )
            
            @test isapprox(result_onset_1.data[1, :Ch1], peak_time, rtol = 0.01)
        end
        
        @testset "Fractional area latency" begin
            # Linear signal: y = t, from 0 to 1
            # Cumulative area from 0 to t: ∫₀ᵗ t dt = t²/2
            # Total area: ∫₀¹ t dt = 1/2 = 0.5
            # For fraction 0.5, we want cumulative area = 0.5 * 0.5 = 0.25
            # So t²/2 = 0.25, which means t² = 0.5, so t = √0.5 ≈ 0.707
            t_start = 0.0
            t_end = 1.0
            total_area = 0.5 * (t_end - t_start)^2  # Integral of t from 0 to 1 = 0.5
            fraction = 0.5
            # Solve: t²/2 = total_area * fraction = 0.5 * 0.5 = 0.25
            # So t² = 0.5, t = √0.5
            expected_latency = sqrt(0.5)
            
            erp = create_known_signal_erp(1, 1, t -> t, t_start = t_start, t_end = t_end)
            
            file_path = joinpath(known_test_dir, "1_frac_area_test.jld2")
            jldsave(file_path; data = [erp])
            
            result = eegfun.erp_measurements(
                "frac_area_test",
                "fractional_area_latency",
                analysis_window = eegfun.samples((t_start, t_end)),
                fractional_area_fraction = fraction,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_area"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            @test isapprox(result.data[1, :Ch1], expected_latency, rtol = 0.05)  # Allow some tolerance for discretization
        end
        
        @testset "Fractional peak latency" begin
            # Signal with known peak: y = -(t - 0.5)^2 + 1, peak at t=0.5, value=1.0
            # For fraction 0.5, we want value = 0.5
            # Solve: -(t - 0.5)^2 + 1 = 0.5
            # (t - 0.5)^2 = 0.5
            # t = 0.5 ± sqrt(0.5) ≈ 0.5 ± 0.707
            # For onset (before peak): t ≈ 0.5 - 0.707 = -0.207 (clamped to start)
            # For offset (after peak): t ≈ 0.5 + 0.707 = 1.207 (clamped to end)
            # Actually, let's use a simpler approach with a linear ramp to peak
            peak_time = 0.5
            peak_value = 1.0
            fraction = 0.5
            target_value = peak_value * fraction
            
            # Create signal: linear from 0 to peak, then linear down
            erp = create_known_signal_erp(1, 1, 
                t -> t < peak_time ? (peak_value / peak_time) * t : peak_value - (peak_value / (1.0 - peak_time)) * (t - peak_time),
                t_start = 0.0, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_frac_peak_test.jld2")
            jldsave(file_path; data = [erp])
            
            # Test onset direction
            expected_onset = peak_time * fraction  # Linear: value = (peak/peak_time) * t, so t = value * peak_time / peak
            result_onset = eegfun.erp_measurements(
                "frac_peak_test",
                "fractional_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                fractional_peak_fraction = fraction,
                fractional_peak_direction = :onset,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_peak_onset"),
            )
            
            @test result_onset isa eegfun.ErpMeasurementsResult
            @test isapprox(result_onset.data[1, :Ch1], expected_onset, rtol = 0.05)
            
            # Test offset direction
            # After peak: value = peak - (peak/(1-peak_time)) * (t - peak_time)
            # For target_value: target = peak - (peak/(1-peak_time)) * (t - peak_time)
            # t - peak_time = (peak - target) * (1 - peak_time) / peak
            expected_offset = peak_time + (peak_value - target_value) * (1.0 - peak_time) / peak_value
            result_offset = eegfun.erp_measurements(
                "frac_peak_test",
                "fractional_peak_latency",
                analysis_window = eegfun.samples((0.0, 1.0)),
                fractional_peak_fraction = fraction,
                fractional_peak_direction = :offset,
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_frac_peak_offset"),
            )
            
            @test result_offset isa eegfun.ErpMeasurementsResult
            @test isapprox(result_offset.data[1, :Ch1], expected_offset, rtol = 0.05)
        end
        
        @testset "Baseline correction" begin
            # Signal: constant 5.0, baseline from -0.2 to 0.0 should be 5.0
            # After baseline correction, signal should be ~0.0
            constant_value = 5.0
            erp = create_known_signal_erp(1, 1, t -> constant_value, t_start = -0.2, t_end = 1.0)
            
            file_path = joinpath(known_test_dir, "1_baseline_test.jld2")
            jldsave(file_path; data = [erp])
            
            result = eegfun.erp_measurements(
                "baseline_test",
                "mean_amplitude",
                analysis_window = eegfun.samples((0.1, 0.2)),
                baseline_window = eegfun.samples((-0.2, 0.0)),
                input_dir = known_test_dir,
                output_dir = joinpath(known_test_dir, "output_baseline"),
            )
            
            @test result isa eegfun.ErpMeasurementsResult
            # After baseline correction, mean should be close to 0
            @test isapprox(result.data[1, :Ch1], 0.0, atol = 0.1)
        end
        
        # Cleanup
        rm(known_test_dir, recursive = true)
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
