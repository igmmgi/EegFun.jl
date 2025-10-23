using Test
using DataFrames
using eegfun
using JLD2
using Statistics



@testset "rereference" begin

    dat = create_test_data(n_channels = 2)
    # 1) Average reference (:avg) subtracts mean of all channels from each selected channel
    eegfun.rereference!(dat, :avg)
    # After average reference, the mean of all channels should be approximately zero
    mean_of_channels = mean([dat.data.Ch1, dat.data.Ch2])
    @test all(abs.(mean_of_channels) .< 1e-10)
    @test eegfun.reference(dat.analysis_info) == :avg

    # 2) Single channel reference uses Ch1 as reference
    dat = create_test_data(n_channels = 2)
    dat_ref = eegfun.rereference(dat, :Ch1)
    # After single channel reference, the reference channel should be zero
    @test all(dat_ref.data.Ch1 .== 0.0)
    # Other channels should be original minus reference
    @test all(dat_ref.data.Ch2 .== dat.data.Ch2 .- dat.data.Ch1)
    @test eegfun.reference(dat_ref.analysis_info) == :Ch1

    # 4) EpochData: per-epoch reference computed independently
    dat = create_test_epoch_data(n_epochs = 2, n_channels = 2, conditions = 1)
    eegfun.rereference!(dat, :avg)
    # After average reference, the mean of all channels should be approximately zero for each epoch
    for epoch in dat.data
        mean_of_channels = mean([epoch.Ch1, epoch.Ch2])
        @test all(abs.(mean_of_channels) .< 1e-10)
    end

    # 5) Non-mutating version returns a new object; original unchanged
    dat = create_test_epoch_data(n_epochs = 2, n_channels = 2, conditions = 1)
    dat_ref = eegfun.rereference(dat, :Ch1)
    @test :Ch1 ∈ propertynames(dat_ref.data[1]) && :Ch1 ∈ propertynames(dat.data[1])
    @test !all(dat_ref.data[1].Ch1 .== dat.data[1].Ch1)


    # 7) Channel included in both reference and selection becomes zero
    dat = create_test_epoch_data(n_epochs = 2, n_channels = 2, conditions = 1)
    eegfun.rereference!(dat, [:Ch1], eegfun.channels([:Ch1]))
    @test all(dat.data[1].Ch1 .== 0.0)

    # 9) Missing reference channel should throw
    dat = create_test_epoch_data(n_epochs = 2, n_channels = 2, conditions = 1)
    @test_throws Any eegfun.rereference!(copy(dat), [:Z], eegfun.channels([:A]))

end



@testset "Batch Rereference" begin
    # Create temporary test directory
    test_dir = mktempdir()

    @testset "Basic rereferencing" begin
        # Create test ERP files
        for participant = 1:3
            erps = [create_test_erp_data(participant, 1), create_test_erp_data(participant, 2)]

            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            save(file_path, "erps", erps)
        end

        output_dir = joinpath(test_dir, "rereferenced")

        # Test average reference
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 4  # One for each participant

        # Load and verify one output file
        output_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        @test isfile(output_file)

        rereferenced_erps = load(output_file, "erps")
        @test length(rereferenced_erps) == 2  # Two conditions

        # Verify data structure is preserved
        for erp in rereferenced_erps
            @test erp isa eegfun.ErpData
            @test nrow(erp.data) == 2501  # n_timepoints (-0.5 to 2.0s at 1000Hz)
            @test "time" in names(erp.data)
            @test "condition" in names(erp.data)
            @test "Ch1" in names(erp.data)
            @test "Ch2" in names(erp.data)
            @test "Ch3" in names(erp.data)
        end
    end

    @testset "Different reference types" begin
        @testset "Average reference" begin
            output_dir = joinpath(test_dir, "rereferenced_avg")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :avg,
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end

        @testset "Mastoid reference" begin
            output_dir = joinpath(test_dir, "rereferenced_mastoid")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :Ch2,
                output_dir = output_dir,
            )

            # Ch2 reference should work
            @test result !== nothing
            @test result.success > 0  # Should have successful files
            @test result.errors == 0  # Should have no errors
        end

        @testset "Single channel reference" begin
            output_dir = joinpath(test_dir, "rereferenced_cz")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:Ch1],
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end

        @testset "Multiple channel reference" begin
            output_dir = joinpath(test_dir, "rereferenced_multiple")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:Ch1, :Ch2],
                output_dir = output_dir,
            )

            @test isdir(output_dir)
            @test length(readdir(output_dir)) == 4
        end
    end

    @testset "Epoch data processing" begin
        # Create test epoch files
        for participant = 1:2
            epochs = create_test_epoch_data(conditions = 2, n_channels = 3)  # This returns Vector{EpochData}

            file_path = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
            save(file_path, "epochs", epochs)
        end

        output_dir = joinpath(test_dir, "rereferenced_epochs")

        result = eegfun.rereference(
            "epochs_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) >= 2  # At least one for each participant

        # Load and verify epoch data structure
        output_file = joinpath(output_dir, "1_epochs_cleaned.jld2")
        rereferenced_epochs = load(output_file, "epochs")
        @test length(rereferenced_epochs) == 2  # Two conditions

        for epoch_data in rereferenced_epochs
            @test epoch_data isa eegfun.EpochData
            @test length(epoch_data.data) == 10  # 10 epochs per condition (default)
        end
    end

    @testset "Participant and condition filtering" begin
        output_dir = joinpath(test_dir, "rereferenced_filtered")

        # Test participant filtering
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = [1, 2],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 3  # Only participants 1 and 2, plus empty file

        # Test condition filtering
        output_dir2 = joinpath(test_dir, "rereferenced_conditions")
        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir2,
        )

        @test isdir(output_dir2)
        output_files2 = readdir(output_dir2)
        @test length(output_files2) == 4  # All participants, but only condition 1

        # Verify only condition 1 was processed
        output_file = joinpath(output_dir2, "1_erps_cleaned.jld2")
        rereferenced_erps = load(output_file, "erps")
        @test length(rereferenced_erps) == 1  # Only condition 1
        @test rereferenced_erps[1].data.condition[1] == 1
    end

    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.rereference("erps_cleaned", input_dir = "/nonexistent/dir")
        end

        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "rereferenced_none")

            result = eegfun.rereference("nonexistent_pattern", input_dir = test_dir, output_dir = output_dir)

            @test result === nothing
        end

        @testset "Files with no recognized data variable" begin
            # Create file with unrecognized variable and non-matching participant
            unrecognized_file = joinpath(test_dir, "999_unrecognized_erps_cleaned.jld2")
            save(unrecognized_file, "other_data", "test")

            output_dir = joinpath(test_dir, "rereferenced_unrecognized")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                participants = 999,  # This should match the file but fail processing
                output_dir = output_dir,
            )

            # Should return BatchResult with errors, not nothing
            @test result !== nothing
            @test result.success == 0
            @test result.errors > 0
        end
    end

    @testset "Edge cases" begin
        @testset "Empty data files" begin
            # Create file with empty data
            empty_file = joinpath(test_dir, "empty_erps_cleaned.jld2")
            save(empty_file, "erps", eegfun.ErpData[])

            output_dir = joinpath(test_dir, "rereferenced_empty")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                participants = 999,  # Non-existent participant
                output_dir = output_dir,
            )

            # Should return BatchResult since both empty_erps_cleaned.jld2 and 999_unrecognized_erps_cleaned.jld2 have no participant number
            @test result !== nothing
            @test result.success == 1  # Only the empty file gets processed successfully
            @test result.errors == 1   # One error for the unrecognized file
        end

        @testset "Invalid reference channels" begin
            output_dir = joinpath(test_dir, "rereferenced_invalid")

            # This should handle invalid channels gracefully in batch processing
            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = [:InvalidChannel],
                output_dir = output_dir,
            )

            # Should have errors due to invalid channels
            @test result !== nothing
            @test result.errors > 0
        end
    end

    @testset "Data integrity verification" begin
        output_dir = joinpath(test_dir, "rereferenced_integrity")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = 1,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Load original and rereferenced data
        original_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        original_erps = load(original_file, "erps")
        original_erp = original_erps[1]  # Condition 1

        rereferenced_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        rereferenced_erps = load(rereferenced_file, "erps")
        rereferenced_erp = rereferenced_erps[1]  # Condition 1

        # Verify metadata is preserved
        @test rereferenced_erp.sample_rate == original_erp.sample_rate
        @test rereferenced_erp.n_epochs == original_erp.n_epochs

        # Verify layout is preserved (check channel labels)
        @test rereferenced_erp.layout.data.label == original_erp.layout.data.label

        # Verify analysis info reflects rereferencing
        @test rereferenced_erp.analysis_info.reference == :avg
        @test original_erp.analysis_info.reference == :none

        # Verify time column is unchanged
        @test all(rereferenced_erp.data.time .== original_erp.data.time)

        # Verify condition information is preserved
        @test all(rereferenced_erp.data.condition .== original_erp.data.condition)
        @test all(rereferenced_erp.data.condition_name .== original_erp.data.condition_name)
        @test all(rereferenced_erp.data.participant .== original_erp.data.participant)
    end

    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_rereferenced")

            result = eegfun.rereference(
                "erps_cleaned",
                input_dir = test_dir,
                reference_selection = :avg,
                output_dir = custom_dir,
            )

            @test isdir(custom_dir)
            @test length(readdir(custom_dir)) == 5
        end

        @testset "Auto-generated output directory" begin
            result = eegfun.rereference("erps_cleaned", input_dir = test_dir, reference_selection = :avg)

            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "rereferenced_erps_cleaned_avg")
            @test isdir(expected_dir)
        end
    end

    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "rereferenced_logging")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Check that log file was created
        log_file = joinpath(output_dir, "rereference.log")
        @test isfile(log_file)

        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch rereferencing started", log_content)
        @test occursin("Found 5 JLD2 files", log_content)
        @test occursin("Reference settings: avg", log_content)
        @test occursin("Batch operation complete", log_content)
    end

    @testset "Reference calculation verification" begin
        output_dir = joinpath(test_dir, "rereferenced_verification")

        result = eegfun.rereference(
            "erps_cleaned",
            input_dir = test_dir,
            participants = 1,
            conditions = [1],
            reference_selection = :avg,
            output_dir = output_dir,
        )

        # Load original and rereferenced data
        original_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        original_erps = load(original_file, "erps")
        original_erp = original_erps[1]

        rereferenced_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        rereferenced_erps = load(rereferenced_file, "erps")
        rereferenced_erp = rereferenced_erps[1]

        # Verify average reference calculation
        # Average reference should make the mean of all channels zero
        for ch in [:Ch1, :Ch2, :Ch3]
            if hasproperty(rereferenced_erp.data, ch)
                # The mean across all channels at each time point should be approximately zero
                all_channels =
                    [rereferenced_erp.data[!, :Ch1], rereferenced_erp.data[!, :Ch2], rereferenced_erp.data[!, :Ch3]]
                mean_across_channels = mean.(zip(all_channels...))
                @test all(abs.(mean_across_channels) .< 1e-10)
            end
        end
    end

    @testset "Different file patterns" begin
        # Test with different file patterns
        for pattern in ["erps_cleaned", "nonexistent", "original"]
            output_dir = joinpath(test_dir, "rereferenced_$pattern")

            result =
                eegfun.rereference(pattern, input_dir = test_dir, reference_selection = :avg, output_dir = output_dir)

            if pattern in ["erps_cleaned"]
                @test isdir(output_dir)
                @test length(readdir(output_dir)) > 0
            else
                # These patterns shouldn't match any files
                @test result === nothing
            end
        end
    end

    # Cleanup
    rm(test_dir, recursive = true)
end
