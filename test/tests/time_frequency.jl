using Test
using JLD2
using DataFrames

@testset "Time-Frequency Batch Processing" begin

    @testset "Batch tf_morlet" begin
        mktempdir() do tmpdir
            # Create test epoch files (2 participants, 2 conditions each)
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_morlet_output")

            # Run batch tf_morlet
            EegFun.tf_morlet("epochs"; input_dir = tmpdir, output_dir = output_dir, cycles = 7, frequencies = range(2, 20, length = 10))

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2  # 2 conditions
        end
    end

    @testset "Batch tf_stft" begin
        mktempdir() do tmpdir
            # Create test epoch files
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_stft_output")

            # Run batch tf_stft with fixed window
            EegFun.tf_stft(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                window_length = 0.3,
                frequencies = range(2, 20, length = 10),
            )

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2
        end
    end

    @testset "Batch tf_multitaper" begin
        mktempdir() do tmpdir
            # Create test epoch files
            for i = 1:2
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_multitaper_output")

            # Run batch tf_multitaper
            EegFun.tf_multitaper("epochs"; input_dir = tmpdir, output_dir = output_dir, cycles = 5, frequencies = range(2, 20, length = 10))

            # Check output files exist
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs.jld2"))

            # Load and verify structure
            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test result isa Vector{EegFun.TimeFreqData}
            @test length(result) == 2
        end
    end

    @testset "Participant filtering" begin
        mktempdir() do tmpdir
            # Create files for 3 participants
            for i = 1:3
                epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
                jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epochs)
            end

            output_dir = joinpath(tmpdir, "tf_filtered")

            # Process only participants 1 and 3
            EegFun.tf_morlet(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                participant_selection = EegFun.participants([1, 3]),
                cycles = 7,
                frequencies = range(2, 20, length = 10),
            )

            # Only participants 1 and 3 should be processed
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs.jld2"))
            @test isfile(joinpath(output_dir, "3_epochs.jld2"))
        end
    end

    @testset "Condition selection" begin
        mktempdir() do tmpdir
            # Create file with 3 conditions
            epochs = EegFun.create_test_epoch_data_vector(conditions = 1:3, n_epochs = 5, n_channels = 3, fs = 256, n = 256)
            jldsave(joinpath(tmpdir, "1_epochs.jld2"); data = epochs)

            output_dir = joinpath(tmpdir, "tf_cond_select")

            # Select only conditions 1 and 3
            EegFun.tf_morlet(
                "epochs";
                input_dir = tmpdir,
                output_dir = output_dir,
                condition_selection = EegFun.conditions([1, 3]),
                cycles = 7,
                frequencies = range(2, 20, length = 10),
            )

            result = load(joinpath(output_dir, "1_epochs.jld2"), "data")
            @test length(result) == 2  # Only 2 conditions
        end
    end

    @testset "No matching files" begin
        mktempdir() do tmpdir
            # No files in directory
            result = EegFun.tf_morlet("nonexistent"; input_dir = tmpdir, cycles = 7)
            @test isnothing(result)
        end
    end

end
