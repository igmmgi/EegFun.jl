"""
Test resampling functionality for continuous, epoched, and ERP data.
"""

using Test
using DataFrames
using JLD2

@testset "Resample tests" begin

    #=============================================================================
        CONTINUOUS DATA RESAMPLING
    =============================================================================#

    @testset "Continuous data resampling" begin

        @testset "Basic downsampling by factor 2" begin
            # Create continuous data at 512 Hz
            n_samples = 1000
            sample_rate = 512
            time = collect(0:(n_samples-1)) ./ sample_rate

            data = DataFrame(
                time = time,
                C3 = sin.(2π .* 10 .* time),  # 10 Hz signal
                C4 = cos.(2π .* 10 .* time),
                triggers = zeros(Int, n_samples),
            )
            data.triggers[100] = 1  # Add a trigger
            data.triggers[500] = 2

            continuous = eegfun.ContinuousData(
                "test_data",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                sample_rate,
                eegfun.AnalysisInfo(),
            )

            # Resample by factor of 2
            resampled = eegfun.resample(continuous, 2)

            # Check sample rate
            @test resampled.sample_rate == 256
            @test continuous.sample_rate == 512  # Original unchanged

            # Check number of samples (exactly 500 from regular downsampling)
            @test nrow(resampled.data) == 500

            # Check that we kept the right samples (every 2nd)
            @test resampled.data.time[1] == continuous.data.time[1]
            @test resampled.data.time[2] == continuous.data.time[3]
            @test resampled.data.time[end] == continuous.data.time[end-1]

            # Check that triggers are preserved
            @test sum(resampled.data.triggers .== 1) == 1
            @test sum(resampled.data.triggers .== 2) == 1

            # Check columns are preserved
            @test names(resampled.data) == names(continuous.data)
        end

        @testset "Mutating version" begin
            # Create continuous data
            n_samples = 1000
            sample_rate = 512
            time = collect(0:(n_samples-1)) ./ sample_rate

            data = DataFrame(time = time, C3 = randn(n_samples), C4 = randn(n_samples))

            continuous = eegfun.ContinuousData(
                "test_data",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                sample_rate,
                eegfun.AnalysisInfo(),
            )

            original_nrow = nrow(continuous.data)

            # Resample in-place
            eegfun.resample!(continuous, 2)

            # Check modifications
            @test continuous.sample_rate == 256
            @test nrow(continuous.data) == original_nrow ÷ 2
        end

        @testset "Factor of 4 downsampling" begin
            # Create data at 1024 Hz
            n_samples = 1024
            sample_rate = 1024
            time = collect(0:(n_samples-1)) ./ sample_rate

            data = DataFrame(time = time, C3 = randn(n_samples))

            continuous = eegfun.ContinuousData(
                "test_data",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                sample_rate,
                eegfun.AnalysisInfo(),
            )

            # Downsample by 4
            resampled = eegfun.resample(continuous, 4)

            @test resampled.sample_rate == 256
            @test nrow(resampled.data) == 256
            @test resampled.data.time[1] == 0.0
            @test resampled.data.time[2] ≈ 4 / 1024
        end

        @testset "Factor of 1 (no change)" begin
            continuous = eegfun.ContinuousData(
                "test_data",
                DataFrame(time = [0.0, 0.001], C3 = [1.0, 2.0]),
                eegfun.Layout(DataFrame(), nothing, nothing),
                1000,
                eegfun.AnalysisInfo(),
            )

            original_nrow = nrow(continuous.data)
            eegfun.resample!(continuous, 1)

            @test continuous.sample_rate == 1000
            @test nrow(continuous.data) == original_nrow
        end

        @testset "Metadata preservation" begin
            # Create data with multiple metadata columns
            n_samples = 500
            data = DataFrame(
                time = collect(0:(n_samples-1)) ./ 500,
                C3 = randn(n_samples),
                C4 = randn(n_samples),
                triggers = zeros(Int, n_samples),
                trial = ones(Int, n_samples),
                condition = fill("A", n_samples),
                rt = fill(0.5, n_samples),
            )
            data.triggers[100] = 1

            continuous = eegfun.ContinuousData(
                "test_data",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                500,
                eegfun.AnalysisInfo(),
            )

            resampled = eegfun.resample(continuous, 5)

            # All columns should be present
            @test names(resampled.data) == names(continuous.data)

            # Metadata values should be preserved (for kept samples)
            @test all(resampled.data.trial .== 1)
            @test all(resampled.data.condition .== "A")
            @test all(resampled.data.rt .== 0.5)

            # Triggers should be preserved
            @test sum(resampled.data.triggers .== 1) == 1
        end
    end


    #=============================================================================
        EPOCHED DATA RESAMPLING
    =============================================================================#

    @testset "Epoched data resampling" begin

        @testset "Basic epoch downsampling" begin
            # Create 3 epochs at 512 Hz
            n_samples = 256
            n_epochs = 3
            sample_rate = 512

            epochs = []
            for i = 1:n_epochs
                time = collect(-0.2:(1/sample_rate):0.3)[1:n_samples]
                epoch = DataFrame(
                    time = time,
                    C3 = randn(n_samples),
                    C4 = randn(n_samples),
                    trial = fill(i, n_samples),
                    rt = fill(0.5 + 0.1*i, n_samples),
                )
                push!(epochs, epoch)
            end

            epoch_data = eegfun.EpochData(
                "test_data",
                1,
                "condition_1",
                epochs,
                eegfun.Layout(DataFrame(), nothing, nothing),
                sample_rate,
                eegfun.AnalysisInfo(),
            )

            # Resample by factor of 2
            resampled = eegfun.resample(epoch_data, 2)

            # Check sample rate
            @test resampled.sample_rate == 256
            @test epoch_data.sample_rate == 512  # Original unchanged

            # Check number of epochs
            @test length(resampled.data) == 3

            # Check samples per epoch
            @test nrow(resampled.data[1]) == 128
            @test nrow(resampled.data[2]) == 128
            @test nrow(resampled.data[3]) == 128

            # Check metadata preservation
            @test resampled.data[1].trial[1] == 1
            @test resampled.data[2].trial[1] == 2
            @test resampled.data[3].trial[1] == 3

            @test all(resampled.data[1].rt .≈ 0.6)
            @test all(resampled.data[2].rt .≈ 0.7)
            @test all(resampled.data[3].rt .≈ 0.8)
        end

        @testset "Mutating epoch resampling" begin
            # Create epochs
            epochs = []
            for i = 1:2
                epoch = DataFrame(time = collect(0:511) ./ 512, C3 = randn(512))
                push!(epochs, epoch)
            end

            epoch_data = eegfun.EpochData(
                "test_data",
                1,
                "condition_1",
                epochs,
                eegfun.Layout(DataFrame(), nothing, nothing),
                512,
                eegfun.AnalysisInfo(),
            )

            original_n_samples = nrow(epoch_data.data[1])

            # Resample in-place
            eegfun.resample!(epoch_data, 4)

            @test epoch_data.sample_rate == 128
            @test nrow(epoch_data.data[1]) == original_n_samples ÷ 4
            @test nrow(epoch_data.data[2]) == original_n_samples ÷ 4
        end

        @testset "Epoch column preservation" begin
            # Create epochs with triggers
            epochs = []
            for i = 1:2
                epoch = DataFrame(
                    time = collect(0:255) ./ 256,
                    C3 = randn(256),
                    trigger = zeros(Int, 256),
                    condition = fill("A", 256),
                )
                epoch.trigger[50] = i
                push!(epochs, epoch)
            end

            epoch_data = eegfun.EpochData(
                "test_data",
                1,
                "condition_1",
                epochs,
                eegfun.Layout(DataFrame(), nothing, nothing),
                256,
                eegfun.AnalysisInfo(),
            )

            resampled = eegfun.resample(epoch_data, 2)

            # Check columns preserved
            @test names(resampled.data[1]) == names(epoch_data.data[1])

            # Check trigger preserved (may not be at exact same index)
            @test sum(resampled.data[1].trigger .== 1) >= 1
            @test sum(resampled.data[2].trigger .== 2) >= 1
        end
    end


    #=============================================================================
        ERP DATA RESAMPLING
    =============================================================================#

    @testset "ERP data resampling" begin

        @testset "Basic ERP downsampling" begin
            # Create ERP at 512 Hz
            n_samples = 512
            sample_rate = 512
            time = collect(-0.2:(1/sample_rate):0.8)[1:n_samples]

            data = DataFrame(time = time, C3 = sin.(2π .* 10 .* time), C4 = cos.(2π .* 10 .* time))

            erp = eegfun.ErpData(
                "test_data",
                1,
                "condition_1",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                sample_rate,
                eegfun.AnalysisInfo(),
                50,  # 50 epochs averaged
            )

            # Resample by factor of 2
            resampled = eegfun.resample(erp, 2)

            # Check sample rate
            @test resampled.sample_rate == 256
            @test erp.sample_rate == 512  # Original unchanged

            # Check number of samples
            @test nrow(resampled.data) == 256

            # Check n_epochs preserved
            @test resampled.n_epochs == 50
            @test erp.n_epochs == 50

            # Check time vector
            @test resampled.data.time[1] == erp.data.time[1]
            @test resampled.data.time[2] == erp.data.time[3]
        end

        @testset "ERP with metadata columns" begin
            # Create ERP with condition info
            data = DataFrame(time = collect(0:511) ./ 512, C3 = randn(512))

            erp = eegfun.ErpData(
                "test_data",
                1,
                "Target",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                512,
                eegfun.AnalysisInfo(),
                30,
            )

            resampled = eegfun.resample(erp, 4)

            @test resampled.sample_rate == 128
            @test nrow(resampled.data) == 128
            @test resampled.condition_name == "Target"  # condition_name is in struct, not DataFrame
        end
    end


    #=============================================================================
        ERROR HANDLING
    =============================================================================#

    @testset "Error handling" begin

        @testset "Invalid factors" begin
            continuous = eegfun.ContinuousData(
                "test_data",
                DataFrame(time = [0.0, 0.001], C3 = [1.0, 2.0]),
                eegfun.Layout(DataFrame(), nothing, nothing),
                1000,
                eegfun.AnalysisInfo(),
            )

            # Zero factor
            @test_throws Exception eegfun.resample(continuous, 0)

            # Negative factor
            @test_throws Exception eegfun.resample(continuous, -1)
        end

        @testset "Non-divisible sample rate" begin
            # 500 Hz cannot be evenly divided by 3
            continuous = eegfun.ContinuousData(
                "test_data",
                DataFrame(time = [0.0, 0.002], C3 = [1.0, 2.0]),
                eegfun.Layout(DataFrame(), nothing, nothing),
                500,
                eegfun.AnalysisInfo(),
            )

            @test_throws Exception eegfun.resample(continuous, 3)
        end
    end


    #=============================================================================
        BATCH PROCESSING
    =============================================================================#

    @testset "Batch processing" begin

        @testset "Batch resample continuous data" begin
            mktempdir() do tmpdir
                # Create test files
                for i = 1:3
                    data = DataFrame(time = collect(0:511) ./ 512, C3 = randn(512), C4 = randn(512))

                    continuous = eegfun.ContinuousData(
                        "test_data",
                        data,
                        eegfun.Layout(DataFrame(), nothing, nothing),
                        512,
                        eegfun.AnalysisInfo(),
                    )

                    jldsave(joinpath(tmpdir, "$(i)_continuous.jld2"); data = continuous)
                end

                # Create output directory
                output_dir = joinpath(tmpdir, "resampled")

                # Batch resample
                eegfun.resample("continuous", 2, input_dir = tmpdir, output_dir = output_dir)

                # Check output files exist
                @test isfile(joinpath(output_dir, "1_continuous.jld2"))
                @test isfile(joinpath(output_dir, "2_continuous.jld2"))
                @test isfile(joinpath(output_dir, "3_continuous.jld2"))

                # Load and check one file
                resampled = load(joinpath(output_dir, "1_continuous.jld2"), "data")
                @test resampled isa eegfun.ContinuousData
                @test resampled.sample_rate == 256
                @test nrow(resampled.data) == 256
            end
        end

        @testset "Batch resample epochs" begin
            mktempdir() do tmpdir
                # Create test files
                for i = 1:2
                    epochs = []
                    for j = 1:3
                        epoch = DataFrame(time = collect(0:255) ./ 256, C3 = randn(256))
                        push!(epochs, epoch)
                    end

                    epoch_data = eegfun.EpochData(
                        "test_data",
                        1,
                        "condition_1",
                        epochs,
                        eegfun.Layout(DataFrame(), nothing, nothing),
                        256,
                        eegfun.AnalysisInfo(),
                    )

                    jldsave(joinpath(tmpdir, "$(i)_epochs.jld2"); data = epoch_data)
                end

                output_dir = joinpath(tmpdir, "resampled")

                # Batch resample
                eegfun.resample("epochs", 2, input_dir = tmpdir, output_dir = output_dir)

                # Check output
                @test isfile(joinpath(output_dir, "1_epochs.jld2"))

                # Load and verify
                resampled = load(joinpath(output_dir, "1_epochs.jld2"), "data")
                @test resampled isa eegfun.EpochData
                @test resampled.sample_rate == 128
                @test nrow(resampled.data[1]) == 128
            end
        end

        @testset "Participant filtering" begin
            mktempdir() do tmpdir
                # Create files for participants 1-5
                for i = 1:5
                    data = DataFrame(time = collect(0:255) ./ 256, C3 = randn(256))
                    continuous = eegfun.ContinuousData(
                        "test_data",
                        data,
                        eegfun.Layout(DataFrame(), nothing, nothing),
                        256,
                        eegfun.AnalysisInfo(),
                    )
                    jldsave(joinpath(tmpdir, "$(i)_continuous.jld2"); data = continuous)
                end

                output_dir = joinpath(tmpdir, "resampled")

                # Process only participants 2 and 4
                eegfun.resample(
                    "continuous",
                    2,
                    input_dir = tmpdir,
                    participant_selection = eegfun.participants([2, 4]),
                    output_dir = output_dir,
                )

                # Check only selected participants were processed
                @test !isfile(joinpath(output_dir, "1_continuous.jld2"))
                @test isfile(joinpath(output_dir, "2_continuous.jld2"))
                @test !isfile(joinpath(output_dir, "3_continuous.jld2"))
                @test isfile(joinpath(output_dir, "4_continuous.jld2"))
                @test !isfile(joinpath(output_dir, "5_continuous.jld2"))
            end
        end
    end


    #=============================================================================
        ANALYSIS INFO PRESERVATION
    =============================================================================#

    @testset "Analysis info preservation" begin

        @testset "Continuous data analysis info" begin
            data = DataFrame(time = collect(0:511) ./ 512, C3 = randn(512))

            analysis_info = eegfun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 40.0)

            continuous = eegfun.ContinuousData(
                "test_data",
                data,
                eegfun.Layout(DataFrame(), nothing, nothing),
                512,
                analysis_info,
            )

            resampled = eegfun.resample(continuous, 2)

            # Analysis info should be preserved
            @test resampled.analysis_info.reference == :avg
            @test resampled.analysis_info.hp_filter == 0.1
            @test resampled.analysis_info.lp_filter == 40.0
        end
    end

end # @testset "Resample tests"
