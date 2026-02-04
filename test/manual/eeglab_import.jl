"""
Test suite for EEGLAB .set file import functionality
"""

using Test
using EegFun
using DataFrames

@testset "EEGLAB Import" begin
    # Test file path
    test_file = "./resources/data/eeglab/epochs.set"

    if !isfile(test_file)
        @warn "Test file not found: $test_file"
        @warn "Skipping EEGLAB import tests"
        return
    end

    @testset "Basic Loading" begin
        # Read data
        epochs = read_eeglab(test_file, verbose = false)

        # Test data type
        @test epochs isa EpochData

        # Test dimensions
        @test length(epochs.data) == 91  # 91 trials

        # Each epoch should be a DataFrame
        @test epochs.data[1] isa DataFrame

        # Should have time column
        @test :time in names(epochs.data[1])

        # Test sample rate
        @test epochs.sample_rate == 125

        # Test file path stored
        @test occursin("Probe.set", epochs.file)
    end

    @testset "Channel Information" begin
        epochs = read_eeglab(test_file, verbose = false)

        # Get channel names (excluding time)
        ch_names = filter(x -> x != :time, names(epochs.data[1]))

        # Test number of channels
        @test length(ch_names) == 71

        # Test some expected channels
        @test :OZ in ch_names || :Oz in ch_names
        @test :FP1 in ch_names || :Fp1 in ch_names

        # Test layout
        @test epochs.layout isa Layout
        @test size(epochs.layout.data, 1) == 71
    end

    @testset "Time Vector" begin
        epochs = read_eeglab(test_file, verbose = false)

        # Get time vector from first epoch
        times = epochs.data[1].time

        # Test time range (approximately -1 to 2 seconds)
        @test minimum(times) ≈ -1.0 atol = 0.01
        @test maximum(times) ≈ 1.992 atol = 0.01

        # Test number of timepoints
        @test length(times) == 375
    end

    @testset "Data Integrity" begin
        epochs = read_eeglab(test_file, verbose = false)

        # All epochs should have same dimensions
        n_rows = nrow(epochs.data[1])
        n_cols = ncol(epochs.data[1])

        for epoch in epochs.data
            @test nrow(epoch) == n_rows
            @test ncol(epoch) == n_cols
        end

        # Data should be numeric
        first_epoch = epochs.data[1]
        ch_name = names(first_epoch)[1]  # First channel
        @test eltype(first_epoch[!, ch_name]) <: Real

        # Should not have NaN or Inf
        @test !any(isnan.(first_epoch[!, ch_name]))
        @test !any(isinf.(first_epoch[!, ch_name]))
    end

    @testset "Integration with EegFun Functions" begin
        epochs = read_eeglab(test_file, verbose = false)

        # Test averaging works
        erp = average_epochs(epochs)
        @test erp isa ErpData

        # Test filtering works
        epochs_copy = copy(epochs)
        lowpass_filter!(epochs_copy, 30.0)
        @test epochs_copy isa EpochData

        # Test baseline works
        epochs_copy2 = copy(epochs)
        baseline!(epochs_copy2, (-1.0, 0.0))
        @test epochs_copy2 isa EpochData
    end
end
