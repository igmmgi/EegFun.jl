"""
Test suite for EEGLAB .set file import functionality
"""

using Test
using EegFun
using DataFrames

@testset "EEGLAB Import" begin
    # Test file path
    test_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "eeglab", "epochs.set")

    if !isfile(test_file)
        EegFun.@minimal_warning "Test file not found: $test_file"
        EegFun.@minimal_warning "Skipping EEGLAB import tests"
        return
    end

    # read_eeglab returns (EpochData, InfoIca) when ICA is present
    result = EegFun.read_eeglab(test_file)

    @testset "Basic Loading" begin
        # Result is a tuple when ICA data is present
        @test result isa Tuple
        epochs, ica_info = result

        @test epochs isa EegFun.EpochData
        @test ica_info isa EegFun.InfoIca

        # Test dimensions
        @test length(epochs.data) == 267

        # Each epoch should be a DataFrame
        @test epochs.data[1] isa DataFrame

        # Should have time column
        @test :time in propertynames(epochs.data[1])

        # Test sample rate
        @test epochs.sample_rate == 125

        # Test file path stored
        @test occursin("epochs.set", epochs.file)
    end

    @testset "Channel Information" begin
        epochs = result[1]

        # Get channel names (excluding metadata columns)
        ch_labels = EegFun.channel_labels(epochs)

        # Test number of channels
        @test length(ch_labels) == 71

        # Test some expected channels
        @test :OZ in ch_labels || :Oz in ch_labels
        @test :FP1 in ch_labels || :Fp1 in ch_labels

        # Test layout
        @test epochs.layout isa EegFun.Layout
        @test size(epochs.layout.data, 1) == 71
    end

    @testset "Time Vector" begin
        epochs = result[1]

        # Get time vector from first epoch
        times = epochs.data[1].time

        # Test time range (approximately -1 to 2 seconds)
        @test minimum(times) ≈ -1.0 atol = 0.01
        @test maximum(times) ≈ 1.992 atol = 0.01

        # Test number of timepoints per epoch
        @test length(times) == 375
    end

    @testset "Data Integrity" begin
        epochs = result[1]

        # All epochs should have same dimensions
        n_rows = nrow(epochs.data[1])
        n_cols = ncol(epochs.data[1])

        for epoch in epochs.data
            @test nrow(epoch) == n_rows
            @test ncol(epoch) == n_cols
        end

        # Data should be numeric
        first_epoch = epochs.data[1]
        ch_labels = EegFun.channel_labels(epochs)
        ch_name = string(ch_labels[1])
        @test eltype(first_epoch[!, ch_name]) <: Real

        # Should not have NaN or Inf
        @test !any(isnan.(first_epoch[!, ch_name]))
        @test !any(isinf.(first_epoch[!, ch_name]))
    end

    @testset "ICA Data" begin
        ica_info = result[2]

        @test ica_info isa EegFun.InfoIca
        @test size(ica_info.unmixing) == (71, 71)
        @test size(ica_info.mixing) == (71, 71)
        @test size(ica_info.sphere) == (71, 71)
        @test length(ica_info.ica_label) == 71
    end

    @testset "Integration with EegFun Functions" begin
        epochs = result[1]

        # Test averaging works
        erp = EegFun.average_epochs(epochs)
        @test erp isa EegFun.ErpData

        # Test filtering works
        epochs_copy = copy(epochs)
        EegFun.lowpass_filter!(epochs_copy, 30.0)
        @test epochs_copy isa EegFun.EpochData

        # Test baseline works
        epochs_copy2 = copy(epochs)
        EegFun.baseline!(epochs_copy2, (-1.0, 0.0))
        @test epochs_copy2 isa EegFun.EpochData
    end
end

@testset "EEGLAB ERP Import" begin

    n_channels = 3
    n_timepoints = 100
    sample_rate = 250

    # Build a minimal EEGLAB-like dict
    mock_eeg = Dict{String,Any}(
        "nbchan" => n_channels,
        "pnts" => n_timepoints,
        "srate" => sample_rate,
        "trials" => 1,
        "xmin" => -0.2,
        "setname" => "test_erp",
        "data" => randn(Float32, n_channels, n_timepoints),
        "chanlocs" => Dict{String,Any}("labels" => ["Fz", "Cz", "Pz"], "theta" => [0.0, 0.0, 0.0], "radius" => [0.5, 0.3, 0.1]),
        "times" => collect(range(-200.0, step = 1000.0 / sample_rate, length = n_timepoints)),  # in ms
    )

    # This should NOT throw MethodError (the bug was missing n_epochs)
    erp = EegFun._eeglab_to_erpdata(mock_eeg, "test.set", true)

    @test erp isa EegFun.ErpData
    @test erp.n_epochs == 1
    @test erp.sample_rate == sample_rate
    @test nrow(erp.data) == n_timepoints
    @test length(EegFun.channel_labels(erp)) == n_channels
    @test erp.condition_name == "test_erp"
end
