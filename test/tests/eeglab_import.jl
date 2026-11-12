"""
Test suite for EEGLAB .set file import functionality
"""

using Test
using EegFun
using DataFrames

@testset "EEGLAB Import" begin
    # Test file path
    test_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "eeglab", "epochs.set")
    @test isfile(test_file)

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

        # Should have time, sample, and epoch columns
        @test :time in propertynames(epochs.data[1])
        @test :sample in propertynames(epochs.data[1])
        @test :epoch in propertynames(epochs.data[1])
        @test epochs.data[1].epoch[1] == 1
        @test epochs.data[2].epoch[1] == 2

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
    @test :sample in propertynames(erp.data)
    @test length(EegFun.channel_labels(erp)) == n_channels
    @test erp.condition_name == "test_erp"
end

@testset "EEGLAB Continuous Triggers, Epoch Latency=0 Trigger, and ICA Sphering" begin
    # 1. Continuous triggers read event.type (with struct-of-arrays and array-of-dicts)
    mock_cont = Dict{String,Any}(
        "nbchan" => 3,
        "pnts" => 20,
        "srate" => 100,
        "trials" => 1,
        "xmin" => 0.0,
        "setname" => "cont_test",
        "data" => ones(Float32, 3, 20),
        "chanlocs" => Dict{String,Any}("labels" => ["Fz", "Cz", "Pz"], "theta" => [0.0, 0.0, 0.0], "radius" => [0.5, 0.3, 0.1]),
        "event" => [
            Dict{String,Any}("latency" => 3, "type" => "S 12", "position" => 99),
            Dict{String,Any}("latency" => 10, "type" => 25, "position" => 88),
        ],
    )
    cont = EegFun._eeglab_to_continuousdata(mock_cont, "cont.set", true)
    @test cont isa EegFun.ContinuousData
    @test :sample in propertynames(cont.data)
    @test cont.data.trigger[3] == 12
    @test cont.data.trigger[10] == 25

    # 2. Epoch _extract_triggers picks event closest to eventlatency == 0
    mock_ep_events = [
        Dict{String,Any}("eventlatency" => [-250.0, 0.0, 350.0], "eventtype" => ["cue", "S 42", "response"], "eventposition" => [1, 2, 3]),
        Dict{String,Any}("eventlatency" => [-100.0, 2.0], "eventtype" => [99, 77], "eventposition" => [1, 2]),
    ]
    trigs, trig_info = EegFun._extract_triggers(mock_ep_events, 2)
    @test trigs == [42, 77]
    @test trig_info == ["S 42", "77"]

    # 3. _extract_ica_info computes unmixing = W * sphere, mean length = size(unmixing, 2), and subsets layout to icachansind
    W = [1.0 2.0; 3.0 4.0]
    sph = [2.0 0.0; 0.0 0.5]
    mock_ica_eeg = Dict{String,Any}("icaweights" => W, "icasphere" => sph, "icawinv" => inv(W * sph), "icachansind" => [1, 3])
    lay3 = EegFun.Layout(DataFrame(label = [:Fz, :Cz, :Pz], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]), nothing, nothing, nothing)
    ica = EegFun._extract_ica_info(mock_ica_eeg, "ica_test.set", lay3)
    @test ica.unmixing ≈ W * sph
    @test length(ica.mean) == 2
    @test EegFun.channel_labels(ica.layout) == [:Fz, :Pz]
    ica.layout.data[1, :label] = :Mutated
    @test lay3.data[1, :label] == :Fz

    # 4. Single-channel chanlocs and single-event struct in _eeglab_extract_continuous_events!
    @test EegFun._extract_channel_names(Dict{String,Any}("labels" => "Cz")) == [:Cz]
    tv = zeros(Int, 10)
    tiv = fill("", 10)
    EegFun._eeglab_extract_continuous_events!(tv, tiv, Dict{String,Any}("latency" => 4.0, "type" => "S 9"), 10)
    @test tv[4] == 9
end
