using Test
using EegFun
using DataFrames

@testset "Read Raw Data" begin
    # Base path for resources
    data_dir = joinpath(dirname(dirname(@__DIR__)), "resources", "data")

    @testset "BDF Import" begin
        test_file = joinpath(data_dir, "bdf", "example1.bdf")
        @test isfile(test_file)

        result = EegFun.read_raw_data(test_file)
        @test result isa EegFun.BiosemiDataFormat.BiosemiData
        @test size(result.data, 1) > 0
        @test length(result.time) > 0
        @test result.header.sample_rate[1] > 0

        # Test ContinuousData construction and layout copy isolation
        lay = EegFun.read_layout(joinpath(dirname(dirname(@__DIR__)), "resources", "layouts", "biosemi", "biosemi72.csv"))
        dat = EegFun.create_eegfun_data(result, lay)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) > 0
        @test size(EegFun.all_data(dat), 1) > 0
        dat.layout.data[1, :label] = :Mutated
        @test lay.data[1, :label] != :Mutated
    end

    @testset "BrainVision Import" begin
        test_file = joinpath(data_dir, "brainvision", "example1.vhdr")
        @test isfile(test_file)

        result = EegFun.read_raw_data(test_file)
        @test result isa EegFun.BrainVisionDataFormat.BrainVisionData
        @test size(result.data, 1) > 0
        @test result.header.Fs > 0

        # Test ContinuousData construction & numeric trigger code extraction ("S253" -> 253)
        dat = EegFun.create_eegfun_data(result)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) > 0
        @test size(EegFun.all_data(dat), 1) > 0
        @test 253 in dat.data.trigger
        @test any(s -> occursin("253", s), dat.data.trigger_info)

        # Verify 0-based marker sample index is not shifted across all markers
        m_proto = result.markers[1]
        m0 = EegFun.BrainVisionDataFormat.BrainVisionMarker(
            [
                f == :sample ? 0 : (f == :value ? "New Segment" : (f == :type ? "New Segment" : getfield(m_proto, f))) for
                f in fieldnames(typeof(m_proto))
            ]...,
        )
        m5 = EegFun.BrainVisionDataFormat.BrainVisionMarker(
            [
                f == :sample ? 5 : (f == :value ? "S 10" : (f == :type ? "Stimulus" : getfield(m_proto, f))) for
                f in fieldnames(typeof(m_proto))
            ]...,
        )
        t_col, ti_col = EegFun._extract_triggers_from_markers([m0, m5], 10)
        @test t_col[5] == 10
        @test ti_col[5] == "S 10"
    end

    @testset "EDF Import" begin
        test_file = joinpath(data_dir, "edf", "test.edf")
        @test isfile(test_file)

        result = EegFun.read_raw_data(test_file)
        @test result isa EegFun.EuropeanDataFormat.EdfData
        @test size(result.data, 1) > 0
        @test length(result.time) > 0
        @test result.header.sample_rate[1] > 0

        # Test ContinuousData construction & unit scaling ("mV" -> μV * 1e3)
        dat = EegFun.create_eegfun_data(result)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) > 0
        @test size(EegFun.all_data(dat), 1) > 0
        @test :Status ∉ EegFun.channel_labels(dat)
        dc_idx = findfirst(==("DC00.1"), strip.(result.header.channel_labels))
        @test !isnothing(dc_idx)
        @test lowercase(strip(result.header.channel_unit[dc_idx])) == "mv"
        @test isapprox(dat.data[!, Symbol("DC00.1")], result.data[:, dc_idx] .* 1e3; rtol = 1e-6)
    end

    @testset "FIF Import" begin
        test_file = joinpath(data_dir, "fif", "test_raw.fif")
        @test isfile(test_file)

        result = EegFun.read_raw_data(test_file)
        @test result isa EegFun.FunctionalImageFormat.FifData
        @test size(result.data, 1) > 0
        @test result.header.sfreq > 0

        # Test ContinuousData construction & Volt -> μV scaling
        dat = EegFun.create_eegfun_data(result)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) > 0
        @test size(EegFun.all_data(dat), 1) > 0
        if Symbol("EEG 001") in EegFun.channel_labels(dat)
            @test maximum(abs, dat.data[!, Symbol("EEG 001")]) > 1.0  # scaled from ~1e-4 V to ~100 μV
        end

        test_epo_file = joinpath(data_dir, "fif", "test_epochs.fif")
        @test isfile(test_epo_file)

        result_epo = EegFun.read_raw_data(test_epo_file)
        @test result_epo isa EegFun.FunctionalImageFormat.FifEpochs
        @test size(result_epo.data, 3) == 4  # 4 evoked conditions

        # Multi-condition Evoked FIF returns Vector{ErpData}
        epo_dat = EegFun.create_eegfun_data(result_epo)
        @test epo_dat isa Vector{EegFun.ErpData}
        @test length(epo_dat) == 4
        @test all(e -> EegFun.sample_rate(e) > 0, epo_dat)
        @test epo_dat[1].condition_name == "Left Auditory"
        @test "trigger" in names(epo_dat[1].data)
    end

    @testset "XDF Import" begin
        test_file = joinpath(data_dir, "xdf", "test1.xdf")
        @test isfile(test_file)

        result = EegFun.read_raw_data(test_file)
        @test result isa EegFun.ExtensibleDataFormat.XdfData

        eeg_streams = [s for s in values(result.streams) if s.header.type == "EEG"]
        @test !isempty(eeg_streams)
        @test size(eeg_streams[1].time_series, 1) > 0

        # Test ContinuousData construction
        dat = EegFun.create_eegfun_data(result)
        @test dat isa EegFun.ContinuousData
        @test EegFun.sample_rate(dat) == 250
        @test size(EegFun.all_data(dat), 1) > 0
        @test "trigger" in names(EegFun.all_data(dat))

        # Test trigger_count
        info = EegFun.trigger_count(dat)
        @test info isa EegFun.TriggerInfo
    end
end
