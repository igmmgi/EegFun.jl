"""
Test suite for FieldTrip .mat file import functionality
"""

using Test
using EegFun
using DataFrames

@testset "FieldTrip Import" begin
    # Test file path
    test_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "fieldtrip", "continuous.mat")

    if !isfile(test_file)
        EegFun.@minimal_warning "Test file not found: $test_file"
        EegFun.@minimal_warning "Skipping FieldTrip import tests"
        return
    end

    layout_file = joinpath(dirname(dirname(@__DIR__)), "resources", "layouts", "biosemi", "biosemi64.csv")
    layout = EegFun.read_layout(layout_file)

    @testset "Basic Loading Continuous" begin
        dat = EegFun.read_fieldtrip(test_file, layout)

        @test dat isa EegFun.ContinuousData
        @test size(dat.data, 1) > 0
        @test :time in propertynames(dat.data)

        ch_labels = EegFun.channel_labels(dat)
        @test length(ch_labels) > 0
        @test :Fp1 in ch_labels || :FP1 in ch_labels
    end

    @testset "Basic Loading Epochs" begin
        epoch_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "fieldtrip", "epochs.mat")
        if isfile(epoch_file)
            dat = EegFun.read_fieldtrip(epoch_file, layout)
            @test dat isa EegFun.EpochData
            @test length(dat.data) > 0
            @test :time in propertynames(dat.data[1])
        end
    end

    @testset "read_fieldtrip_csv" begin
        using CSV
        mktempdir() do tmpdir
            # 1. Missing file error
            @test_throws ErrorException EegFun.read_fieldtrip_csv(tmpdir)

            # 2. Write valid epochs_data.csv
            csv_path = joinpath(tmpdir, "epochs_data.csv")
            df = DataFrame(
                trial = [1, 1, 1, 2, 2, 2],
                time = [0.0, 0.01, 0.02, 0.0, 0.01, 0.02],
                Cz = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
                Pz = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            )
            CSV.write(csv_path, df)

            ep = EegFun.read_fieldtrip_csv(tmpdir; fsample = 100)
            @test ep isa EegFun.EpochData
            @test length(ep.data) == 2
            @test :Cz in EegFun.channel_labels(ep)
            @test :Pz in EegFun.channel_labels(ep)
            @test ep.sample_rate == 100

            # Test automatic sample rate estimation
            ep_auto = EegFun.read_fieldtrip_csv(tmpdir)
            @test ep_auto.sample_rate == 100

            # 3. Invalid CSV header
            bad_csv = joinpath(tmpdir, "epochs_data.csv")
            CSV.write(bad_csv, DataFrame(a = [1], b = [2]))
            @test_throws ErrorException EegFun.read_fieldtrip_csv(tmpdir)
        end
    end
end
