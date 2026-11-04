"""
Test suite for FieldTrip .mat file import functionality
"""

using Test
using EegFun

@testset "FieldTrip Import" begin
    # Test file path
    test_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "fieldtrip", "continuous.mat")
    @test isfile(test_file)

    layout_file = joinpath(dirname(dirname(@__DIR__)), "resources", "layouts", "biosemi", "biosemi64.csv")
    @test isfile(layout_file)
    layout = EegFun.read_layout(layout_file)

    @testset "Basic Loading Continuous" begin
        dat = EegFun.read_fieldtrip(test_file, layout)

        @test dat isa EegFun.ContinuousData
        @test size(dat.data, 1) > 0
        @test :time in propertynames(dat.data)
        @test :sample in propertynames(dat.data)

        ch_labels = EegFun.channel_labels(dat)
        @test length(ch_labels) > 0
        @test :Fp1 in ch_labels || :FP1 in ch_labels

        # Verify layout copy isolation
        dat.layout.data[1, :label] = :Mutated
        @test layout.data[1, :label] != :Mutated
    end

    @testset "Basic Loading Epochs" begin
        epoch_file = joinpath(dirname(dirname(@__DIR__)), "resources", "data", "fieldtrip", "epochs.mat")
        @test isfile(epoch_file)

        dat = EegFun.read_fieldtrip(epoch_file, layout)
        @test dat isa EegFun.EpochData
        @test length(dat.data) > 0
        @test :time in propertynames(dat.data[1])
        @test :sample in propertynames(dat.data[1])
        @test :epoch in propertynames(dat.data[1])
        @test dat.data[1].epoch[1] == 1
    end
end
