using Test
using DataFrames
using EegFun

@testset "channel_delete" begin

    @testset "Mutating: ContinuousData" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        original_channels = EegFun.channel_labels(dat)
        @test :Ch1 ∈ original_channels

        EegFun.channel_delete!(dat, :Ch1)
        @test :Ch1 ∉ EegFun.channel_labels(dat)
        @test :Ch1 ∉ propertynames(dat.data)
        @test :Ch2 ∈ EegFun.channel_labels(dat)
        @test :Ch2 ∈ propertynames(dat.data)
        # Layout should be updated too
        @test :Ch1 ∉ dat.layout.data.label
    end

    @testset "Mutating: multiple channels" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        EegFun.channel_delete!(dat, [:Ch1, :Ch3])
        @test :Ch1 ∉ EegFun.channel_labels(dat)
        @test :Ch3 ∉ EegFun.channel_labels(dat)
        @test :Ch2 ∈ EegFun.channel_labels(dat)
        @test :Ch1 ∉ dat.layout.data.label
        @test :Ch3 ∉ dat.layout.data.label
    end

    @testset "Mutating: EpochData" begin
        dat = EegFun.create_test_epoch_data(n = 200, n_epochs = 3)
        EegFun.channel_delete!(dat, :Ch1)
        @test :Ch1 ∉ EegFun.channel_labels(dat)
        @test :Ch1 ∉ dat.layout.data.label
        # All epoch dataframes should have the channel removed
        for df in dat.data
            @test :Ch1 ∉ propertynames(df)
            @test :Ch2 ∈ propertynames(df)
        end
    end

    @testset "Mutating: ErpData" begin
        dat = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 3)
        EegFun.channel_delete!(dat, :Ch2)
        @test :Ch2 ∉ EegFun.channel_labels(dat)
        @test :Ch2 ∉ propertynames(dat.data)
        @test :Ch1 ∈ propertynames(dat.data)
        @test :Ch2 ∉ dat.layout.data.label
    end

    @testset "Non-mutating: returns copy" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        result = EegFun.channel_delete(dat, :Ch1)
        # Original should be unchanged
        @test :Ch1 ∈ EegFun.channel_labels(dat)
        @test :Ch1 ∈ propertynames(dat.data)
        # Result should have the channel removed
        @test :Ch1 ∉ EegFun.channel_labels(result)
        @test :Ch1 ∉ propertynames(result.data)
    end

    @testset "Non-mutating: EpochData" begin
        dat = EegFun.create_test_epoch_data(n = 200, n_epochs = 2)
        result = EegFun.channel_delete(dat, [:Ch1, :Ch3])
        # Original untouched
        @test :Ch1 ∈ propertynames(dat.data[1])
        @test :Ch3 ∈ propertynames(dat.data[1])
        # Result has channels removed
        for df in result.data
            @test :Ch1 ∉ propertynames(df)
            @test :Ch3 ∉ propertynames(df)
            @test :Ch2 ∈ propertynames(df)
        end
    end

    @testset "Nonexistent channel: no-op" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        n_cols_before = ncol(dat.data)
        n_layout_before = nrow(dat.layout.data)
        EegFun.channel_delete!(dat, :DoesNotExist)
        # Nothing should have changed
        @test ncol(dat.data) == n_cols_before
        @test nrow(dat.layout.data) == n_layout_before
    end

    @testset "Mix of existing and nonexistent channels" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        EegFun.channel_delete!(dat, [:Ch1, :DoesNotExist])
        @test :Ch1 ∉ EegFun.channel_labels(dat)
        @test :Ch2 ∈ EegFun.channel_labels(dat)
    end

    @testset "Metadata preserved after deletion" begin
        dat = EegFun.create_test_continuous_data(n = 200)
        original_time = copy(dat.data.time)
        original_sr = dat.sample_rate
        EegFun.channel_delete!(dat, :Ch1)
        @test dat.data.time == original_time
        @test dat.sample_rate == original_sr
    end

    @testset "Delete all EEG channels" begin
        dat = EegFun.create_test_continuous_data(n = 200, n_channels = 2)
        EegFun.channel_delete!(dat, [:Ch1, :Ch2])
        @test isempty(EegFun.channel_labels(dat))
        # Metadata columns should still exist
        @test :time ∈ propertynames(dat.data)
    end

end
