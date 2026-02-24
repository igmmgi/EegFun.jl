using Test
using DataFrames

@testset "TF Channel Difference" begin

    @testset "Mutating: basic difference" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 5.0, noise = 0.0)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :diff_12,
        )

        # Should have :diff_12 in both power and phase
        @test hasproperty(tf.data_power, :diff_12)
        @test hasproperty(tf.data_phase, :diff_12)

        # With same offset and no noise, difference should be ~0
        @test all(abs.(tf.data_power[!, :diff_12]) .< 1e-10)
    end

    @testset "Mutating: group difference" begin
        tf = EegFun.create_test_tf_data(n_channels = 4, power_offset = 3.0, noise = 0.0)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1, :Ch2]),
            channel_selection2 = EegFun.channels([:Ch3, :Ch4]),
            channel_out = :laterality,
        )

        @test hasproperty(tf.data_power, :laterality)
        @test hasproperty(tf.data_phase, :laterality)
    end

    @testset "Non-mutating: returns copy" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        result = EegFun.channel_difference(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :test_diff,
        )

        # Original should not be modified
        @test !hasproperty(tf.data_power, :test_diff)

        # Result should have the diff column
        @test hasproperty(result.data_power, :test_diff)
        @test hasproperty(result.data_phase, :test_diff)
    end

    @testset "Non-mutating: Vector{TimeFreqData}" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        results = EegFun.channel_difference(
            tfs,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :diff,
        )

        @test length(results) == 2
        @test all(hasproperty(r.data_power, :diff) for r in results)
        @test all(hasproperty(r.data_phase, :diff) for r in results)

        # Originals untouched
        @test !hasproperty(tfs[1].data_power, :diff)
    end

    @testset "Both DataFrames stay in sync" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        EegFun.channel_difference!(
            tf,
            channel_selection1 = EegFun.channels([:Ch1]),
            channel_selection2 = EegFun.channels([:Ch2]),
            channel_out = :sync_test,
        )

        # Same number of rows
        @test nrow(tf.data_power) == nrow(tf.data_phase)

        # Same columns
        @test propertynames(tf.data_power) == propertynames(tf.data_phase)
    end
end
