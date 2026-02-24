using Test
using DataFrames

@testset "TF Channel Average" begin

    @testset "Mutating: basic average all channels" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 3.0, noise = 0.0)

        EegFun.channel_average!(tf)

        # Should have added :avg column to both power and phase
        @test hasproperty(tf.data_power, :avg)
        @test hasproperty(tf.data_phase, :avg)

        # :avg should be the mean of Ch1, Ch2, Ch3 (all 3.0 with noise=0)
        @test all(tf.data_power[!, :avg] .≈ 3.0)
    end

    @testset "Mutating: specific channel group" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 2.0, noise = 0.0)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:frontal])

        @test hasproperty(tf.data_power, :frontal)
        @test hasproperty(tf.data_phase, :frontal)
        @test all(tf.data_power[!, :frontal] .≈ 2.0)
    end

    @testset "Mutating: reduce mode" begin
        tf = EegFun.create_test_tf_data(n_channels = 3, power_offset = 1.0, noise = 0.0)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:avg_group], reduce = true)

        # Should only have metadata + averaged columns
        @test hasproperty(tf.data_power, :time)
        @test hasproperty(tf.data_power, :freq)
        @test hasproperty(tf.data_power, :avg_group)
        @test !hasproperty(tf.data_power, :Ch1)
        @test !hasproperty(tf.data_power, :Ch3)

        # Same for phase
        @test hasproperty(tf.data_phase, :avg_group)
        @test !hasproperty(tf.data_phase, :Ch1)
    end

    @testset "Non-mutating: returns copy" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        result = EegFun.channel_average(tf)

        # Original should not be modified
        @test !hasproperty(tf.data_power, :avg)

        # Result should have the averaged column
        @test hasproperty(result.data_power, :avg)
        @test hasproperty(result.data_phase, :avg)
    end

    @testset "Non-mutating: Vector{TimeFreqData}" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        results = EegFun.channel_average(tfs)

        @test length(results) == 2
        @test all(hasproperty(r.data_power, :avg) for r in results)
        @test all(hasproperty(r.data_phase, :avg) for r in results)

        # Originals untouched
        @test !hasproperty(tfs[1].data_power, :avg)
    end

    @testset "Both DataFrames stay in sync" begin
        tf = EegFun.create_test_tf_data(n_channels = 3)

        EegFun.channel_average!(tf, channel_selections = [EegFun.channels([:Ch1, :Ch2])], output_labels = [:test_avg])

        # Same number of rows
        @test nrow(tf.data_power) == nrow(tf.data_phase)

        # Same columns added
        @test propertynames(tf.data_power) == propertynames(tf.data_phase)
    end
end
