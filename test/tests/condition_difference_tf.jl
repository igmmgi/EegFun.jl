using Test
using DataFrames

@testset "TF Condition Difference" begin

    @testset "Basic difference" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        @test length(diff) == 1
        @test diff[1].condition == 1  # pair index
        @test diff[1].condition_name == "difference_1_2"
        @test diff[1] isa EegFun.TimeFreqData
    end

    @testset "Multiple difference pairs" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1),
            EegFun.create_test_tf_data(condition = 2),
            EegFun.create_test_tf_data(condition = 3),
            EegFun.create_test_tf_data(condition = 4),
        ]

        diff = EegFun.condition_difference(tfs, [(1, 2), (3, 4)])

        @test length(diff) == 2
        @test diff[1].condition_name == "difference_1_2"
        @test diff[2].condition_name == "difference_3_4"
    end

    @testset "Vector syntax" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [[1, 2]])

        @test length(diff) == 1
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "Power data integrity" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1, power_offset = 6.0, noise = 0.0),
            EegFun.create_test_tf_data(condition = 2, power_offset = 4.0, noise = 0.0),
        ]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        # Difference of 6.0 - 4.0 should be 2.0
        for ch in [:Ch1, :Ch2, :Ch3]
            @test all(diff[1].data_power[!, ch] .≈ 2.0)
        end
    end

    @testset "Phase data integrity" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        # Phase should be subtracted too
        for ch in [:Ch1, :Ch2, :Ch3]
            expected_phase = tfs[1].data_phase[!, ch] .- tfs[2].data_phase[!, ch]
            @test all(abs.(diff[1].data_phase[!, ch] .- expected_phase) .< 1e-10)
        end
    end

    @testset "Metadata preservation" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        diff = EegFun.condition_difference(tfs, [(1, 2)])

        @test diff[1].sample_rate == tfs[1].sample_rate
        @test diff[1].method == tfs[1].method
        @test nrow(diff[1].data_power) == nrow(tfs[1].data_power)
        @test nrow(diff[1].data_phase) == nrow(tfs[1].data_phase)

        # time and freq columns should be preserved
        @test diff[1].data_power.time == tfs[1].data_power.time
        @test diff[1].data_power.freq == tfs[1].data_power.freq
    end

    @testset "Missing condition skipping" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # Pair (3, 4) should be skipped, but (1, 2) should succeed
        diff = EegFun.condition_difference(tfs, [(1, 2), (3, 4)])

        @test length(diff) == 1
        @test diff[1].condition_name == "difference_1_2"
    end

    @testset "All conditions missing" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # No valid pairs → should throw
        @test_throws Exception EegFun.condition_difference(tfs, [(5, 6)])
    end
end
