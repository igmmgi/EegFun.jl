using Test
using DataFrames

@testset "TF Condition Average" begin

    @testset "Basic average" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1, condition_name = "cond_1"),
            EegFun.create_test_tf_data(condition = 2, condition_name = "cond_2"),
        ]

        avg = EegFun.condition_average(tfs, [[1, 2]])

        @test length(avg) == 1
        @test avg[1].condition == 1  # group index
        @test avg[1].condition_name == "average_1_2"
        @test avg[1] isa EegFun.TimeFreqData
    end

    @testset "Multiple average groups" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1),
            EegFun.create_test_tf_data(condition = 2),
            EegFun.create_test_tf_data(condition = 3),
            EegFun.create_test_tf_data(condition = 4),
        ]

        avg = EegFun.condition_average(tfs, [[1, 2], [3, 4]])

        @test length(avg) == 2
        @test avg[1].condition_name == "average_1_2"
        @test avg[2].condition_name == "average_3_4"
    end

    @testset "Average all conditions" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1),
            EegFun.create_test_tf_data(condition = 2),
            EegFun.create_test_tf_data(condition = 3),
        ]

        avg = EegFun.condition_average(tfs, [[1, 2, 3]])

        @test length(avg) == 1
        @test avg[1].condition_name == "average_1_2_3"
    end

    @testset "Power data integrity" begin
        tfs = [
            EegFun.create_test_tf_data(condition = 1, power_offset = 4.0, noise = 0.0),
            EegFun.create_test_tf_data(condition = 2, power_offset = 6.0, noise = 0.0),
        ]

        avg = EegFun.condition_average(tfs, [[1, 2]])

        # Average of 4.0 and 6.0 should be 5.0
        for ch in [:Ch1, :Ch2, :Ch3]
            @test all(avg[1].data_power[!, ch] .≈ 5.0)
        end
    end

    @testset "Phase data integrity" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        avg = EegFun.condition_average(tfs, [[1, 2]])

        # Phase should be averaged too
        for ch in [:Ch1, :Ch2, :Ch3]
            expected_phase = (tfs[1].data_phase[!, ch] .+ tfs[2].data_phase[!, ch]) ./ 2
            @test all(abs.(avg[1].data_phase[!, ch] .- expected_phase) .< 1e-10)
        end
    end

    @testset "Metadata preservation" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        avg = EegFun.condition_average(tfs, [[1, 2]])

        @test avg[1].sample_rate == tfs[1].sample_rate
        @test avg[1].method == tfs[1].method
        @test nrow(avg[1].data_power) == nrow(tfs[1].data_power)
        @test nrow(avg[1].data_phase) == nrow(tfs[1].data_phase)

        # time and freq columns should be preserved
        @test avg[1].data_power.time == tfs[1].data_power.time
        @test avg[1].data_power.freq == tfs[1].data_power.freq
    end

    @testset "Missing condition skipping" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # Group [3, 4] should be skipped, but [1, 2] should succeed
        avg = EegFun.condition_average(tfs, [[1, 2], [3, 4]])

        @test length(avg) == 1
        @test avg[1].condition_name == "average_1_2"
    end

    @testset "All conditions missing" begin
        tfs = [EegFun.create_test_tf_data(condition = 1), EegFun.create_test_tf_data(condition = 2)]

        # No valid groups → should throw
        @test_throws Exception EegFun.condition_average(tfs, [[5, 6]])
    end
end
