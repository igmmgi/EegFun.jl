using Test
using Random

@testset "TF Decoding" begin

    # ────────────────────────────────────────────────────────────
    # 1. _prepare_decoding_data (TimeFreqEpochData)
    # ────────────────────────────────────────────────────────────
    @testset "_prepare_decoding_data (TF)" begin
        Random.seed!(200)
        n_channels = 3
        n_epochs = 8
        frequencies = [4.0, 8.0, 12.0]
        time_points = [0.0, 0.1, 0.2, 0.3]

        cond1 = EegFun.create_test_tf_epoch_data(
            condition = 1,
            condition_name = "cond_1",
            n_channels = n_channels,
            n_epochs = n_epochs,
            frequencies = frequencies,
            time_points = time_points,
            power_offset = 5.0,
            noise = 0.1,
        )
        cond2 = EegFun.create_test_tf_epoch_data(
            condition = 2,
            condition_name = "cond_2",
            n_channels = n_channels,
            n_epochs = n_epochs,
            frequencies = frequencies,
            time_points = time_points,
            power_offset = 3.0,
            noise = 0.1,
        )

        data_arrays, n_trials = EegFun._prepare_decoding_data([cond1, cond2])

        @test length(data_arrays) == 2
        @test length(n_trials) == 2
        @test n_trials[1] == n_epochs
        @test n_trials[2] == n_epochs

        # Shape: [(channels × frequencies) × time × trials]
        n_features = n_channels * length(frequencies)
        n_times = length(time_points)
        @test size(data_arrays[1]) == (n_features, n_times, n_epochs)
        @test size(data_arrays[2]) == (n_features, n_times, n_epochs)

        # Data should be finite
        @test all(isfinite, data_arrays[1])
        @test all(isfinite, data_arrays[2])
    end

    @testset "_prepare_decoding_data (TF) - unequal trials" begin
        Random.seed!(201)
        cond1 = EegFun.create_test_tf_epoch_data(condition = 1, n_epochs = 10, power_offset = 5.0)
        cond2 = EegFun.create_test_tf_epoch_data(condition = 2, n_epochs = 6, power_offset = 3.0)

        data_arrays, n_trials = EegFun._prepare_decoding_data([cond1, cond2])

        @test n_trials[1] == 10
        @test n_trials[2] == 6
        @test size(data_arrays[1], 3) == 10
        @test size(data_arrays[2], 3) == 6
        # Feature dimension should be same for both
        @test size(data_arrays[1], 1) == size(data_arrays[2], 1)
    end


    # ────────────────────────────────────────────────────────────
    # 2. decode_libsvm (TimeFreqEpochData)
    # ────────────────────────────────────────────────────────────
    @testset "decode_libsvm (TF) - basic" begin
        Random.seed!(210)
        # Create TF epoch data with clear condition difference
        cond1 = EegFun.create_test_tf_epoch_data(
            condition = 1,
            condition_name = "cond_1",
            n_channels = 3,
            n_epochs = 20,
            frequencies = [4.0, 8.0, 12.0],
            time_points = [0.0, 0.1, 0.2, 0.3],
            power_offset = 10.0,
            noise = 0.5,
        )
        cond2 = EegFun.create_test_tf_epoch_data(
            condition = 2,
            condition_name = "cond_2",
            n_channels = 3,
            n_epochs = 20,
            frequencies = [4.0, 8.0, 12.0],
            time_points = [0.0, 0.1, 0.2, 0.3],
            power_offset = 0.0,
            noise = 0.5,
        )

        result = EegFun.decode_libsvm([cond1, cond2]; n_iterations = 2, n_folds = 2)

        @test result isa EegFun.DecodedData
        @test length(result.times) == 4  # time_points
        @test length(result.average_score) == 4
        @test result.parameters.n_classes == 2
        @test result.parameters.n_iterations == 2
        @test result.parameters.n_folds == 2

        # With clear separation, accuracy should be above chance (0.5) at most time points
        @test mean(result.average_score) > 0.4
    end

    @testset "decode_libsvm (TF) - vector of participants" begin
        Random.seed!(211)

        # 2 participants, each with 2 conditions
        participants = Vector{Vector{EegFun.TimeFreqEpochData}}()
        for p = 1:2
            cond1 = EegFun.create_test_tf_epoch_data(participant = p, condition = 1, n_epochs = 12, power_offset = 8.0, noise = 1.0)
            cond2 = EegFun.create_test_tf_epoch_data(participant = p, condition = 2, n_epochs = 12, power_offset = 0.0, noise = 1.0)
            push!(participants, [cond1, cond2])
        end

        results = EegFun.decode_libsvm(participants; n_iterations = 2, n_folds = 2)

        @test length(results) == 2
        @test all(r -> r isa EegFun.DecodedData, results)
        @test results[1].times == results[2].times
    end


end # @testset "TF Decoding"
