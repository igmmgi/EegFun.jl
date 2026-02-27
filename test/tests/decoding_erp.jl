using Test
using Random
using Statistics

@testset "ERP Decoding" begin

    # ────────────────────────────────────────────────────────────
    # 1. _prepare_decoding_data (EpochData)
    # ────────────────────────────────────────────────────────────
    @testset "_prepare_decoding_data (ERP)" begin
        Random.seed!(300)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 8, n_channels = 3, n = 100, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 8, n_channels = 3, n = 100, fs = 1000)

        data_arrays, n_trials = EegFun._prepare_decoding_data([cond1, cond2])

        @test length(data_arrays) == 2
        @test length(n_trials) == 2
        @test n_trials[1] == 8
        @test n_trials[2] == 8

        # Shape: [channels × time × trials]
        @test size(data_arrays[1]) == (3, 100, 8)
        @test size(data_arrays[2]) == (3, 100, 8)

        # Data should be finite
        @test all(isfinite, data_arrays[1])
        @test all(isfinite, data_arrays[2])
    end

    @testset "_prepare_decoding_data (ERP) - unequal trials" begin
        Random.seed!(301)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 50, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 6, n_channels = 3, n = 50, fs = 1000)

        data_arrays, n_trials = EegFun._prepare_decoding_data([cond1, cond2])

        @test n_trials[1] == 10
        @test n_trials[2] == 6
        @test size(data_arrays[1], 3) == 10
        @test size(data_arrays[2], 3) == 6
        # Channel dimension should be same for both
        @test size(data_arrays[1], 1) == size(data_arrays[2], 1)
    end

    # ────────────────────────────────────────────────────────────
    # 2. _equalize_trials
    # ────────────────────────────────────────────────────────────
    @testset "_equalize_trials" begin
        Random.seed!(302)
        # Create arrays with different trial counts
        data1 = randn(3, 10, 12)  # 12 trials
        data2 = randn(3, 10, 8)   # 8 trials
        data_arrays = [data1, data2]
        n_trials = [12, 8]

        eq_data, eq_trials = EegFun._equalize_trials(data_arrays, n_trials)

        @test eq_trials == [8, 8]
        @test size(eq_data[1], 3) == 8
        @test size(eq_data[2], 3) == 8
    end

    # ────────────────────────────────────────────────────────────
    # 3. decode_libsvm (EpochData)
    # ────────────────────────────────────────────────────────────
    @testset "decode_libsvm (ERP) - basic" begin
        Random.seed!(310)
        # Create epoch data with distinguishable conditions (amplitude differs by condition)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 20, n_channels = 3, n = 50, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 20, n_channels = 3, n = 50, fs = 1000)

        result = EegFun.decode_libsvm([cond1, cond2]; n_iterations = 2, n_folds = 2)

        @test result isa EegFun.DecodedData
        @test length(result.times) == 50
        @test length(result.average_score) == 50
        @test result.parameters.n_classes == 2
        @test result.parameters.n_iterations == 2
        @test result.parameters.n_folds == 2
        @test result.parameters.chance_level ≈ 0.5
        @test result.parameters.class_coding == :binary

        # All scores should be valid probabilities
        @test all(0 .<= result.average_score .<= 1)

        # Should have standard error
        @test !isnothing(result.stderror)
        @test length(result.stderror) == 50

        # Should have confusion matrix
        @test !isnothing(result.confusion_matrix)
        @test size(result.confusion_matrix) == (50, 2, 2)

        # Condition names should be set
        @test result.condition_names == ["condition_1", "condition_2"]

        # Channels should be set
        @test length(result.channels) == 3
    end

    @testset "decode_libsvm (ERP) - channel selection" begin
        Random.seed!(311)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)

        result = EegFun.decode_libsvm([cond1, cond2]; channel_selection = EegFun.channels([:Ch1, :Ch2]), n_iterations = 2, n_folds = 2)

        @test length(result.channels) == 2
        @test result.channels == [:Ch1, :Ch2]
    end

    @testset "decode_libsvm (ERP) - interval selection" begin
        Random.seed!(312)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 12, n_channels = 3, n = 100, fs = 1000)
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 12, n_channels = 3, n = 100, fs = 1000)

        # Select a subset of time points  
        result = EegFun.decode_libsvm([cond1, cond2]; interval_selection = EegFun.times(0.01, 0.05), n_iterations = 2, n_folds = 2)

        @test length(result.times) < 100  # Should be smaller than full time range
        @test all(0.01 .<= result.times .<= 0.05)
    end

    @testset "decode_libsvm (ERP) - error cases" begin
        Random.seed!(313)
        cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 10, n_channels = 3, n = 50, fs = 1000)

        # Need at least 2 conditions
        @test_throws Exception EegFun.decode_libsvm([cond1]; n_iterations = 2, n_folds = 2)

        # Need at least 2 folds
        cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 10, n_channels = 3, n = 50, fs = 1000)
        @test_throws Exception EegFun.decode_libsvm([cond1, cond2]; n_iterations = 2, n_folds = 1)

        # Empty epochs
        @test_throws Exception EegFun.decode_libsvm(EegFun.EpochData[]; n_iterations = 2, n_folds = 2)
    end

    # ────────────────────────────────────────────────────────────
    # 4. Batch decode_libsvm (Vector of participants)
    # ────────────────────────────────────────────────────────────
    @testset "decode_libsvm (ERP) - vector of participants" begin
        Random.seed!(320)

        # 2 participants, each with 2 conditions
        participants = Vector{Vector{EegFun.EpochData}}()
        for p = 1:2
            cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
            cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
            push!(participants, [cond1, cond2])
        end

        results = EegFun.decode_libsvm(participants; n_iterations = 2, n_folds = 2)

        @test length(results) == 2
        @test all(r -> r isa EegFun.DecodedData, results)
        @test results[1].times == results[2].times
    end

    # ────────────────────────────────────────────────────────────
    # 5. grand_average (DecodedData)
    # ────────────────────────────────────────────────────────────
    @testset "grand_average (DecodedData)" begin
        Random.seed!(330)

        # Decode 3 participants
        all_decoded = EegFun.DecodedData[]
        for p = 1:3
            cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
            cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 12, n_channels = 3, n = 50, fs = 1000)
            decoded = EegFun.decode_libsvm([cond1, cond2]; n_iterations = 2, n_folds = 2, show_progress = false)
            push!(all_decoded, decoded)
        end

        ga = EegFun.grand_average(all_decoded)

        @test ga isa EegFun.DecodedData
        @test ga.file == "grand_average"
        @test length(ga.times) == length(all_decoded[1].times)
        @test length(ga.average_score) == length(ga.times)
        @test !isnothing(ga.stderror)
        @test all(0 .<= ga.average_score .<= 1)

        # Grand average should be mean of individual scores
        expected_avg = mean([d.average_score for d in all_decoded])
        @test ga.average_score ≈ expected_avg

        # Single participant should return itself
        single_ga = EegFun.grand_average([all_decoded[1]])
        @test single_ga === all_decoded[1]
    end

end # @testset "ERP Decoding"
