using Test
using Random
using Statistics

@testset "Decoding Statistics" begin

    # ────────────────────────────────────────────────────────────
    # Helper to create decoded data for statistics testing
    # ────────────────────────────────────────────────────────────
    function _create_test_decoded_list(; n_participants = 5, n_iterations = 5, n_folds = 2)
        Random.seed!(400)  # Seed once — each participant gets different random data
        decoded_list = EegFun.DecodedData[]
        for p = 1:n_participants
            cond1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 16, n_channels = 3, n = 20, fs = 1000)
            cond2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 16, n_channels = 3, n = 20, fs = 1000)
            decoded = EegFun.decode_libsvm([cond1, cond2]; n_iterations = n_iterations, n_folds = n_folds, show_progress = false)
            push!(decoded_list, decoded)
        end
        return decoded_list
    end

    # ────────────────────────────────────────────────────────────
    # 1. test_against_chance
    # ────────────────────────────────────────────────────────────
    @testset "test_against_chance - uncorrected" begin
        decoded_list = _create_test_decoded_list(n_participants = 5)

        stats = EegFun.test_against_chance(decoded_list; alpha = 0.05, correction_method = :none)

        @test stats isa EegFun.DecodingStatisticsResult
        @test length(stats.times) == length(decoded_list[1].times)
        @test length(stats.t_statistics) == length(stats.times)
        @test length(stats.p_values) == length(stats.times)
        @test length(stats.significant_mask) == length(stats.times)
        @test stats.df ≈ 4.0  # n_participants - 1
        @test stats.alpha == 0.05
        @test stats.correction_method == :none
        @test isnothing(stats.clusters)  # No clusters for simple t-test

        # P-values should be between 0 and 1 (NaN is valid for zero-variance time points)
        @test all(p -> isnan(p) || (0 <= p <= 1), stats.p_values)
    end

    @testset "test_against_chance - Bonferroni" begin
        decoded_list = _create_test_decoded_list(n_participants = 5)

        stats = EegFun.test_against_chance(decoded_list; alpha = 0.05, correction_method = :bonferroni)

        @test stats isa EegFun.DecodingStatisticsResult
        @test stats.correction_method == :bonferroni
        @test isnothing(stats.clusters)

        # Bonferroni should be more conservative (fewer or equal significant points)
        stats_uncorr = EegFun.test_against_chance(decoded_list; alpha = 0.05, correction_method = :none)
        @test sum(stats.significant_mask) <= sum(stats_uncorr.significant_mask)
    end

    @testset "test_against_chance - error cases" begin
        # Empty list
        @test_throws Exception EegFun.test_against_chance(EegFun.DecodedData[])

        # Invalid correction method
        decoded_list = _create_test_decoded_list(n_participants = 3)
        @test_throws Exception EegFun.test_against_chance(decoded_list; correction_method = :invalid)
    end

    # ────────────────────────────────────────────────────────────
    # 2. test_against_chance_cluster
    # ────────────────────────────────────────────────────────────
    @testset "test_against_chance_cluster - basic" begin
        decoded_list = _create_test_decoded_list(n_participants = 5)

        stats = EegFun.test_against_chance_cluster(decoded_list; alpha = 0.05, n_permutations = 20, show_progress = false)

        @test stats isa EegFun.DecodingStatisticsResult
        @test length(stats.times) == length(decoded_list[1].times)
        @test length(stats.t_statistics) == length(stats.times)
        @test length(stats.p_values) == length(stats.times)
        @test stats.correction_method == :cluster_permutation
        @test !isnothing(stats.clusters)

        # P-values should be between 0 and 1 (NaN is valid for zero-variance time points)
        @test all(p -> isnan(p) || (0 <= p <= 1), stats.p_values)

        # All clusters should have valid fields
        for cluster in stats.clusters
            @test cluster isa EegFun.TemporalCluster
            @test cluster.id > 0
            @test !isempty(cluster.time_indices)
            @test 0 <= cluster.p_value <= 1
            @test cluster.time_range[1] <= cluster.time_range[2]
        end
    end

    @testset "test_against_chance_cluster - max statistic" begin
        decoded_list = _create_test_decoded_list(n_participants = 5)

        stats = EegFun.test_against_chance_cluster(
            decoded_list;
            alpha = 0.05,
            n_permutations = 10,
            cluster_statistic = :max,
            show_progress = false,
        )

        @test stats isa EegFun.DecodingStatisticsResult
        @test stats.correction_method == :cluster_permutation
    end

    @testset "test_against_chance_cluster - error cases" begin
        # Empty list
        @test_throws Exception EegFun.test_against_chance_cluster(EegFun.DecodedData[])

        # Invalid cluster statistic
        decoded_list = _create_test_decoded_list(n_participants = 3)
        @test_throws Exception EegFun.test_against_chance_cluster(decoded_list; cluster_statistic = :invalid)
    end

    # ────────────────────────────────────────────────────────────
    # 3. _find_temporal_clusters
    # ────────────────────────────────────────────────────────────
    @testset "_find_temporal_clusters" begin
        times = collect(0.0:0.1:0.9)  # 10 time points

        # No significant points
        mask = falses(10)
        clusters = EegFun._find_temporal_clusters(mask, times)
        @test isempty(clusters)

        # Single cluster
        mask = falses(10)
        mask[3:5] .= true
        clusters = EegFun._find_temporal_clusters(mask, times)
        @test length(clusters) >= 1
        # Check that the cluster covers the correct indices
        all_indices = vcat([c.time_indices for c in clusters]...)
        @test 3 ∈ all_indices
        @test 4 ∈ all_indices
        @test 5 ∈ all_indices

        # Two separate clusters
        mask = falses(10)
        mask[2:3] .= true
        mask[7:9] .= true
        clusters = EegFun._find_temporal_clusters(mask, times)
        @test length(clusters) >= 2
    end

    # ────────────────────────────────────────────────────────────
    # 4. _compute_cluster_statistics
    # ────────────────────────────────────────────────────────────
    @testset "_compute_cluster_statistics" begin
        # Create a simple cluster
        cluster = EegFun.TemporalCluster(1, [3, 4, 5], (0.2, 0.4), 0.0, 1.0, false)
        t_statistics = [0.0, 0.5, 2.0, 3.0, 1.5, 0.1, -0.5, 0.0, 0.0, 0.0]

        # Sum statistic
        stats_sum = EegFun._compute_cluster_statistics([cluster], t_statistics, :sum)
        @test length(stats_sum) == 1
        @test stats_sum[1] ≈ 2.0 + 3.0 + 1.5  # Sum of t-stats at indices 3, 4, 5

        # Max statistic
        stats_max = EegFun._compute_cluster_statistics([cluster], t_statistics, :max)
        @test length(stats_max) == 1
        @test stats_max[1] ≈ 3.0  # Max abs t-stat at indices 3, 4, 5

        # Empty clusters
        empty_stats = EegFun._compute_cluster_statistics(EegFun.TemporalCluster[], t_statistics)
        @test isempty(empty_stats)
    end

end # @testset "Decoding Statistics"
