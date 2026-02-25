using Test
using Random
using SparseArrays
using CairoMakie

# ============================================================
# Helper: build a Vector{TimeFreqData} with N participants × 2 conditions
# ============================================================
function _build_tf_test_data(;
    n_participants = 6,
    n_channels = 3,
    frequencies = [4.0, 8.0, 12.0, 16.0],
    time_points = [0.0, 0.1, 0.2, 0.3, 0.4],
    offset_cond1 = 5.0,
    offset_cond2 = 3.0,
    noise = 0.05,
    fs = 256,
)
    tfs = EegFun.TimeFreqData[]
    for p = 1:n_participants
        for (cond, offset) in [(1, offset_cond1), (2, offset_cond2)]
            tf = EegFun.create_test_tf_data(
                file = "participant$p",
                condition = cond,
                condition_name = "cond_$cond",
                n_channels = n_channels,
                frequencies = frequencies,
                time_points = time_points,
                power_offset = offset,
                noise = noise,
                fs = fs,
            )
            push!(tfs, tf)
        end
    end
    return tfs
end

# ============================================================
@testset "TF Statistics" begin
    # ============================================================

    # ────────────────────────────────────────────────────────────
    # 1. prepare_stats (TF)
    # ────────────────────────────────────────────────────────────
    @testset "prepare_stats (TF)" begin
        Random.seed!(42)
        n_participants = 6
        n_channels = 3
        freqs = [4.0, 8.0, 12.0, 16.0]
        times = [0.0, 0.1, 0.2, 0.3, 0.4]

        tfs = _build_tf_test_data(n_participants = n_participants, n_channels = n_channels, frequencies = freqs, time_points = times)

        prepared = EegFun.prepare_stats(tfs; design = :paired)

        @test prepared isa EegFun.TFStatisticalData
        @test length(prepared.data) == 2
        @test prepared.data[1] isa EegFun.TimeFreqData
        @test prepared.data[2] isa EegFun.TimeFreqData

        # 4D data arrays: [participants × electrodes × frequencies × time]
        @test size(prepared.analysis.data[1]) == (n_participants, n_channels, length(freqs), length(times))
        @test size(prepared.analysis.data[2]) == size(prepared.analysis.data[1])

        @test prepared.analysis.design == :paired
        @test prepared.analysis.frequencies == freqs
        @test prepared.analysis.time_points == times

        # Condition 1 mean should be ~5.0, condition 2 ~3.0
        mean1 = sum(prepared.analysis.data[1]) / length(prepared.analysis.data[1])
        mean2 = sum(prepared.analysis.data[2]) / length(prepared.analysis.data[2])
        @test mean1 > mean2
        @test abs(mean1 - 5.0) < 0.5
        @test abs(mean2 - 3.0) < 0.5

        # Grand averages should reflect condition names
        @test prepared.data[1].condition_name == "cond_1"
        @test prepared.data[2].condition_name == "cond_2"
    end

    @testset "prepare_stats (TF) - frequency selection" begin
        Random.seed!(43)
        tfs = _build_tf_test_data(frequencies = [2.0, 4.0, 8.0, 12.0, 20.0, 30.0])

        prepared = EegFun.prepare_stats(tfs; design = :paired, frequency_selection = (4.0, 20.0))

        @test prepared.analysis.frequencies == [4.0, 8.0, 12.0, 20.0]
        @test size(prepared.analysis.data[1], 3) == 4  # 4 frequencies selected
    end

    @testset "prepare_stats (TF) - time selection" begin
        Random.seed!(44)
        tfs = _build_tf_test_data(time_points = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

        prepared = EegFun.prepare_stats(tfs; design = :paired, interval_selection = (0.0, 0.3))

        @test prepared.analysis.time_points == [0.0, 0.1, 0.2, 0.3]
        @test size(prepared.analysis.data[1], 4) == 4  # 4 time points
    end

    @testset "prepare_stats (TF) - requires 2 conditions" begin
        Random.seed!(45)
        # Build data with only 1 condition
        tfs_one = [EegFun.create_test_tf_data(file = "p$i", condition = 1, condition_name = "cond_1") for i = 1:4]
        @test_throws Exception EegFun.prepare_stats(tfs_one; design = :paired)
    end


    # ────────────────────────────────────────────────────────────
    # 2. _compute_t_matrix_tf
    # ────────────────────────────────────────────────────────────
    @testset "_compute_t_matrix_tf" begin
        Random.seed!(46)
        n_participants = 8
        n_electrodes = 2
        n_freqs = 3
        n_time = 4

        # Condition A: constant 5.0, Condition B: constant 3.0 → large positive t
        data1 = fill(5.0, n_participants, n_electrodes, n_freqs, n_time) .+ 0.01 .* randn(n_participants, n_electrodes, n_freqs, n_time)
        data2 = fill(3.0, n_participants, n_electrodes, n_freqs, n_time) .+ 0.01 .* randn(n_participants, n_electrodes, n_freqs, n_time)

        t_matrix, df, p_matrix = EegFun._compute_t_matrix_tf(data1, data2, :paired)

        # Dimensions: [electrodes × freqs × time]
        @test size(t_matrix) == (n_electrodes, n_freqs, n_time)
        @test size(p_matrix) == (n_electrodes, n_freqs, n_time)

        # df for paired = n_participants - 1
        @test df == Float64(n_participants - 1)

        # All t-values should be large and positive (cond A > cond B)
        @test all(t_matrix .> 10.0)

        # All p-values should be very small
        @test all(p_matrix .< 0.001)
    end

    @testset "_compute_t_matrix_tf - no difference" begin
        Random.seed!(47)
        n = 10
        data_same = fill(5.0, n, 2, 3, 4)

        t_matrix, df, p_matrix = EegFun._compute_t_matrix_tf(data_same, data_same, :paired)

        # Zero difference → NaN t-values (0/0)
        @test all(isnan, t_matrix)
    end

    @testset "_compute_t_matrix_tf - prepared dispatch" begin
        Random.seed!(48)
        tfs = _build_tf_test_data(n_participants = 5)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        t_matrix, df, p_matrix = EegFun._compute_t_matrix_tf(prepared)

        @test size(t_matrix, 1) == 3  # n_channels
        @test size(t_matrix, 2) == 4  # n_freqs
        @test size(t_matrix, 3) == 5  # n_time
        @test df == 4.0  # 5 participants - 1
    end


    # ────────────────────────────────────────────────────────────
    # 3. Thresholding (TF)
    # ────────────────────────────────────────────────────────────
    @testset "parametric thresholding (TF)" begin
        t_matrix = zeros(2, 3, 4)
        t_matrix[1, 1, 1] = 3.0   # large positive
        t_matrix[2, 2, 3] = -3.0  # large negative
        t_matrix[1, 3, 4] = 1.5   # moderate positive

        critical_t = fill(2.0, 2, 3, 4)

        mask_pos, mask_neg = EegFun._threshold_t_matrix_parametric_tf(t_matrix, critical_t, :both)

        @test mask_pos[1, 1, 1] == true   # 3.0 > 2.0
        @test mask_neg[2, 2, 3] == true   # -3.0 < -2.0
        @test mask_pos[1, 3, 4] == false  # 1.5 < 2.0
        @test count(mask_pos) == 1
        @test count(mask_neg) == 1
    end

    @testset "parametric thresholding (TF) - right tail" begin
        t_matrix = zeros(2, 2, 2)
        t_matrix[1, 1, 1] = 3.0
        t_matrix[2, 2, 2] = -3.0

        critical_t = fill(2.0, 2, 2, 2)
        mask_pos, mask_neg = EegFun._threshold_t_matrix_parametric_tf(t_matrix, critical_t, :right)

        @test mask_pos[1, 1, 1] == true
        @test mask_neg[2, 2, 2] == false  # ignored for right tail
        @test count(mask_neg) == 0
    end

    @testset "parametric thresholding (TF) - in-place" begin
        t_matrix = zeros(2, 2, 3)
        t_matrix[1, 1, 1] = 5.0
        t_matrix[2, 2, 3] = -5.0

        critical_t = fill(2.0, 2, 2, 3)
        mask_pos = falses(2, 2, 3)
        mask_neg = falses(2, 2, 3)

        EegFun._threshold_t_matrix_parametric_tf!(mask_pos, mask_neg, t_matrix, critical_t, :both)

        @test mask_pos[1, 1, 1] == true
        @test mask_neg[2, 2, 3] == true
    end


    # ────────────────────────────────────────────────────────────
    # 4. 3D Clustering
    # ────────────────────────────────────────────────────────────
    @testset "3D BFS clustering - temporal" begin
        # 2 electrodes × 3 frequencies × 5 time points
        mask = falses(2, 3, 5)
        # Create a cluster connected in time at electrode 1, freq 2
        mask[1, 2, 2] = true
        mask[1, 2, 3] = true
        mask[1, 2, 4] = true
        # Create an isolated point
        mask[2, 1, 1] = true

        electrodes = [:Ch1, :Ch2]
        frequencies = [4.0, 8.0, 12.0]
        time_points = [0.0, 0.1, 0.2, 0.3, 0.4]
        spatial_conn = sparse(Int[], Int[], Bool[], 2, 2)

        clusters = EegFun._find_clusters_connected_components_tf(mask, electrodes, frequencies, time_points, spatial_conn, :temporal)

        @test length(clusters) == 2  # two separate clusters

        # Find the larger cluster (3 pixels)
        big = findfirst(c -> length(c.pixels) == 3, clusters)
        @test big |> !isnothing
        @test clusters[big].freq_range == (8.0, 8.0)
        @test length(clusters[big].time_indices) == 3

        # The isolated point is its own cluster
        small = findfirst(c -> length(c.pixels) == 1, clusters)
        @test small |> !isnothing
    end

    @testset "3D BFS clustering - spectral" begin
        mask = falses(1, 4, 1)
        mask[1, 1, 1] = true
        mask[1, 2, 1] = true
        mask[1, 3, 1] = true
        # gap at f=4 → separate cluster
        mask[1, 4, 1] = false

        electrodes = [:Ch1]
        frequencies = [4.0, 8.0, 12.0, 16.0]
        time_points = [0.0]
        spatial_conn = sparse(Int[], Int[], Bool[], 1, 1)

        clusters = EegFun._find_clusters_connected_components_tf(mask, electrodes, frequencies, time_points, spatial_conn, :spectral)

        @test length(clusters) == 1
        @test clusters[1].freq_range == (4.0, 12.0)
    end

    @testset "3D BFS clustering - full connectivity" begin
        # Two electrodes spatially connected, spanning freq and time
        mask = falses(2, 2, 2)
        mask[1, 1, 1] = true
        mask[2, 1, 1] = true  # spatial neighbour
        mask[1, 2, 1] = true  # spectral neighbour
        mask[1, 1, 2] = true  # temporal neighbour

        electrodes = [:Ch1, :Ch2]
        frequencies = [4.0, 8.0]
        time_points = [0.0, 0.1]
        # Ch1 and Ch2 are spatial neighbours
        spatial_conn = sparse([1, 2], [2, 1], true, 2, 2)

        clusters = EegFun._find_clusters_connected_components_tf(mask, electrodes, frequencies, time_points, spatial_conn, :full)

        # All 4 points should form one connected cluster
        @test length(clusters) == 1
        @test length(clusters[1].pixels) == 4
    end

    @testset "_find_clusters_tf (positive + negative)" begin
        mask_pos = falses(2, 2, 3)
        mask_neg = falses(2, 2, 3)
        mask_pos[1, 1, 1] = true
        mask_pos[1, 1, 2] = true
        mask_neg[2, 2, 3] = true

        electrodes = [:Ch1, :Ch2]
        frequencies = [4.0, 8.0]
        time_points = [0.0, 0.1, 0.2]
        spatial_conn = sparse(Int[], Int[], Bool[], 2, 2)

        pos_clusters, neg_clusters =
            EegFun._find_clusters_tf(mask_pos, mask_neg, electrodes, frequencies, time_points, spatial_conn, :temporal)

        @test length(pos_clusters) == 1
        @test pos_clusters[1].polarity == :positive
        @test length(neg_clusters) == 1
        @test neg_clusters[1].polarity == :negative
    end

    @testset "_compute_cluster_statistics_tf" begin
        # Create a simple t_matrix and cluster
        t_matrix = ones(2, 2, 3) .* 2.0  # all t-values = 2.0

        cluster = EegFun.TFCluster(
            1,
            [:Ch1],
            [1],
            [1, 2, 3],
            (4.0, 4.0),
            (0.0, 0.2),
            0.0,
            1.0,
            false,
            :positive,
            [CartesianIndex(1, 1, 1), CartesianIndex(1, 1, 2), CartesianIndex(1, 1, 3)],
        )

        updated, stats = EegFun._compute_cluster_statistics_tf([cluster], t_matrix, [:Ch1, :Ch2])

        @test length(stats) == 1
        @test stats[1] ≈ 6.0  # 3 pixels × 2.0
        @test updated[1].cluster_stat ≈ 6.0
    end

    @testset "_prefilter_mask_by_neighbors_tf!" begin
        # 3 electrodes, 2 freqs, 2 time points
        mask = trues(3, 2, 2)

        # Ch1-Ch2 connected, Ch3 isolated
        spatial_conn = sparse([1, 2], [2, 1], true, 3, 3)

        EegFun._prefilter_mask_by_neighbors_tf!(mask, spatial_conn, 2)

        # Ch3 has no spatial neighbours → removed (self count = 1, need >= 2)
        @test all(mask[3, :, :] .== false)
        # Ch1 and Ch2 should survive (each has 1 neighbour + self = 2)
        @test all(mask[1, :, :] .== true)
        @test all(mask[2, :, :] .== true)
    end

    @testset "empty mask produces no clusters" begin
        mask = falses(2, 3, 4)
        electrodes = [:Ch1, :Ch2]
        frequencies = [4.0, 8.0, 12.0]
        time_points = [0.0, 0.1, 0.2, 0.3]
        spatial_conn = sparse(Int[], Int[], Bool[], 2, 2)

        clusters = EegFun._find_clusters_connected_components_tf(mask, electrodes, frequencies, time_points, spatial_conn, :full)
        @test isempty(clusters)
    end


    # ────────────────────────────────────────────────────────────
    # 5. analytic_test (TF)
    # ────────────────────────────────────────────────────────────
    @testset "analytic_test (TF)" begin
        Random.seed!(50)
        tfs = _build_tf_test_data(n_participants = 8, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.analytic_test(prepared; alpha = 0.05)

        @test result isa EegFun.TFAnalyticResult
        @test result.test_info.type == :paired
        @test result.test_info.alpha == 0.05
        @test result.test_info.tail == :both
        @test result.test_info.correction_method == :no

        # Dimensions
        @test size(result.stat_matrix.t) == (3, 4, 5)  # n_channels × n_freqs × n_time
        @test size(result.stat_matrix.p) == (3, 4, 5)
        @test size(result.masks.positive) == (3, 4, 5)
        @test size(result.masks.negative) == (3, 4, 5)

        # Strong effect: cond1 (5.0) > cond2 (3.0) → positive significance everywhere
        @test count(result.masks.positive) > 0
        @test count(result.masks.negative) == 0  # no negative effect

        # Critical t should be positive
        @test result.critical_t > 0

        # Electrodes, frequencies, time_points
        @test length(result.electrodes) == 3
        @test result.frequencies == [4.0, 8.0, 12.0, 16.0]
        @test result.time_points == [0.0, 0.1, 0.2, 0.3, 0.4]
    end

    @testset "analytic_test (TF) - Bonferroni" begin
        Random.seed!(51)
        tfs = _build_tf_test_data(n_participants = 8, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result_no = EegFun.analytic_test(prepared; alpha = 0.05, correction_method = :no)
        result_bf = EegFun.analytic_test(prepared; alpha = 0.05, correction_method = :bonferroni)

        # Bonferroni should have equal or fewer significant points
        @test count(result_bf.masks.positive) <= count(result_no.masks.positive)
    end

    @testset "analytic_test (TF) - right tail" begin
        Random.seed!(52)
        tfs = _build_tf_test_data(n_participants = 8, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.analytic_test(prepared; tail = :right)

        # Right tail: can have positive significance but no negative
        @test count(result.masks.negative) == 0
    end

    @testset "analytic_test (TF) - left tail" begin
        Random.seed!(53)
        tfs = _build_tf_test_data(n_participants = 8, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.analytic_test(prepared; tail = :left)

        # Left tail with cond1 > cond2: no significance expected
        @test count(result.masks.positive) == 0
    end


    # ────────────────────────────────────────────────────────────
    # 6. permutation_test (TF)
    # ────────────────────────────────────────────────────────────
    @testset "permutation_test (TF) - basic" begin
        Random.seed!(60)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 50, cluster_type = :temporal)

        @test result isa EegFun.TFClusterPermutationResult
        @test result.test_info.type == :paired
        @test result.test_info.cluster_info.cluster_type == :temporal
        @test result.test_info.cluster_info.n_permutations == 50

        # Check dimensions
        @test size(result.stat_matrix.t) == (3, 4, 5)
        @test size(result.masks.positive) == (3, 4, 5)
        @test size(result.masks.negative) == (3, 4, 5)

        # Permutation distribution should have 50 entries
        @test length(result.permutation_distribution.positive) == 50
        @test length(result.permutation_distribution.negative) == 50

        # Result metadata
        @test length(result.electrodes) == 3
        @test result.frequencies == [4.0, 8.0, 12.0, 16.0]
        @test result.time_points == [0.0, 0.1, 0.2, 0.3, 0.4]

        # Clusters should be TFCluster vectors
        @test result.clusters.positive isa Vector{EegFun.TFCluster}
        @test result.clusters.negative isa Vector{EegFun.TFCluster}
    end

    @testset "permutation_test (TF) - full connectivity" begin
        Random.seed!(61)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 30, cluster_type = :full)

        @test result isa EegFun.TFClusterPermutationResult
        @test result.test_info.cluster_info.cluster_type == :full
    end

    @testset "permutation_test (TF) - no effect" begin
        Random.seed!(62)
        # Both conditions have same offset → no real effect
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 5.0, noise = 0.5)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 30, cluster_type = :temporal)

        @test result isa EegFun.TFClusterPermutationResult
        # With no effect, significant clusters should be rare/nonexistent
        n_sig = count(c -> c.is_significant, result.clusters.positive) + count(c -> c.is_significant, result.clusters.negative)
        @test n_sig == 0
    end


    # ────────────────────────────────────────────────────────────
    # 7. Type show methods
    # ────────────────────────────────────────────────────────────
    @testset "TFStatisticalData show" begin
        Random.seed!(70)
        tfs = _build_tf_test_data(n_participants = 4)
        prepared = EegFun.prepare_stats(tfs; design = :paired)

        io = IOBuffer()
        show(io, prepared)
        output = String(take!(io))

        @test contains(output, "TFStatisticalData")
        @test contains(output, "paired")
    end

    @testset "TFClusterPermutationResult show" begin
        Random.seed!(71)
        tfs = _build_tf_test_data(n_participants = 4, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.permutation_test(prepared; n_permutations = 20, cluster_type = :temporal)

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "TFClusterPermutationResult")
        @test contains(output, "paired")
    end

    @testset "TFAnalyticResult show" begin
        Random.seed!(72)
        tfs = _build_tf_test_data(n_participants = 4, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "TFAnalyticResult")
    end


    # ────────────────────────────────────────────────────────────
    # 8. nonparametric thresholding (TF)
    # ────────────────────────────────────────────────────────────
    @testset "nonparametric threshold common (TF)" begin
        Random.seed!(80)
        # 2 electrodes × 3 freqs × 4 time × 100 permutations
        perm_t = randn(2, 3, 4, 100)

        thresh_pos, thresh_neg = EegFun._compute_nonparametric_threshold_common_tf(perm_t, 0.05, :both)

        @test thresh_pos > 0
        @test thresh_neg > 0
        @test thresh_pos == thresh_neg  # symmetric for :both
    end

    @testset "nonparametric thresholding in-place (TF)" begin
        t_matrix = zeros(2, 2, 2)
        t_matrix[1, 1, 1] = 5.0
        t_matrix[2, 2, 2] = -5.0

        mask_pos = falses(2, 2, 2)
        mask_neg = falses(2, 2, 2)

        EegFun._threshold_t_matrix_nonparametric_tf!(mask_pos, mask_neg, t_matrix, 2.0, 2.0, :both)

        @test mask_pos[1, 1, 1] == true   # 5.0 > 2.0
        @test mask_neg[2, 2, 2] == true   # -5.0 < -2.0
        @test count(mask_pos) == 1
        @test count(mask_neg) == 1
    end


    # ────────────────────────────────────────────────────────────
    # 9. Critical t-values (TF)
    # ────────────────────────────────────────────────────────────
    @testset "critical t-values (TF)" begin
        # With df=9, alpha=0.05, two-tailed: critical t ≈ 2.262
        crit = EegFun._compute_critical_t_values_tf(9.0, (2, 3, 4), 0.05, :both)

        @test size(crit) == (2, 3, 4)
        @test all(crit .≈ crit[1, 1, 1])  # uniform
        @test 2.2 < crit[1, 1, 1] < 2.4   # ~2.262

        # Invalid df
        crit_nan = EegFun._compute_critical_t_values_tf(NaN, (2, 3, 4))
        @test all(isnan, crit_nan)
    end


    # ────────────────────────────────────────────────────────────
    # 10. plot_tf_stats
    # ────────────────────────────────────────────────────────────
    @testset "plot_tf_stats - analytic result" begin
        Random.seed!(90)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        # Default: t-values with contour
        out = EegFun.plot_tf_stats(result; display_plot = false)
        @test out.fig isa Figure
        @test length(out.axes) == 3  # 3 channels

        # Single channel
        out2 = EegFun.plot_tf_stats(result; channel_selection = EegFun.channels(:Ch1), display_plot = false)
        @test length(out2.axes) == 1
    end

    @testset "plot_tf_stats - permutation result" begin
        Random.seed!(91)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.permutation_test(prepared; n_permutations = 20, cluster_type = :temporal)

        out = EegFun.plot_tf_stats(result; display_plot = false)
        @test out.fig isa Figure
        @test length(out.axes) == 3
    end

    @testset "plot_tf_stats - content modes" begin
        Random.seed!(92)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        for mode in [:tvalues, :difference, :power_a, :power_b]
            out = EegFun.plot_tf_stats(result; content = mode, display_plot = false)
            @test out.fig isa Figure
        end
    end

    @testset "plot_tf_stats - significance modes" begin
        Random.seed!(93)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        for mode in [:contour, :stipple, :opacity, :none]
            out = EegFun.plot_tf_stats(result; significance = mode, display_plot = false)
            @test out.fig isa Figure
        end
    end

    @testset "plot_tf_stats - invalid arguments" begin
        Random.seed!(94)
        tfs = _build_tf_test_data(n_participants = 6)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        @test_throws Exception EegFun.plot_tf_stats(result; content = :invalid, display_plot = false)
        @test_throws Exception EegFun.plot_tf_stats(result; significance = :invalid, display_plot = false)
    end

    @testset "plot_tf_stats - internal helpers" begin
        Random.seed!(95)
        tfs = _build_tf_test_data(n_participants = 6, offset_cond1 = 5.0, offset_cond2 = 3.0, noise = 0.05)
        prepared = EegFun.prepare_stats(tfs; design = :paired)
        result = EegFun.analytic_test(prepared)

        # _extract_tf_stats_matrix
        mat = EegFun._extract_tf_stats_matrix(result, 1, :tvalues)
        @test size(mat) == (4, 5)  # n_freqs × n_time

        mat_diff = EegFun._extract_tf_stats_matrix(result, 1, :difference)
        @test size(mat_diff) == (4, 5)
        # Difference should be positive (cond A > cond B)
        @test sum(mat_diff) > 0

        # _extract_tf_significance_mask
        sig = EegFun._extract_tf_significance_mask(result, 1)
        @test sig isa BitMatrix
        @test size(sig) == (4, 5)
    end


    # ────────────────────────────────────────────────────────────
    # 11. Auto-detect prepare_stats (file-based)
    # ────────────────────────────────────────────────────────────
    @testset "prepare_stats - auto-detect file-based" begin
        TF_TEST_DIR = joinpath(@__DIR__, "..", "..", "resources", "data", "julia", "tf")
        prepared = EegFun.prepare_stats("tf_data", :paired; input_dir = TF_TEST_DIR)
        @test prepared isa EegFun.TFStatisticalData
        @test length(prepared.data) == 2  # 2 conditions
        @test size(prepared.analysis.data[1], 1) == 2  # 2 participants
    end


end # @testset "TF Statistics"
