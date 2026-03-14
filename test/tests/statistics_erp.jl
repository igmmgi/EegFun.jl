using Test
using Random
using SparseArrays

# ============================================================
# Helper: build a Vector{ErpData} with N participants × 2 conditions
# Each participant gets a unique file name so prepare_stats can match them
# ============================================================
function _build_erp_test_data(;
    n_participants = 6,
    n_channels = 3,
    fs = 256,
    offset_cond1 = 2.0,
    offset_cond2 = 0.0,
    noise = 0.05,
    signal_freq = 10.0,
)
    erps = EegFun.ErpData[]
    time = collect(-0.2:(1/fs):0.8)
    n_time = length(time)

    for p = 1:n_participants
        for (cond, offset) in [(1, offset_cond1), (2, offset_cond2)]
            df = DataFrame(
                time = time,
                sample = 1:n_time,
                condition = fill(cond, n_time),
                condition_name = fill("condition_$cond", n_time),
                participant = fill(p, n_time),
            )

            channel_labels = [Symbol("Ch$i") for i = 1:n_channels]
            Random.seed!(100 * p + 10 * cond)
            for ch in channel_labels
                # Each channel gets a sinusoidal signal with condition-dependent amplitude + noise
                # This survives baseline correction because the signal varies over time
                df[!, ch] = offset .* sin.(2π * signal_freq .* time) .+ noise .* randn(n_time)
            end

            layout = EegFun.create_test_layout(n_channels = n_channels)
            erp = EegFun.ErpData("participant$p", cond, "condition_$cond", df, layout, fs, EegFun.AnalysisInfo(:none, 0.0, 0.0), 20)
            push!(erps, erp)
        end
    end
    return erps
end

# ============================================================
@testset "ERP Statistics" begin
    # ============================================================

    # ────────────────────────────────────────────────────────────
    # 1. prepare_stats (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "prepare_stats" begin
        Random.seed!(42)
        n_participants = 6
        n_channels = 3

        erps = _build_erp_test_data(n_participants = n_participants, n_channels = n_channels)

        prepared = EegFun.prepare_stats(erps; design = :paired)

        @test prepared isa EegFun.StatisticalData
        @test length(prepared.data) == 2
        @test prepared.data[1] isa EegFun.ErpData
        @test prepared.data[2] isa EegFun.ErpData

        # 3D data arrays: [participants × electrodes × time]
        @test size(prepared.analysis.data[1], 1) == n_participants
        @test size(prepared.analysis.data[1], 2) == n_channels
        @test size(prepared.analysis.data[2]) == size(prepared.analysis.data[1])

        @test prepared.analysis.design == :paired

        # Condition 1 should have larger amplitude (more variance) than condition 2
        var1 = sum(prepared.analysis.data[1] .^ 2) / length(prepared.analysis.data[1])
        var2 = sum(prepared.analysis.data[2] .^ 2) / length(prepared.analysis.data[2])
        @test var1 > var2

        # Grand averages should reflect condition names (prefixed with grand_avg_)
        @test contains(prepared.data[1].condition_name, "condition_1")
        @test contains(prepared.data[2].condition_name, "condition_2")

        # Per-condition SE should be computed over full display interval
        n_display_time = nrow(prepared.data[1].data)
        @test size(prepared.se_cond1) == (n_channels, n_display_time)
        @test size(prepared.se_cond2) == (n_channels, n_display_time)
        @test all(x -> x >= 0, prepared.se_cond1)
        @test all(x -> x >= 0, prepared.se_cond2)
    end

    @testset "prepare_stats - requires 2 conditions" begin
        Random.seed!(43)
        # Build data with only 1 condition
        erps_one = [EegFun.create_test_erp_data(participant = i, condition = 1) for i = 1:4]
        @test_throws Exception EegFun.prepare_stats(erps_one; design = :paired)
    end

    @testset "prepare_stats - channel selection" begin
        Random.seed!(44)
        erps = _build_erp_test_data(n_channels = 5)

        # Select only first 2 channels
        prepared = EegFun.prepare_stats(erps; design = :paired, channel_selection = EegFun.channels([:Ch1, :Ch2]))

        @test size(prepared.analysis.data[1], 2) == 2  # only 2 channels
    end


    # ────────────────────────────────────────────────────────────
    # 2. _compute_t_matrix (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "_compute_t_matrix - paired" begin
        Random.seed!(46)
        n_participants = 8
        n_electrodes = 2
        n_time = 5

        # Condition A: constant 5.0, Condition B: constant 3.0 → large positive t
        data1 = fill(5.0, n_participants, n_electrodes, n_time) .+ 0.01 .* randn(n_participants, n_electrodes, n_time)
        data2 = fill(3.0, n_participants, n_electrodes, n_time) .+ 0.01 .* randn(n_participants, n_electrodes, n_time)

        t_matrix, df, p_matrix, se_matrix = EegFun._compute_t_matrix(data1, data2, :paired)

        # Dimensions: [electrodes × time]
        @test size(t_matrix) == (n_electrodes, n_time)
        @test size(p_matrix) == (n_electrodes, n_time)
        @test size(se_matrix) == (n_electrodes, n_time)

        # df for paired = n_participants - 1
        @test df == Float64(n_participants - 1)

        # All t-values should be large and positive (cond A > cond B)
        @test all(t_matrix .> 10.0)

        # All p-values should be very small
        @test all(p_matrix .< 0.001)

        # SE should be positive
        @test all(se_matrix .> 0)
    end

    @testset "_compute_t_matrix - no difference" begin
        Random.seed!(47)
        n = 10
        data_same = fill(5.0, n, 3, 4)

        t_matrix, df, p_matrix, se_matrix = EegFun._compute_t_matrix(data_same, data_same, :paired)

        # Zero difference → NaN t-values (0/0)
        @test all(isnan, t_matrix)

        # SE should be zero (no variability in differences)
        @test all(se_matrix .== 0.0)
    end

    @testset "_compute_t_matrix - SE correctness" begin
        # Known data: 4 participants, 1 electrode, 1 time point
        # Differences: [2, 4, 6, 8] → mean=5, std=sqrt(20/3), se=std/sqrt(4)
        using Statistics: std as std_func
        data1 = reshape([12.0, 14.0, 16.0, 18.0], 4, 1, 1)
        data2 = reshape([10.0, 10.0, 10.0, 10.0], 4, 1, 1)

        t_matrix, df, p_matrix, se_matrix = EegFun._compute_t_matrix(data1, data2, :paired)

        diffs = [2.0, 4.0, 6.0, 8.0]
        expected_se = std_func(diffs) / sqrt(4)
        @test se_matrix[1, 1] ≈ expected_se
        @test t_matrix[1, 1] ≈ mean(diffs) / expected_se
    end

    @testset "_compute_t_matrix - prepared dispatch" begin
        Random.seed!(48)
        erps = _build_erp_test_data(n_participants = 5, n_channels = 3)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        t_matrix, df, p_matrix, se_matrix = EegFun._compute_t_matrix(prepared)

        @test size(t_matrix, 1) == 3   # n_channels
        @test df == 4.0                # 5 participants - 1
        @test size(t_matrix) == size(p_matrix)
        @test size(se_matrix) == size(t_matrix)
    end


    # ────────────────────────────────────────────────────────────
    # 3. Thresholding (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "parametric thresholding" begin
        t_matrix = zeros(3, 5)
        t_matrix[1, 1] = 3.0    # large positive
        t_matrix[2, 3] = -3.0   # large negative
        t_matrix[3, 5] = 1.5    # moderate positive

        critical_t = fill(2.0, 3, 5)

        mask_pos, mask_neg = EegFun._threshold_t_matrix_parametric(t_matrix, critical_t, :both)

        @test mask_pos[1, 1] == true    # 3.0 > 2.0
        @test mask_neg[2, 3] == true    # -3.0 < -2.0
        @test mask_pos[3, 5] == false   # 1.5 < 2.0
        @test count(mask_pos) == 1
        @test count(mask_neg) == 1
    end

    @testset "parametric thresholding - right tail" begin
        t_matrix = zeros(2, 2)
        t_matrix[1, 1] = 3.0
        t_matrix[2, 2] = -3.0

        critical_t = fill(2.0, 2, 2)
        mask_pos, mask_neg = EegFun._threshold_t_matrix_parametric(t_matrix, critical_t, :right)

        @test mask_pos[1, 1] == true
        @test mask_neg[2, 2] == false   # ignored for right tail
        @test count(mask_neg) == 0
    end

    @testset "parametric thresholding - in-place" begin
        t_matrix = zeros(2, 3)
        t_matrix[1, 1] = 5.0
        t_matrix[2, 3] = -5.0

        critical_t = fill(2.0, 2, 3)
        mask_pos = falses(2, 3)
        mask_neg = falses(2, 3)

        EegFun._threshold_t_matrix_parametric!(mask_pos, mask_neg, t_matrix, critical_t, :both)

        @test mask_pos[1, 1] == true
        @test mask_neg[2, 3] == true
    end

    @testset "critical t-values" begin
        # With df=9, alpha=0.05, two-tailed: critical t ≈ 2.262
        crit = EegFun._compute_critical_t_values(9.0, (3, 5), 0.05, :both)

        @test size(crit) == (3, 5)
        @test all(crit .≈ crit[1, 1])   # uniform
        @test 2.2 < crit[1, 1] < 2.4    # ~2.262

        # Invalid df
        crit_nan = EegFun._compute_critical_t_values(NaN, (3, 5))
        @test all(isnan, crit_nan)
    end


    # ────────────────────────────────────────────────────────────
    # 4. Clustering (ERP - 2D)
    # ────────────────────────────────────────────────────────────
    @testset "2D BFS clustering - temporal" begin
        # 3 electrodes × 5 time points
        mask = falses(3, 5)
        # Create a cluster connected in time at electrode 1
        mask[1, 2] = true
        mask[1, 3] = true
        mask[1, 4] = true
        # Create an isolated point
        mask[3, 1] = true

        electrodes = [:Ch1, :Ch2, :Ch3]
        time_points = [0.0, 0.1, 0.2, 0.3, 0.4]
        spatial_conn = sparse(Int[], Int[], Bool[], 3, 3)

        clusters = EegFun._find_clusters_connected_components(mask, electrodes, time_points, spatial_conn, :temporal)

        @test length(clusters) == 2  # two separate clusters

        # Find the larger cluster (3 points)
        big = findfirst(c -> length(c.time_indices) == 3, clusters)
        @test !isnothing(big)

        # The isolated point is its own cluster
        small = findfirst(c -> length(c.time_indices) == 1, clusters)
        @test !isnothing(small)
    end

    @testset "2D BFS clustering - spatiotemporal" begin
        # 2 electrodes × 3 time points
        mask = falses(2, 3)
        mask[1, 1] = true
        mask[2, 1] = true   # spatial neighbour at same time
        mask[1, 2] = true   # temporal neighbour

        electrodes = [:Ch1, :Ch2]
        time_points = [0.0, 0.1, 0.2]
        # Ch1 and Ch2 are spatial neighbours
        spatial_conn = sparse([1, 2], [2, 1], true, 2, 2)

        clusters = EegFun._find_clusters_connected_components(mask, electrodes, time_points, spatial_conn, :spatiotemporal)

        # All 3 points should form one connected cluster
        @test length(clusters) == 1
    end

    @testset "_find_clusters (positive + negative)" begin
        mask_pos = falses(2, 3)
        mask_neg = falses(2, 3)
        mask_pos[1, 1] = true
        mask_pos[1, 2] = true
        mask_neg[2, 3] = true

        electrodes = [:Ch1, :Ch2]
        time_points = [0.0, 0.1, 0.2]
        spatial_conn = sparse(Int[], Int[], Bool[], 2, 2)

        pos_clusters, neg_clusters = EegFun._find_clusters(mask_pos, mask_neg, electrodes, time_points, spatial_conn, :temporal)

        @test length(pos_clusters) == 1
        @test pos_clusters[1].polarity == :positive
        @test length(neg_clusters) == 1
        @test neg_clusters[1].polarity == :negative
    end

    @testset "_compute_cluster_statistics" begin
        # Create a simple t_matrix and cluster
        t_matrix = ones(2, 3) .* 2.0  # all t-values = 2.0

        cluster = EegFun.Cluster(
            1,                      # id
            [:Ch1],                 # electrodes
            [1, 2, 3],             # time_indices
            (0.0, 0.2),           # time_range
            0.0,                   # cluster_stat (placeholder)
            1.0,                   # p_value (placeholder)
            false,                 # is_significant
            :positive,             # polarity
            [(1, 1), (1, 2), (1, 3)],  # members: Ch1 at t=1,2,3
        )

        updated, stats = EegFun._compute_cluster_statistics([cluster], t_matrix, [:Ch1, :Ch2])

        @test length(stats) == 1
        @test stats[1] ≈ 6.0   # 3 member points × 2.0
        @test updated[1].cluster_stat ≈ 6.0
    end

    @testset "empty mask produces no clusters" begin
        mask = falses(2, 5)
        electrodes = [:Ch1, :Ch2]
        time_points = [0.0, 0.1, 0.2, 0.3, 0.4]
        spatial_conn = sparse(Int[], Int[], Bool[], 2, 2)

        clusters = EegFun._find_clusters_connected_components(mask, electrodes, time_points, spatial_conn, :temporal)
        @test isempty(clusters)
    end

    @testset "_prefilter_mask_by_neighbors" begin
        # 3 electrodes × 3 time points
        mask = trues(3, 3)

        # Ch1-Ch2 connected, Ch3 isolated
        spatial_conn = sparse([1, 2], [2, 1], true, 3, 3)

        filtered = EegFun._prefilter_mask_by_neighbors(mask, spatial_conn, 1)

        # Ch3 has no spatial neighbours → removed
        @test all(filtered[3, :] .== false)
        # Ch1 and Ch2 each have 1 neighbour → survive (>= 1)
        @test all(filtered[1, :] .== true)
        @test all(filtered[2, :] .== true)
    end


    # ────────────────────────────────────────────────────────────
    # 5. Inference helpers
    # ────────────────────────────────────────────────────────────
    @testset "_compute_cluster_pvalues" begin
        cluster = EegFun.Cluster(1, [:Ch1], [1, 2], (0.0, 0.1), 5.0, 1.0, false, :positive, [(1, 1), (1, 2)])
        # With permutation max values all < 5.0, p should be very small
        perm_max = [1.0, 2.0, 3.0, 4.0, 2.5]

        updated = EegFun._compute_cluster_pvalues([cluster], [5.0], perm_max, 5, 0.05)

        @test length(updated) == 1
        @test updated[1].p_value < 0.5
        @test updated[1].p_value > 0.0
    end

    @testset "_create_significance_mask" begin
        p_matrix = [0.01 0.06 NaN; 0.04 0.10 0.001]

        mask = EegFun._create_significance_mask(p_matrix, 0.05, :no)
        @test mask[1, 1] == true    # 0.01 < 0.05
        @test mask[1, 2] == false   # 0.06 > 0.05
        @test mask[1, 3] == false   # NaN
        @test mask[2, 1] == true    # 0.04 < 0.05
        @test mask[2, 2] == false   # 0.10 > 0.05
        @test mask[2, 3] == true    # 0.001 < 0.05
    end

    @testset "_apply_bonferroni_correction" begin
        p_values = [0.001, 0.01, 0.04, 0.06]
        t_stats = [3.0, 2.5, 2.0, 1.5]

        mask = EegFun._apply_bonferroni_correction(p_values, 0.05, t_stats)

        # corrected alpha = 0.05/4 = 0.0125
        @test mask[1] == true    # 0.001 < 0.0125
        @test mask[2] == true    # 0.01 < 0.0125
        @test mask[3] == false   # 0.04 > 0.0125
        @test mask[4] == false   # 0.06 > 0.0125
    end


    # ────────────────────────────────────────────────────────────
    # 6. Nonparametric thresholding (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "nonparametric threshold common" begin
        Random.seed!(80)
        # 3 electrodes × 5 time × 100 permutations
        perm_t = randn(3, 5, 100)

        thresh_pos, thresh_neg = EegFun._compute_nonparametric_threshold_common(perm_t, 0.05, :both)

        @test thresh_pos > 0
        @test thresh_neg > 0
        @test thresh_pos == thresh_neg  # symmetric for :both
    end

    @testset "nonparametric thresholding in-place" begin
        t_matrix = zeros(2, 3)
        t_matrix[1, 1] = 5.0
        t_matrix[2, 3] = -5.0

        mask_pos = falses(2, 3)
        mask_neg = falses(2, 3)

        EegFun._threshold_t_matrix_nonparametric!(mask_pos, mask_neg, t_matrix, 2.0, 2.0, :both)

        @test mask_pos[1, 1] == true    # 5.0 > 2.0
        @test mask_neg[2, 3] == true    # -5.0 < -2.0
        @test count(mask_pos) == 1
        @test count(mask_neg) == 1
    end


    # ────────────────────────────────────────────────────────────
    # 7. analytic_test (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "analytic_test" begin
        Random.seed!(50)
        erps = _build_erp_test_data(n_participants = 8, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.analytic_test(prepared; alpha = 0.05)

        @test result isa EegFun.AnalyticResult
        @test result.test_info.type == :paired
        @test result.test_info.alpha == 0.05
        @test result.test_info.tail == :both
        @test result.test_info.correction_method == :no

        # Dimensions
        @test size(result.stat_matrix.t, 1) == 3   # n_channels
        @test size(result.stat_matrix.p) == size(result.stat_matrix.t)
        @test size(result.masks.positive) == size(result.stat_matrix.t)
        @test size(result.masks.negative) == size(result.stat_matrix.t)

        # With sinusoidal signal at different amplitudes, there should be significant differences
        @test count(result.masks.positive) + count(result.masks.negative) > 0

        # Critical t should be positive
        @test result.critical_t > 0

        # SE should have correct dimensions and be positive
        @test size(result.se) == size(result.stat_matrix.t)
        @test all(x -> x > 0 || isnan(x), result.se)

        # Per-condition SE should cover full display interval
        n_display_time = nrow(result.data[1].data)
        @test size(result.se_cond1) == (3, n_display_time)
        @test size(result.se_cond2) == (3, n_display_time)

        # Electrodes and time_points
        @test length(result.electrodes) == 3
    end

    @testset "analytic_test - Bonferroni" begin
        Random.seed!(51)
        erps = _build_erp_test_data(n_participants = 8, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result_no = EegFun.analytic_test(prepared; alpha = 0.05, correction_method = :no)
        result_bf = EegFun.analytic_test(prepared; alpha = 0.05, correction_method = :bonferroni)

        # Bonferroni should have equal or fewer significant points
        @test count(result_bf.masks.positive) <= count(result_no.masks.positive)
    end

    @testset "analytic_test - right tail" begin
        Random.seed!(52)
        erps = _build_erp_test_data(n_participants = 8, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.analytic_test(prepared; tail = :right)

        # Right tail: can have positive significance but no negative
        @test count(result.masks.negative) == 0
    end

    @testset "analytic_test - left tail" begin
        Random.seed!(53)
        erps = _build_erp_test_data(n_participants = 8, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.analytic_test(prepared; tail = :left)

        # Left tail with cond1 > cond2: no significance expected
        @test count(result.masks.positive) == 0
    end

    @testset "analytic_test - invalid args" begin
        Random.seed!(54)
        erps = _build_erp_test_data(n_participants = 5)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        @test_throws Exception EegFun.analytic_test(prepared; correction_method = :invalid)
        @test_throws Exception EegFun.analytic_test(prepared; tail = :invalid)
    end


    # ────────────────────────────────────────────────────────────
    # 8. permutation_test (ERP)
    # ────────────────────────────────────────────────────────────
    @testset "permutation_test - basic" begin
        Random.seed!(60)
        erps = _build_erp_test_data(n_participants = 6, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 50, cluster_type = :temporal)

        @test result isa EegFun.PermutationResult
        @test result.test_info.type == :paired
        @test result.test_info.cluster_info.cluster_type == :temporal
        @test result.test_info.cluster_info.n_permutations == 50

        # Check dimensions
        @test size(result.stat_matrix.t, 1) == 3   # n_channels
        @test size(result.masks.positive) == size(result.stat_matrix.t)
        @test size(result.masks.negative) == size(result.stat_matrix.t)

        # Permutation distribution should have 50 entries
        @test length(result.permutation_distribution.positive) == 50
        @test length(result.permutation_distribution.negative) == 50

        # Result metadata
        @test length(result.electrodes) == 3

        # SE should have correct dimensions
        @test size(result.se) == size(result.stat_matrix.t)
        @test all(x -> x > 0 || isnan(x), result.se)

        # Clusters should be Cluster vectors
        @test result.clusters.positive isa Vector{EegFun.Cluster}
        @test result.clusters.negative isa Vector{EegFun.Cluster}
    end

    @testset "permutation_test - spatiotemporal" begin
        Random.seed!(61)
        erps = _build_erp_test_data(n_participants = 6, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 30, cluster_type = :spatiotemporal)

        @test result isa EegFun.PermutationResult
        @test result.test_info.cluster_info.cluster_type == :spatiotemporal
    end

    @testset "permutation_test - no effect" begin
        Random.seed!(62)
        # Both conditions have same offset → no real effect
        erps = _build_erp_test_data(n_participants = 6, offset_cond1 = 2.0, offset_cond2 = 2.0, noise = 0.5)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 30, cluster_type = :temporal)

        @test result isa EegFun.PermutationResult
        # With no effect, significant clusters should be rare/nonexistent
        n_sig = count(c -> c.is_significant, result.clusters.positive) + count(c -> c.is_significant, result.clusters.negative)
        @test n_sig == 0
    end

    @testset "permutation_test - nonparametric common" begin
        Random.seed!(63)
        erps = _build_erp_test_data(n_participants = 6, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        result = EegFun.permutation_test(prepared; n_permutations = 30, threshold_method = :nonparametric_common, cluster_type = :temporal)

        @test result isa EegFun.PermutationResult
        @test result.test_info.cluster_info.threshold_method == :nonparametric_common
    end


    # ────────────────────────────────────────────────────────────
    # 9. Type show methods
    # ────────────────────────────────────────────────────────────
    @testset "StatisticalData show" begin
        Random.seed!(70)
        erps = _build_erp_test_data(n_participants = 4)
        prepared = EegFun.prepare_stats(erps; design = :paired)

        io = IOBuffer()
        show(io, prepared)
        output = String(take!(io))

        @test contains(output, "StatisticalData")
        @test contains(output, "paired")
    end

    @testset "PermutationResult show" begin
        Random.seed!(71)
        erps = _build_erp_test_data(n_participants = 4, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)
        result = EegFun.permutation_test(prepared; n_permutations = 20, cluster_type = :temporal)

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "PermutationResult")
        @test contains(output, "paired")
    end

    @testset "AnalyticResult show" begin
        Random.seed!(72)
        erps = _build_erp_test_data(n_participants = 4, offset_cond1 = 2.0, offset_cond2 = 0.0, noise = 0.05)
        prepared = EegFun.prepare_stats(erps; design = :paired)
        result = EegFun.analytic_test(prepared)

        io = IOBuffer()
        show(io, result)
        output = String(take!(io))

        @test contains(output, "AnalyticResult")
    end

end # @testset "ERP Statistics"
