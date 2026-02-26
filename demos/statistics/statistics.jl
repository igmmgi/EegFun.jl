# Demo: Permutation-Based Statistical Testing for ERPs
# Shows how to run permutation tests on ERP data using different
# thresholding and cluster correction approaches.

using EegFun

input_dir = "./resources/data/julia/erps"
file_pattern = "erps_good"

println("Preparing data...")
stat_data = EegFun.prepare_stats(
    file_pattern,  # Pattern to match ERP files
    :paired;       # :paired for within-subject, :independent for between-subject
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]), # Conditions to compare
    channel_selection = EegFun.channels(1:72),       # Select all 72 channels
    interval_selection = EegFun.times(),  # Full time interval
    baseline_interval = EegFun.times((-0.2, 0.0)),   # Baseline: -200 to 0 ms
    analysis_interval = EegFun.times((0.0, 1.0)),    # Analysis interval: 100-1000 ms
)

# ----------------------------------------------------------------------------
# Option: Analytic t-test with NO correction
# ----------------------------------------------------------------------------
result_analytic_no = EegFun.analytic_test(stat_data)

EegFun.plot_erp_stats(
    result_analytic_no,
    channel_selection = EegFun.channels(:PO8),
    plot_erp = true,
    plot_difference = false,
    plot_significance = true,
    plot_critical_t = true,
)

# ----------------------------------------------------------------------------
# Option: Analytic t-test with BONFERRONI correction (v. conservative; really only useful for demonstrating)
# ----------------------------------------------------------------------------
result_analytic_bonf = EegFun.analytic_test(stat_data, correction_method = :bonferroni)

EegFun.plot_erp_stats(
    result_analytic_bonf,
    channel_selection = EegFun.channels(:PO8),
    plot_erp = true,
    plot_difference = false,
    plot_significance = true,
    plot_critical_t = true,
)

# ----------------------------------------------------------------------------
# Option: Parametric Thresholding (Default, Fastest)
# ----------------------------------------------------------------------------
result_permutation_parametric = EegFun.permutation_test(
    stat_data,
    n_permutations = 1000,           # Number of permutations (more = more accurate)
    threshold_method = :parametric,  # Use t-distribution for thresholding
    cluster_type = :spatiotemporal,  # Cluster in space AND time
    min_num_neighbors = 3,           # Pre-filter: need 3+ neighbors
    show_progress = true,            # Show progress bar
)

EegFun.plot_erp_stats(
    result_permutation_parametric,
    channel_selection = EegFun.channels(:PO8),
    plot_erp = true,
    plot_difference = false,
    plot_significance = true,
    plot_critical_t = true,
)


# ----------------------------------------------------------------------------
# Option: Non-Parametric Common Thresholding
# ----------------------------------------------------------------------------
# What actually happens:
# Run ALL permutations first (collect t-matrices)
# Pool all t-values from all permutations
# Compute (1-α) percentile threshold from pooled distribution
# Use this single threshold for all points
# Reuse stored t-matrices for cluster-level inference
# Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster', clusterThreshold='nonparametric_common'

result_permutation_nonparametric_common = EegFun.permutation_test(
    stat_data,
    n_permutations = 1000,
    threshold_method = :nonparametric_common,  # Non-parametric common threshold
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
)

# ----------------------------------------------------------------------------
# Option: Non-Parametric Individual Thresholding
# ----------------------------------------------------------------------------
# What actually happens:
# Run ALL permutations first (collect t-matrices)
# For each electrode × time point:
#    - Extract t-values from all permutations at that point
#    - Compute (1-α) percentile threshold from point-specific distribution
# Threshold observed data using point-specific thresholds
# Reuse stored t-matrices for cluster-level inference
# Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster', clusterThreshold='nonparametric_individual'

result_permutation_nonparametric_individual = EegFun.permutation_test(
    stat_data,
    n_permutations = 1000,
    threshold_method = :nonparametric_individual,  # Non-parametric individual thresholds
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
)


