# Statistics

This demo shows permutation-based statistical testing for ERP data, with multiple thresholding and cluster correction approaches.

## Key Functions

| Function | Purpose |
| --- | --- |
| `prepare_stats` | Load multi-participant ERPs and prepare for statistical comparison |
| `analytic_test` | Parametric t-test (uncorrected or Bonferroni) |
| `permutation_test` | Cluster-based permutation test with multiple thresholding options |
| `plot_erp_stats` | Visualise results with significance shading and critical t-values |

## Thresholding Methods

| Method | Description |
| --- | --- |
| `:parametric` | T-distribution threshold (fastest, default) |
| `:nonparametric_common` | Single threshold from pooled permutation distribution |
| `:nonparametric_individual` | Point-specific thresholds from permutation distributions |

## Cluster-Based Permutation Testing

Addresses the multiple comparisons problem in ERP analysis:

1. Compute test statistic at each time point × channel
2. Identify clusters of contiguous spatiotemporal significant effects
3. Permute condition labels and repeat (e.g., 1000 times)
4. Compare observed cluster mass to permutation distribution

**Key advantages:**
- Controls family-wise error rate without being overly conservative
- Sensitive to spatially and temporally distributed effects

## Workflow Summary

## Prepare Data

- Load group ERP data with `prepare_stats` specifying conditions, channels, and time intervals
- Supports `:paired` (within-subject) and `:independent` (between-subject) designs

## Analytic Tests

- Quick uncorrected or Bonferroni-corrected t-tests for exploration

## Permutation Tests

- Choose thresholding method based on assumptions and computational budget
- Set `cluster_type = :spatiotemporal` and `min_num_neighbors` for spatial constraints
- Visualise results with `plot_erp_stats`


## Code Examples

::: details Show Code

```julia
# Demo: Permutation-Based Statistical Testing for ERPs
# Shows how to run permutation tests on ERP data using different
# thresholding and cluster correction approaches.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

input_dir = EegFun.example_path("data/julia/erps")
file_pattern = "erps"

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
@time result_permutation_parametric = EegFun.permutation_test(
    stat_data,
    n_permutations = 50000,           # Number of permutations (more = more accurate)
    threshold_method = :parametric,  # Use t-distribution for thresholding
    cluster_type = :spatiotemporal,  # Cluster in space AND time
    min_num_neighbors = 3,           # Pre-filter: need 3+ neighbors
    show_progress = true,            # Show progress bar
    use_gpu = false,
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

@time result_permutation_nonparametric_common = EegFun.permutation_test(
    stat_data,
    n_permutations = 10000,
    threshold_method = :nonparametric_common,  # Non-parametric common threshold
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
    use_gpu = true,
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

@time result_permutation_nonparametric_individual = EegFun.permutation_test(
    stat_data,
    n_permutations = 10000,
    threshold_method = :nonparametric_individual,  # Non-parametric individual thresholds
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
    use_gpu = true,
)


```

:::

## See Also

- [API Reference](../../reference/index.md)
