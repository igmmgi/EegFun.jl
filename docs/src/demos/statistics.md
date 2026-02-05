# Statistics

## Overview

## Overview

This demo demonstrates statistical analysis methods for ERP data.

### Statistical Testing

Compare experimental conditions using appropriate statistics:

**Parametric Tests:**
- **T-tests**: Compare two conditions
- **ANOVA**: Multiple conditions or factors
- **Assumptions**: Normality, equal variance

**Non-parametric Tests:**
- **Cluster-based permutation**: Control for multiple comparisons
- **No distribution assumptions**: More robust
- **Spatial-temporal clustering**: Accounts for dependencies

### Cluster-Based Permutation Tests

Powerful method for ERP analysis:
1. Compute test statistic at each time point/channel
2. Find clusters of contiguous significance
3. Permute condition labels and repeat
4. Compare observed clusters to permutation distribution

### Advantages

- Controls family-wise error rate
- Sensitive to spatiotemporal effects
- No stringent parametric assumptions
- Accounts for multiple comparisons elegantly

### Applications

- Condition comparisons in ERP studies
- Group differences
- Time-frequency power comparisons


## Code Examples

::: details Show Code

```julia
"""
Tutorial: Statistical Analysis Options for ERP Data

This script provides an introduction to the statistical analysis
options available in EegFun, loosely based on FieldTrip's approach. 

1. Data preparation for statistical tests
2. Analytic t-tests (with/without multiple comparison correction)
3. Cluster-based permutation tests (with different thresholding methods)
4. Visualization of results
"""

using EegFun
using BenchmarkTools

input_dir = "./resources/data/erps"
file_pattern = "erps_good"

println("Preparing data...")
stat_data = EegFun.prepare_stats(
    file_pattern,  # Pattern to match ERP files
    :paired;       # :paired for within-subject, :independent for between-subject
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]), # Conditions to compare
    channel_selection = EegFun.channels(1:72),       # Select all 72 channels
    sample_selection = EegFun.samples((-0.5, 2.0)),  # Full time window
    baseline_window = EegFun.samples((-0.2, 0.0)),   # Baseline: -200 to 0 ms
    analysis_window = EegFun.samples((0.0, 1.0)),    # Analysis window: 100-1000 ms
)

# ----------------------------------------------------------------------------
# Option: Analytic t-test with NO correction
# ----------------------------------------------------------------------------
result_analytic_no = EegFun.analytic_test(stat_data)

# TODO: add plotting types e.g. :grid, :layout, :topo 
EegFun.plot_analytic_test(
    result_analytic_no,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
)

# ----------------------------------------------------------------------------
# Option: Analytic t-test with BONFERRONI correction (v. conservative; really only useful for demonstrating)
# ----------------------------------------------------------------------------
result_analytic_bonf = EegFun.analytic_test(stat_data, correction_method = :bonferroni)

EegFun.plot_analytic_test(
    result_analytic_bonf,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
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

EegFun.plot_analytic_test(
    result_permutation_parametric,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
)


# ----------------------------------------------------------------------------
# Option: Non-Parametric Common Thresholding
# ----------------------------------------------------------------------------
# What actually happens:
# 1. Run ALL permutations first (collect t-matrices)
# 2. Pool all t-values from all permutations
# 3. Compute (1-α) percentile threshold from pooled distribution
# 4. Use this single threshold for all points
# 5. Reuse stored t-matrices for cluster-level inference
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
# 1. Run ALL permutations first (collect t-matrices)
# 2. For each electrode × time point:
#    - Extract t-values from all permutations at that point
#    - Compute (1-α) percentile threshold from point-specific distribution
# 3. Threshold observed data using point-specific thresholds
# 4. Reuse stored t-matrices for cluster-level inference
# Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster', clusterThreshold='nonparametric_individual'

@btime result_permutation_nonparametric_individual = EegFun.permutation_test(
    stat_data,
    n_permutations = 1000,
    threshold_method = :nonparametric_individual,  # Non-parametric individual thresholds
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
)


```

:::

## See Also

- [API Reference](../reference/index.md)
