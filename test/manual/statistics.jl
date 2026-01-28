"""
Tutorial: Statistical Analysis Options for ERP Data

This script provides a comprehensive introduction to the statistical analysis
options available in EegFun, following FieldTrip's approach. It demonstrates:

1. Data preparation for statistical tests
2. Analytic t-tests (with/without multiple comparison correction)
3. Cluster-based permutation tests (with different thresholding methods)
4. Visualization of results

Each section includes explanations of when and why to use each method.
"""

using EegFun
using BenchmarkTools


input_dir = "/home/ian/Documents/Julia/output_data/filtered_erps_good_lp_30hz"
file_pattern = "erps_good"

println("Preparing data...")
prepared = EegFun.prepare_stats(
    file_pattern,  # Pattern to match ERP files
    :paired;       # :paired for within-subject, :independent for between-subject
    input_dir = input_dir,
    participant_selection = EegFun.participants(3:18),  # Select participants 3-18
    condition_selection = EegFun.conditions([1, 2]),    # Conditions to compare
    channel_selection = EegFun.channels(1:72),          # Select all 72 channels
    sample_selection = EegFun.samples((-0.5, 2.0)),     # Full time window
    baseline_window = EegFun.samples((-0.2, 0.0)),      # Baseline: -200 to 0 ms
    analysis_window = EegFun.samples((0.1, 1.0)),       # Analysis window: 100-1000 ms
)

# ----------------------------------------------------------------------------
# Option 2a: Analytic t-test with NO correction
# ----------------------------------------------------------------------------
result_analytic_no = EegFun.analytic_test(
    prepared,
    alpha = 0.05,           # Significance threshold
    tail = :both,           # Two-tailed test (:both, :left, or :right)
    correction_method = :no, # No multiple comparison correction
)

fig = EegFun.plot_analytic_test(
    result_analytic_no,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
)

# ----------------------------------------------------------------------------
# Option 2b: Analytic t-test with BONFERRONI correction
# ----------------------------------------------------------------------------
result_analytic_bonf = EegFun.analytic_test(
    prepared,
    alpha = 0.05,
    tail = :both,
    correction_method = :bonferroni,  # Bonferroni correction
)

fig = EegFun.plot_analytic_test(
    result_analytic_bonf,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
)


println(result_analytic_bonf)
println("\nNote: Bonferroni correction is very conservative - compare the")
println("      number of significant points to the uncorrected test above.")

# ============================================================================
# PART 3: CLUSTER-BASED PERMUTATION TESTS (Monte Carlo)
# ============================================================================

# ----------------------------------------------------------------------------
# Option 3a: Parametric Thresholding (Default, Fastest)
# ----------------------------------------------------------------------------

result_permutation_parametric = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,        # Number of permutations (more = more accurate)
    threshold = 0.05,             # Significance level
    threshold_method = :parametric,  # Use t-distribution for thresholding
    cluster_type = :spatiotemporal,  # Cluster in space AND time
    min_num_neighbors = 3,       # Pre-filter: need 3+ neighbors
    tail = :both,                 # Two-tailed test
    show_progress = true,          # Show progress bar
)

fig = EegFun.plot_analytic_test(
    result_permutation_parametric,
    channel = :PO8,
    plot_erp = true,
    plot_difference = false,
    show_significance = true,
    show_critical_t = true,
)
display(fig)

println(result_permutation_parametric)

# ----------------------------------------------------------------------------
# Option 3b: Non-Parametric Common Thresholding
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3b: Cluster Permutation Test - NON-PARAMETRIC COMMON")
println("-"^80)
println("""
This method avoids the parametric assumption by using the permutation
distribution itself to determine the threshold.

How it works:
1. Run ALL permutations first (collect t-matrices)
2. Pool all t-values from all permutations
3. Compute (1-α) percentile threshold from pooled distribution
4. Use this single threshold for all points
5. Reuse stored t-matrices for cluster-level inference

✓ More robust: No distributional assumption
✓ Single threshold: Easier to interpret
✗ Slower: Must run permutations twice (once for threshold, once for inference)
✗ Less sensitive: Doesn't account for point-specific variability

Use when: You want robustness but don't need point-specific thresholds.

Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster',
                          clusterThreshold='nonparametric_common'
""")

result_permutation_nonparametric_common = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :nonparametric_common,  # Non-parametric common threshold
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true,
)

println(result_permutation_nonparametric_common)

# ----------------------------------------------------------------------------
# Option 3c: Non-Parametric Individual Thresholding
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3c: Cluster Permutation Test - NON-PARAMETRIC INDIVIDUAL")
println("-"^80)
println("""
This is the most sophisticated method. It computes a different threshold
for each electrode × time point based on that point's permutation distribution.

How it works:
1. Run ALL permutations first (collect t-matrices)
2. For each electrode × time point:
   - Extract t-values from all permutations at that point
   - Compute (1-α) percentile threshold from point-specific distribution
3. Threshold observed data using point-specific thresholds
4. Reuse stored t-matrices for cluster-level inference

✓ Most robust: No distributional assumption
✓ Most sensitive: Accounts for point-specific variability
✗ Slowest: Must run permutations twice
✗ Memory-intensive: Stores all permutation t-matrices
✗ Complex: Point-specific thresholds can be harder to interpret

Use when: You need maximum robustness and have sufficient memory/time.

Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster',
                          clusterThreshold='nonparametric_individual'
""")

result_permutation_nonparametric_individual = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :nonparametric_individual,  # Non-parametric individual thresholds
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true,
)

println(result_permutation_nonparametric_individual)

# ============================================================================
# PART 3.5: CLUSTER TYPE OPTIONS
# ============================================================================
println("\n" * "="^80)
println("PART 3.5: Cluster Type Options")
println("="^80)
println("""
The cluster_type parameter determines how points are grouped into clusters.
Different cluster types are useful for different research questions.
""")

# ----------------------------------------------------------------------------
# Option 3d: Spatial Clustering Only
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3d: SPATIAL Clustering")
println("-"^80)
println("""
Spatial clustering groups adjacent electrodes at the same time point.
This is useful when you want to find spatially contiguous effects
regardless of their temporal duration.

Use when:
- You're interested in spatial patterns (e.g., which brain regions)
- Temporal duration is less important
- You want to find effects that span multiple electrodes simultaneously

Example: Finding which electrode groups show effects at any time point.
""")

result_permutation_spatial = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :parametric,
    cluster_type = :spatial,  # Only spatial adjacency
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true,
)

println(result_permutation_spatial)
println("\nSpatial clusters found:")
println("  Positive clusters: ", length(result_permutation_spatial.clusters.positive))
println("  Negative clusters: ", length(result_permutation_spatial.clusters.negative))

# ----------------------------------------------------------------------------
# Option 3e: Temporal Clustering Only
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3e: TEMPORAL Clustering")
println("-"^80)
println("""
Temporal clustering groups consecutive time points at the same electrode.
This is useful when you want to find temporally contiguous effects
regardless of their spatial extent.

Use when:
- You're interested in temporal dynamics (e.g., when effects occur)
- Spatial extent is less important
- You want to find effects that persist over time at individual electrodes

Example: Finding sustained effects at specific electrodes over time.
""")

result_permutation_temporal = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :parametric,
    cluster_type = :temporal,  # Only temporal adjacency
    min_num_neighbors = 0,  # Not applicable for temporal-only
    tail = :both,
    show_progress = true,
)

println(result_permutation_temporal)
println("\nTemporal clusters found:")
println("  Positive clusters: ", length(result_permutation_temporal.clusters.positive))
println("  Negative clusters: ", length(result_permutation_temporal.clusters.negative))

# ----------------------------------------------------------------------------
# Option 3f: Spatiotemporal Clustering (Default)
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3f: SPATIOTEMPORAL Clustering (Default)")
println("-"^80)
println("""
Spatiotemporal clustering groups points that are adjacent in space OR time.
This is the most common choice for ERP/EEG data because it captures
effects that are contiguous in both dimensions.

Use when:
- You want to find effects that are contiguous in space and/or time
- This is the standard approach for most ERP/EEG analyses
- You want maximum sensitivity to detect neural responses

Example: Finding ERP components that show spatial and temporal clustering.
""")

result_permutation_spatiotemporal = EegFun.permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :parametric,
    cluster_type = :spatiotemporal,  # Both spatial and temporal
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true,
)

println(result_permutation_spatiotemporal)
println("\nSpatiotemporal clusters found:")
println("  Positive clusters: ", length(result_permutation_spatiotemporal.clusters.positive))
println("  Negative clusters: ", length(result_permutation_spatiotemporal.clusters.negative))

# ----------------------------------------------------------------------------
# Comparison of Cluster Types
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Comparison of Cluster Types")
println("-"^80)
println("""
Different cluster types will find different patterns:

Spatial only:
- Finds spatially contiguous effects
- May split temporally separated effects into multiple clusters
- Useful for identifying brain regions

Temporal only:
- Finds temporally contiguous effects
- May split spatially separated effects into multiple clusters
- Useful for identifying time windows

Spatiotemporal (default):
- Finds effects contiguous in space OR time
- Most sensitive (finds more clusters)
- Standard approach for ERP/EEG data
- Recommended for most analyses

The same data can produce different numbers of clusters depending on
the cluster_type, because the clustering algorithm groups points differently.
""")

println("\nCluster counts comparison:")
println(
    "  Spatial:          ",
    length(result_permutation_spatial.clusters.positive),
    " positive, ",
    length(result_permutation_spatial.clusters.negative),
    " negative",
)
println(
    "  Temporal:         ",
    length(result_permutation_temporal.clusters.positive),
    " positive, ",
    length(result_permutation_temporal.clusters.negative),
    " negative",
)
println(
    "  Spatiotemporal:  ",
    length(result_permutation_spatiotemporal.clusters.positive),
    " positive, ",
    length(result_permutation_spatiotemporal.clusters.negative),
    " negative",
)

# ============================================================================
# PART 4: COMPARISON OF METHODS
# ============================================================================
println("\n" * "="^80)
println("PART 4: Comparison of Methods")
println("="^80)
println("""
Let's compare the results from different methods to see how they differ.
""")

println("\n" * "-"^80)
println("Significant Clusters Found:")
println("-"^80)
println("\nParametric thresholding:")
println(
    "  Positive clusters: ",
    count(c -> c.is_significant, result_permutation_parametric.clusters.positive),
    " significant out of ",
    length(result_permutation_parametric.clusters.positive),
    " total",
)
println(
    "  Negative clusters: ",
    count(c -> c.is_significant, result_permutation_parametric.clusters.negative),
    " significant out of ",
    length(result_permutation_parametric.clusters.negative),
    " total",
)

println("\nNon-parametric common thresholding:")
println(
    "  Positive clusters: ",
    count(c -> c.is_significant, result_permutation_nonparametric_common.clusters.positive),
    " significant out of ",
    length(result_permutation_nonparametric_common.clusters.positive),
    " total",
)
println(
    "  Negative clusters: ",
    count(c -> c.is_significant, result_permutation_nonparametric_common.clusters.negative),
    " significant out of ",
    length(result_permutation_nonparametric_common.clusters.negative),
    " total",
)

println("\nNon-parametric individual thresholding:")
println(
    "  Positive clusters: ",
    count(c -> c.is_significant, result_permutation_nonparametric_individual.clusters.positive),
    " significant out of ",
    length(result_permutation_nonparametric_individual.clusters.positive),
    " total",
)
println(
    "  Negative clusters: ",
    count(c -> c.is_significant, result_permutation_nonparametric_individual.clusters.negative),
    " significant out of ",
    length(result_permutation_nonparametric_individual.clusters.negative),
    " total",
)

println("\n" * "-"^80)
println("Interpretation:")
println("-"^80)
println("""
All three methods use the SAME cluster-level inference (permutation-based).
The differences you see are due to:
1. Different initial thresholds (which points are included in clusters)
2. Different cluster formation (which points cluster together)

In general:
- Parametric: Standard, fast, widely used
- Non-parametric common: More robust, single threshold
- Non-parametric individual: Most robust, point-specific thresholds

For most analyses, parametric thresholding is sufficient and recommended.
Use non-parametric methods when:
- You have concerns about distributional assumptions
- You want to be extra conservative
- You have sufficient computational resources
""")

# ============================================================================
# PART 5: VISUALIZATION
# ============================================================================
println("\n" * "="^80)
println("PART 5: Visualization")
println("="^80)
println("""
The plot_analytic_ttest() function can visualize results from both analytic
and cluster permutation tests. It's highly configurable with boolean flags.

Available plot components:
- plot_erp: Show ERP waveforms for both conditions
- plot_difference: Show difference wave (condition A - condition B)
- plot_tvalues: Show t-statistics over time
- show_significance: Show significance markers (grey bars)
- show_critical_t: Show critical t-value lines

Let's create some example plots...
""")

using GLMakie  # or CairoMakie for static plots

# Example 1: Basic ERP waveforms
println("\n" * "-"^80)
println("Example 1: ERP Waveforms Only")
println("-"^80)
fig1 = EegFun.plot_analytic_test(result_analytic_no, channel = :PO7, plot_erp = true, plot_difference = false, show_significance = false)
display(fig1)

# Example 2: Difference wave with significance
println("\n" * "-"^80)
println("Example 2: Difference Wave with Significance Markers")
println("-"^80)
fig2 = EegFun.plot_analytic_test(
    result_analytic_no,
    channel = :PO7,
    plot_erp = false,
    plot_difference = true,
    show_significance = true,
    show_critical_t = true,
)
display(fig2)

# Example 3: Cluster permutation results
println("\n" * "-"^80)
println("Example 3: Cluster Permutation Test Results (Parametric)")
println("-"^80)
println("""
This shows the results from cluster-based permutation testing.
Significance markers indicate points that are part of significant clusters.
""")
fig3 = EegFun.plot_analytic_test(
    result_permutation_parametric,
    channel = :PO7,
    plot_erp = false,
    plot_difference = true,
    show_significance = true,
    show_critical_t = true,
    sig_bar_position = 4.0,
    sig_bar_color = (:gray, 0.6),
)
display(fig3)

# ============================================================================
# PART 6: SUMMARY AND RECOMMENDATIONS
# ============================================================================
println("\n" * "="^80)
println("PART 6: Summary and Recommendations")
println("="^80)
println("""
WHEN TO USE EACH METHOD:

1. ANALYTIC T-TEST (no correction)
   → Quick exploratory analysis
   → Demonstrations/teaching
   → ⚠️  NOT recommended for publication (too many false positives)

2. ANALYTIC T-TEST (Bonferroni correction)
   → When you need fast results
   → When you have few tests
   → ⚠️  Very conservative, may miss real effects

3. CLUSTER PERMUTATION TEST (parametric thresholding) ⭐ RECOMMENDED
   → Standard approach for ERP/EEG data
   → Good balance of speed and robustness
   → Most commonly used in published research
   → ✅ Recommended for most analyses

4. CLUSTER PERMUTATION TEST (non-parametric common)
   → When you want extra robustness
   → When parametric assumptions are questionable
   → When you have moderate computational resources

5. CLUSTER PERMUTATION TEST (non-parametric individual)
   → Maximum robustness needed
   → When you have sufficient memory and time
   → For critical analyses where assumptions matter

KEY PARAMETERS TO CONSIDER:

- n_permutations: More = more accurate (1000 is standard, 10000 for publication)
- min_num_neighbors: Pre-filters noise (3 is common, 0 = no filtering)
- cluster_statistic: :maxsum (default) or :size
- cluster_type: :spatiotemporal (default), :spatial, or :temporal

For more information, see:
- FieldTrip documentation: https://www.fieldtriptoolbox.org/tutorial/stats/statistics/
- EegFun documentation: IMPLEMENTED_FEATURES.md
""")

println("\n" * "="^80)
println("Tutorial Complete!")
println("="^80)
println("\nYou can now use these functions in your own analyses.")
println("Adjust the parameters based on your specific research questions.\n")
