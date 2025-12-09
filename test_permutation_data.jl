#!/usr/bin/env julia
"""
Tutorial: Statistical Analysis Options for ERP Data

This script provides a comprehensive introduction to the statistical analysis
options available in eegfun, following FieldTrip's approach. It demonstrates:

1. Data preparation for statistical tests
2. Analytic t-tests (with/without multiple comparison correction)
3. Cluster-based permutation tests (with different thresholding methods)
4. Visualization of results

Each section includes explanations of when and why to use each method.
"""

using eegfun
using BenchmarkTools

println("\n" * "="^80)
println("TUTORIAL: Statistical Analysis Options for ERP Data")
println("="^80)
println("\nThis tutorial demonstrates the statistical analysis options available")
println("in eegfun, which closely follow FieldTrip's implementation.\n")

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================
println("\n" * "="^80)
println("PART 1: Data Preparation")
println("="^80)
println("""
Before running any statistical test, we need to prepare the data:

1. Load ERP data from files
2. Select participants, channels, and time windows
3. Apply baseline correction (optional but recommended)
4. Organize data into conditions for comparison

The prepare_statistical_test_data() function handles all of this.
""")

# Data directory - adjust this to your data location
input_dir = "/home/ian/Documents/Julia/output_data/filtered_erps_good_lp_30hz"
file_pattern = "erps_good"


println("Preparing data...")
prepared = eegfun.prepare_statistical_test_data(
    file_pattern,  # Pattern to match ERP files
    :paired;       # :paired for within-subject, :independent for between-subject
    input_dir = input_dir,
    participant_selection = eegfun.participants(3:18),  # Select participants 3-18
    condition_selection = eegfun.conditions([1, 2]),    # Conditions to compare
    channel_selection = eegfun.channels(1:72),          # Select all 72 channels
    sample_selection = eegfun.samples((-0.5, 2.0)),     # Full time window
    baseline_window = eegfun.samples((-0.2, 0.0)),      # Baseline: -200 to 0 ms
    analysis_window = eegfun.samples((0.1, 1.0)),       # Analysis window: 100-1000 ms
)

println("✓ Data prepared successfully!")
println("  - Design: ", prepared.design)
println("  - Participants: ", length(prepared.ids1), " in condition 1, ", 
        length(prepared.ids2), " in condition 2")
println("  - Channels: ", length(prepared.electrodes))
println("  - Time points: ", length(prepared.time_points), " (", 
        first(prepared.time_points), " to ", last(prepared.time_points), " s)")

# ============================================================================
# PART 2: ANALYTIC T-TESTS (Parametric, No Permutation)
# ============================================================================
println("\n" * "="^80)
println("PART 2: Analytic T-Tests (Parametric Tests)")
println("="^80)
println("""
Analytic t-tests are fast parametric tests that use the t-distribution.
They are useful for:
- Quick exploratory analysis
- When you need fast results
- Demonstrating basic statistical concepts

However, they require multiple comparison correction when testing many
electrode × time points simultaneously.

Available correction methods:
- :no - No correction (uncorrected, many false positives)
- :bonferroni - Bonferroni correction (conservative, controls family-wise error)
""")

# ----------------------------------------------------------------------------
# Option 2a: Analytic t-test with NO correction
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 2a: Analytic T-Test with NO Correction")
println("-"^80)
println("""
This is the simplest test - just computes t-statistics and p-values at each
electrode × time point without correcting for multiple comparisons.

⚠️  WARNING: Without correction, you will have many false positives!
This is mainly useful for demonstration or when you only test a few points.

Equivalent to FieldTrip: method='analytic', corrMethod='no'
""")

result_analytic_no = eegfun.analytic_ttest(
    prepared,
    alpha = 0.05,           # Significance threshold
    tail = :both,           # Two-tailed test (:both, :left, or :right)
    correction_method = :no # No multiple comparison correction
)

println(result_analytic_no)
println("\nNote: This test found ", 
        count(result_analytic_no.significant_mask_positive) + 
        count(result_analytic_no.significant_mask_negative), 
        " significant points (likely many false positives!)")

fig = eegfun.plot_analytic_ttest(result_analytic_no, prepared, channel = :PO8, plot_erp = true, 
plot_difference = false, show_significance = true, show_critical_t = true)
display(fig)

# ----------------------------------------------------------------------------
# Option 2b: Analytic t-test with BONFERRONI correction
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 2b: Analytic T-Test with BONFERRONI Correction")
println("-"^80)
println("""
Bonferroni correction divides the alpha level by the number of tests.
This is very conservative but controls the family-wise error rate.

✓ Conservative: Very few false positives
✗ May miss real effects: High false negative rate

Equivalent to FieldTrip: method='analytic', corrMethod='bonferroni'
""")

result_analytic_bonf = eegfun.analytic_ttest(
    prepared,
    alpha = 0.05,
    tail = :both,
    correction_method = :bonferroni  # Bonferroni correction
)

fig = eegfun.plot_analytic_ttest(result_analytic_no, channel = :PO9, plot_erp = true, 
plot_difference = false, show_significance = true, show_critical_t = true)
display(fig)


println(result_analytic_bonf)
println("\nNote: Bonferroni correction is very conservative - compare the")
println("      number of significant points to the uncorrected test above.")

# ============================================================================
# PART 3: CLUSTER-BASED PERMUTATION TESTS (Monte Carlo)
# ============================================================================
println("\n" * "="^80)
println("PART 3: Cluster-Based Permutation Tests (Monte Carlo)")
println("="^80)
println("""
Cluster-based permutation tests are the recommended method for ERP/EEG data.
They:
- Control for multiple comparisons at the cluster level
- Are non-parametric (make fewer assumptions)
- Are more powerful than Bonferroni correction
- Account for spatial and temporal correlations

The key idea: Instead of testing each point individually, we:
1. Threshold the data (find "significant" points)
2. Form clusters of adjacent significant points
3. Test whether cluster-level statistics are significant
4. Use permutation to create a null distribution

This method is more appropriate for neuroimaging data because it accounts
for the fact that neighboring electrodes and time points are correlated.
""")

println("\n" * "="^80)
println("Thresholding Methods Explained")
println("="^80)
println("""
Before clustering, we need to decide which points are "significant enough"
to include. There are three thresholding methods:

1. PARAMETRIC (fastest, default)
   - Uses t-distribution to compute critical t-values
   - Assumes data follows t-distribution
   - Fast: No need to run permutations first
   
2. NON-PARAMETRIC COMMON (more robust)
   - Runs all permutations first
   - Computes single threshold from permutation distribution
   - Avoids distributional assumption
   - Slower but more robust
   
3. NON-PARAMETRIC INDIVIDUAL (most robust, most memory-intensive)
   - Runs all permutations first
   - Computes point-specific thresholds
   - Most sensitive to point-specific variability
   - Slowest and most memory-intensive

All three methods use the SAME permutation-based inference for clusters.
The difference is only in the initial thresholding step.
""")

# ----------------------------------------------------------------------------
# Option 3a: Parametric Thresholding (Default, Fastest)
# ----------------------------------------------------------------------------
println("\n" * "-"^80)
println("Option 3a: Cluster Permutation Test - PARAMETRIC Thresholding")
println("-"^80)
println("""
This is the default and fastest method. It uses the t-distribution to
determine which points exceed the significance threshold.

✓ Fast: No need to run permutations for thresholding
✓ Standard approach: Most commonly used
✓ Good balance: Speed vs. robustness

Parameters explained:
- n_permutations: Number of permutations (1000 is standard, more is better)
- threshold: Significance level (0.05 = 5%)
- threshold_method: :parametric (uses t-distribution)
- min_num_neighbors: Pre-filter isolated points (FieldTrip's minNumChannels)
  - Removes points with fewer than N neighboring significant points
  - Helps reduce noise before clustering

Equivalent to FieldTrip: method='montecarlo', corrMethod='cluster', 
                          clusterThreshold='parametric'
""")

result_cluster_parametric = eegfun.cluster_permutation_test(
    prepared,
    n_permutations = 1000,        # Number of permutations (more = more accurate)
    threshold = 0.05,             # Significance level
    threshold_method = :parametric,  # Use t-distribution for thresholding
    cluster_type = :spatiotemporal,  # Cluster in space AND time
    cluster_statistic = :sum,     # Sum of t-values in cluster (default)
    min_num_neighbors = 3,       # Pre-filter: need 3+ neighbors
    tail = :both,                 # Two-tailed test
    show_progress = true          # Show progress bar
)

fig = eegfun.plot_analytic_ttest(result_cluster_parametric, prepared, channel = :PO8, plot_erp = true, plot_difference = false, show_significance = true, show_critical_t = true)
display(fig)

println(result_cluster_parametric)

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

@btime result_cluster_nonparametric_common = eegfun.cluster_permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :nonparametric_common,  # Non-parametric common threshold
    cluster_type = :spatiotemporal,
    cluster_statistic = :sum,
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true
)

println(result_cluster_nonparametric_common)

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

result_cluster_nonparametric_individual = eegfun.cluster_permutation_test(
    prepared,
    n_permutations = 1000,
    threshold = 0.05,
    threshold_method = :nonparametric_common,  # Non-parametric individual thresholds
    cluster_type = :spatiotemporal,
    cluster_statistic = :sum,
    min_num_neighbors = 3,
    tail = :both,
    show_progress = true
)

println(result_cluster_nonparametric_individual)

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
println("  Positive clusters: ", count(c -> c.is_significant, result_cluster_parametric.positive_clusters), 
        " significant out of ", length(result_cluster_parametric.positive_clusters), " total")
println("  Negative clusters: ", count(c -> c.is_significant, result_cluster_parametric.negative_clusters),
        " significant out of ", length(result_cluster_parametric.negative_clusters), " total")

println("\nNon-parametric common thresholding:")
println("  Positive clusters: ", count(c -> c.is_significant, result_cluster_nonparametric_common.positive_clusters),
        " significant out of ", length(result_cluster_nonparametric_common.positive_clusters), " total")
println("  Negative clusters: ", count(c -> c.is_significant, result_cluster_nonparametric_common.negative_clusters),
        " significant out of ", length(result_cluster_nonparametric_common.negative_clusters), " total")

println("\nNon-parametric individual thresholding:")
println("  Positive clusters: ", count(c -> c.is_significant, result_cluster_nonparametric_individual.positive_clusters),
        " significant out of ", length(result_cluster_nonparametric_individual.positive_clusters), " total")
println("  Negative clusters: ", count(c -> c.is_significant, result_cluster_nonparametric_individual.negative_clusters),
        " significant out of ", length(result_cluster_nonparametric_individual.negative_clusters), " total")

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
fig1 = eegfun.plot_analytic_ttest(
    result_analytic_no, 
    prepared, 
    channel = :PO7,
    plot_erp = true,
    plot_difference = false,
    show_significance = false
)
display(fig1)

# Example 2: Difference wave with significance
println("\n" * "-"^80)
println("Example 2: Difference Wave with Significance Markers")
println("-"^80)
fig2 = eegfun.plot_analytic_ttest(
    result_analytic_no,
    prepared,
    channel = :PO7,
    plot_erp = false,
    plot_difference = true,
    show_significance = true,
    show_critical_t = true
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
fig3 = eegfun.plot_analytic_ttest(
    result_cluster_parametric,
    prepared,
    channel = :PO7,
    plot_erp = false,
    plot_difference = true,
    show_significance = true,
    show_critical_t = true,
    sig_bar_position = 4.0,
    sig_bar_color = (:gray, 0.6)
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
- cluster_statistic: :sum (default), :max, :size, or :wcm
- cluster_type: :spatiotemporal (default), :spatial, or :temporal

For more information, see:
- FieldTrip documentation: https://www.fieldtriptoolbox.org/tutorial/stats/statistics/
- eegfun documentation: IMPLEMENTED_FEATURES.md
""")

println("\n" * "="^80)
println("Tutorial Complete!")
println("="^80)
println("\nYou can now use these functions in your own analyses.")
println("Adjust the parameters based on your specific research questions.\n")
