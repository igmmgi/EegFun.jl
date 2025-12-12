# Statistics Module Summary

## Overview

The `statistics.jl` module provides comprehensive statistical testing functions for EEG/ERP data analysis, implementing cluster-based permutation tests and analytic t-tests following the FieldTrip approach.

## Main Functions

### 1. Data Preparation

#### `prepare_condition_comparison`

Prepares ERP data for statistical testing by organizing data into participant × electrode × time arrays.

**Two variants:**
- **Direct data version**: Accepts `Vector{ErpData}`
- **File-based version**: Loads data from JLD2 files matching a pattern

**Key features:**
- Validates design type (`:paired` or `:independent`)
- Ensures exactly 2 conditions for comparison
- Applies baseline correction
- Creates grand averages for plotting
- Handles channel, time, and participant selection

**Parameters:**
- `design::Symbol`: `:paired` (within-subject) or `:independent` (between-subject)
- `condition_selection::Function`: Select exactly 2 conditions
- `channel_selection::Function`: Filter channels
- `sample_selection::Function`: Select time points
- `baseline_window::Function`: Baseline correction window
- `analysis_window::Function`: Analysis time window

---

### 2. Cluster-Based Permutation Tests

#### `cluster_permutation_test`

Performs cluster-based permutation tests on prepared ERP data. This is the main function for robust statistical inference that controls for multiple comparisons.

**Parameters:**

- `prepared::StatisticalTestData`: Prepared data from `prepare_condition_comparison`
- `n_permutations::Int`: Number of permutations (default: 1000)
- `threshold::Float64`: P-value threshold (default: 0.05)
- `threshold_method::Symbol`: 
  - `:parametric` (default): Uses t-distribution for thresholding
  - `:nonparametric_common`: Single threshold from pooled permutation distribution
  - `:nonparametric_individual`: Point-specific thresholds from permutation distributions
- `cluster_type::Symbol`:
  - `:spatiotemporal` (default): Clusters across both space and time
  - `:spatial`: Clusters only across electrodes (at each time point)
  - `:temporal`: Clusters only across time (at each electrode)
- `cluster_statistic::Symbol`:
  - `:sum` (default): Sum of t-values in cluster
  - `:max`: Maximum t-value in cluster
  - `:size`: Number of points in cluster
  - `:wcm`: Weighted cluster mass
- `min_num_neighbors::Int`: Minimum neighboring significant channels required (FieldTrip's `minNumChannels`, default: 0)
- `tail::Symbol`: `:both` (default), `:left`, or `:right`
- `random_seed::Union{Int, Nothing}`: Random seed for reproducibility
- `show_progress::Bool`: Show progress bar (default: true)

**Returns:**
- `ClusterPermutationResult`: Complete results structure with clusters, masks, and statistics

**Algorithm:**
1. Compute observed t-statistics
2. Threshold data (parametric or non-parametric)
3. Find clusters of contiguous significant points
4. Compute cluster-level statistics
5. Run permutations to build null distribution
6. Compute cluster p-values by comparing observed to null distribution

---

### 3. Analytic T-Tests

#### `analytic_ttest`

Performs analytic (parametric) t-tests with optional multiple comparison correction.

**Parameters:**
- `prepared::StatisticalTestData`: Prepared data
- `alpha::Float64`: Significance threshold (default: 0.05)
- `tail::Symbol`: `:both` (default), `:left`, or `:right`
- `correction_method::Symbol`:
  - `:no` (default): No correction
  - `:bonferroni`: Bonferroni correction

**Returns:**
- `AnalyticTTestResult`: Results with t-statistics, p-values, and significance masks

**Note:** Analytic tests are faster but less robust than permutation tests. Bonferroni correction is very conservative.

---

## Supporting Functions

### Core Statistics

- **`compute_t_matrix`**: Computes t-statistics and p-values for all electrode × time points
  - Optimized version for raw arrays (used in permutation loop)
  - Original version for `StatisticalTestData`
  - Supports paired and independent designs

- **`compute_critical_t_values`**: Computes critical t-values for parametric thresholding

- **`_compute_p_matrix`**: Internal function to compute p-values from t-statistics

### Thresholding

- **`threshold_t_matrix_parametric`**: Thresholds t-matrix using parametric critical values
- **`threshold_t_matrix_parametric!`**: In-place version
- **`threshold_t_matrix_nonparametric`**: Thresholds using non-parametric thresholds
- **`threshold_t_matrix_nonparametric!`**: In-place version

### Non-Parametric Thresholding

- **`compute_nonparametric_threshold_common`**: Computes single threshold from pooled permutation distribution
- **`compute_nonparametric_threshold_individual`**: Computes point-specific thresholds from permutation distributions
- **`collect_permutation_t_matrices`**: Collects all permutation t-matrices for non-parametric methods

### Connectivity and Clustering

- **`build_connectivity_matrix`**: Builds spatial connectivity matrix from layout neighbors
  - Supports `:spatial`, `:temporal`, and `:spatiotemporal` cluster types
  - Uses `Layout.neighbours` for spatial connectivity

- **`prefilter_mask_by_neighbors`**: Pre-filters masks to remove isolated points (FieldTrip's `minNumChannels`)
- **`prefilter_mask_by_neighbors!`**: In-place version

- **`find_clusters`**: Main function to find clusters in thresholded data
- **`find_clusters_connected_components`**: Finds connected components using BFS
- **`set_cluster_polarity`**: Sets cluster polarity (positive/negative)

### Cluster Statistics

- **`compute_cluster_statistics`**: Computes cluster-level statistics
  - Supports `:sum`, `:max`, `:size`, `:wcm`
- **`_compute_cluster_statistics`**: Internal helper with electrode mapping
- **`_compute_cluster_statistics_only`**: Computes statistics without updating cluster objects

### Permutation Logic

- **`run_permutations`**: Main permutation loop
  - Generates null distribution of maximum cluster statistics
  - Supports parametric and non-parametric thresholding
  - Can reuse pre-computed permutation t-matrices

- **`shuffle_labels!`**: Shuffles condition labels for permutations
- **`shuffle_labels`**: Wrapper that creates shuffled data
- **`generate_swap_mask`**: Generates random swap mask for paired designs

### P-Value Computation

- **`compute_cluster_pvalues`**: Computes cluster p-values from permutation distribution
  - Compares observed cluster statistics to null distribution
  - Handles both positive and negative clusters

### Utility Functions

- **`validate_permutation_inputs`**: Validates input parameters
- **`_create_significance_mask`**: Creates significance mask with optional correction

---

## Basic T-Test Functions

Located in `statistics_ttest.jl`:

- **`paired_ttest`**: Paired (dependent samples) t-test
- **`independent_ttest`**: Independent (between-subjects) t-test

Both return `TTestResult` with `df`, `t`, and `p` values.

---

## Key Design Decisions

### Performance Optimizations

1. **Pre-allocated buffers**: Functions accept optional buffers to avoid allocations in permutation loops
2. **Single-pass loops**: Combined mean and std calculations to reduce memory usage
3. **Vectorized operations**: Where possible, uses vectorized operations
4. **In-place functions**: Mutating versions (`!`) for critical loops

### FieldTrip Compatibility

- Follows FieldTrip's cluster permutation test approach
- Supports same thresholding methods
- Uses same cluster statistics
- Implements `minNumChannels` pre-filtering

### Code Organization

- **Phase-based structure**: Code organized into logical phases (preparation, thresholding, clustering, etc.)
- **Internal functions**: Helper functions prefixed with `_` for internal use
- **Clear separation**: Basic t-tests in separate file (`statistics_ttest.jl`)

---

## Result Types

### `ClusterPermutationResult`

Contains:
- `test_info::TestInfo`: Test configuration
- `data::Vector{ErpData}`: Grand averages for plotting
- `stat_matrix::StatMatrix`: T-statistics matrix
- `masks::Masks`: Significance masks
- `clusters::Clusters`: Positive and negative clusters
- `permutation_distribution::PermutationDistribution`: Null distribution
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points
- `critical_t`: Critical t-values (varies by threshold method)

### `AnalyticTTestResult`

Contains:
- `test_info::TestInfo`: Test configuration
- `data::Vector{ErpData}`: Grand averages
- `stat_matrix::StatMatrix`: T-statistics and p-values
- `masks::Masks`: Significance masks
- `electrodes::Vector{Symbol}`: Electrode labels
- `time_points::Vector{Float64}`: Time points
- `critical_t::Float64`: Critical t-value

---

## Usage Workflow

1. **Prepare data**:
   ```julia
   prepared = prepare_condition_comparison(
       "erps_pattern", :paired;
       input_dir = "data/",
       condition_selection = conditions([1, 2]),
       baseline_window = samples((-0.2, 0.0)),
       analysis_window = samples((0.1, 1.0))
   )
   ```

2. **Run cluster permutation test**:
   ```julia
   result = cluster_permutation_test(
       prepared,
       n_permutations = 1000,
       threshold = 0.05,
       threshold_method = :parametric,
       cluster_type = :spatiotemporal,
       cluster_statistic = :sum,
       min_num_neighbors = 3
   )
   ```

3. **Access results**:
   ```julia
   # Significant clusters
   sig_pos = [c for c in result.clusters.positive if c.is_significant]
   sig_neg = [c for c in result.clusters.negative if c.is_significant]
   
   # Significance mask
   sig_mask = result.masks.positive .| result.masks.negative
   ```

---

## Notes

- All functions are optimized for performance, especially the permutation loop
- The implementation follows FieldTrip's approach for maximum compatibility
- Cluster permutation tests are recommended for most analyses
- Analytic tests with Bonferroni are very conservative
- Non-parametric thresholding methods are more robust but slower

