# Cluster-Based Permutation Tests: Implementation Guide

## Overview

Implementation of cluster-based permutation tests for ERP/EEG data following the FieldTrip approach. Supports both paired t-tests (within-subject) and independent sample t-tests (between-subject) with spatial and/or temporal clustering.

## Supported Designs

1. **Paired t-test**: Compare two conditions within participants (e.g., Condition A vs Condition B)
   - Data: `EpochData` or `ErpData` with condition labels
   - Design: Within-subject (repeated measures)

2. **Independent sample t-test**: Compare two groups (e.g., Group 1 vs Group 2)
   - Data: `EpochData` or `ErpData` with group labels
   - Design: Between-subject (independent groups)

## Algorithm Overview

### Step 1: Data Preparation
- Load/accept `EpochData` or `ErpData`
- Extract condition/group labels
- Organize data into design matrix structure
- Handle missing data/unbalanced designs

### Step 2: Compute Test Statistics
- For each electrode × time point:
  - Compute t-statistic (paired or independent)
  - Store t-values in electrode × time matrix

### Step 3: Threshold
- Apply uncorrected threshold (e.g., p < 0.05, or t > threshold)
- Create binary mask of significant points

### Step 4: Find Clusters
- **Spatial clustering**: Use electrode neighbors from `Layout`
- **Temporal clustering**: Consecutive time points
- **Spatiotemporal clustering**: Both spatial and temporal adjacency
- Identify connected components (clusters) in thresholded data

### Step 5: Compute Cluster Statistics
- For each cluster:
  - Sum of t-values within cluster (cluster-level statistic)
  - Alternative: max t-value, or cluster mass
- Store cluster statistics

### Step 6: Permutation
- Shuffle condition/group labels (respecting design)
- Recompute t-statistics
- Re-threshold and find clusters
- Store maximum cluster statistic from each permutation
- Repeat N times (typically 1000-10000)

### Step 7: Compute Cluster p-values
- Compare observed cluster statistics to permutation distribution
- For each cluster: p = (number of permutations with max cluster stat >= observed) / N
- Apply cluster-level correction

### Step 8: Return Results
- Cluster locations (electrodes, time windows)
- Cluster statistics (sum of t-values)
- Cluster p-values
- Permutation distribution
- Visualization-ready data

## Detailed Steps

### Step 1: Data Preparation

**Input:**
- `data::Union{Vector{EpochData}, Vector{ErpData}}`
- Condition/group labels
- Design type (`:paired` or `:independent`)

**Process:**
1. Extract data matrices: `[participants × electrodes × time]` or `[participants × electrodes × time × epochs]`
2. Organize by condition/group
3. Validate design:
   - Paired: same participants in both conditions
   - Independent: different participants in each group
4. Handle missing data (warn or exclude)

**Output:**
- Organized data arrays
- Design matrix
- Participant/condition mappings

### Step 2: Compute Test Statistics

**For each electrode × time point:**

**Paired t-test:**
```
t = mean(diff) / (std(diff) / sqrt(n))
where diff = condition_A - condition_B for each participant
```

**Independent t-test:**
```
t = (mean_A - mean_B) / sqrt(pooled_variance * (1/n_A + 1/n_B))
where pooled_variance = ((n_A-1)*var_A + (n_B-1)*var_B) / (n_A + n_B - 2)
```

**Implementation approach:**
- **Write from scratch** (recommended)
  - Simple formulas, easy to implement
  - Only need t-statistic value (not full test object)
  - Performance optimized for permutation loop
  - Use existing `Statistics.jl` and `StatsBase.jl` helpers:
    - `mean()`, `std()`, `var()` from Statistics
    - `pooledvar()` from StatsBase for independent t-test
  - No new dependencies needed

**Process:**
1. For each (electrode, time) combination:
   - Extract data for both conditions/groups
   - Compute t-statistic using formulas above
   - Store in `t_matrix[electrode, time]`

**Output:**
- `t_matrix`: Array of t-values [electrodes × time]

### Step 3: Threshold

**Options:**
- **P-value threshold**: Convert t to p, threshold at p < 0.05 (two-tailed)
- **T-value threshold**: Direct threshold on |t| > critical_t
- **Default**: p < 0.05 (two-tailed)

**Process:**
1. Convert t-values to p-values (or use t-threshold)
2. Create binary mask: `significant_mask[electrode, time] = (p < threshold)`
3. Apply two-tailed test: check both positive and negative clusters

**Output:**
- `significant_mask`: Boolean array [electrodes × time]

### Step 4: Find Clusters

**Spatial adjacency:**
- Use `Layout.neighbours` (already computed in eegfun)
- Two electrodes are neighbors if connected in neighbor graph

**Temporal adjacency:**
- Consecutive time points (sample i and i+1)

**Spatiotemporal clustering:**
- Point (electrode_i, time_j) is adjacent to:
  - Spatial: All neighbors of electrode_i at time_j
  - Temporal: electrode_i at time_j±1
  - Both must be significant

**Algorithm:**
1. Start with first significant point
2. Use graph traversal (BFS/DFS) to find all connected significant points
3. Label cluster
4. Move to next unlabeled significant point
5. Repeat until all points labeled

**Output:**
- `clusters`: Vector of cluster objects containing:
  - Cluster ID
  - Electrodes in cluster
  - Time points in cluster
  - Cluster statistic (sum of t-values)

### Step 5: Compute Cluster Statistics

**For each cluster:**
- **Cluster statistic** = sum of t-values within cluster
- Alternative options:
  - Max t-value in cluster
  - Cluster mass (sum of |t|)
  - Cluster size (number of points)

**Default:** Sum of t-values (FieldTrip default)

**Output:**
- `cluster_stats`: Vector of cluster statistics

### Step 6: Permutation

**Paired design:**
- Shuffle condition labels within each participant
- For each participant: randomly assign A/B labels
- Maintains within-subject structure

**Independent design:**
- Shuffle group labels across participants
- Randomly assign participants to groups
- Maintains group sizes

**Process:**
1. For each permutation (1 to N):
   - Shuffle labels according to design
   - Recompute t-statistics (Step 2)
   - Re-threshold (Step 3)
   - Find clusters (Step 4)
   - Compute cluster statistics (Step 5)
   - Store maximum cluster statistic (across all clusters, both positive and negative)

**Output:**
- `permutation_distribution`: Vector of max cluster statistics from N permutations

### Step 7: Compute Cluster p-values

**For each observed cluster:**
```
p_value = (number of permutations with max_cluster_stat >= observed_cluster_stat) / N
```

**Handle positive and negative clusters separately:**
- Positive clusters: compare to max positive cluster stats
- Negative clusters: compare to max negative cluster stats

**Significance threshold:**
- Typically p < 0.05 (cluster-level corrected)

**Output:**
- `cluster_pvalues`: Vector of p-values for each cluster
- `significant_clusters`: Clusters with p < 0.05

### Step 8: Return Results

**Output Structure:**
```julia
struct ClusterPermutationResult
    # Test information
    test_type::Symbol  # :paired or :independent
    design::DesignMatrix
    
    # Observed statistics
    t_matrix::Array{Float64, 2}  # [electrodes × time]
    significant_mask::BitArray{2}  # [electrodes × time]
    
    # Clusters
    clusters::Vector{Cluster}
    cluster_stats::Vector{Float64}
    cluster_pvalues::Vector{Float64}
    significant_clusters::Vector{Int}  # Indices of significant clusters
    
    # Permutation info
    n_permutations::Int
    permutation_distribution::Vector{Float64}
    
    # Metadata
    electrodes::Vector{Symbol}
    time_points::Vector{Float64}
    threshold::Float64
end

struct Cluster
    id::Int
    electrodes::Vector{Symbol}
    time_range::Tuple{Float64, Float64}  # (start, end) in seconds
    time_indices::UnitRange{Int}
    cluster_stat::Float64
    p_value::Float64
    is_significant::Bool
end
```

## Visualization Requirements

### 1. Cluster Topography
- Plot significant clusters on electrode layout
- Color-code by cluster
- Show time window for each cluster
- Overlay on head plot

### 2. Cluster Time Course
- Plot ERP waveforms for each condition/group
- Highlight significant time windows
- Mark cluster boundaries
- Show t-values over time

### 3. Spatiotemporal Cluster Map
- 2D plot: electrodes × time
- Color-code by t-value
- Outline significant clusters
- Show cluster IDs and p-values

### 4. Permutation Distribution
- Histogram of permutation distribution
- Mark observed cluster statistics
- Show significance threshold

### 5. Summary Table
- List all clusters with:
  - Cluster ID
  - Electrodes
  - Time window
  - Cluster statistic
  - p-value
  - Significant (yes/no)

## Implementation Checklist

### Phase 1: Core Statistics
- [ ] Data preparation function
- [ ] T-test computation (paired and independent)
- [ ] Thresholding function
- [ ] Basic result structure

### Phase 2: Clustering
- [ ] Spatial neighbor extraction from Layout
- [ ] Temporal adjacency definition
- [ ] Cluster finding algorithm (graph traversal)
- [ ] Cluster statistic computation

### Phase 3: Permutation
- [ ] Label shuffling (paired design)
- [ ] Label shuffling (independent design)
- [ ] Permutation loop
- [ ] Distribution storage

### Phase 4: Inference
- [ ] P-value computation
- [ ] Significance determination
- [ ] Result structure population

### Phase 5: Visualization
- [ ] Cluster topography plot
- [ ] Time course with clusters
- [ ] Spatiotemporal map
- [ ] Permutation distribution plot
- [ ] Summary table

## Design Decisions Needed

1. **Default threshold**: p < 0.05 or t-value threshold?
2. **Cluster statistic**: Sum of t-values (default) or other?
3. **Number of permutations**: Default 1000? Configurable?
4. **Handle both positive and negative**: Separate or combined?
5. **Multiple comparisons**: Only cluster correction, or also offer FDR/FWE?
6. **Output format**: Return full result struct or also save to file?

## FieldTrip Reference

**Key FieldTrip functions to reference:**
- `ft_timelockstatistics` - Main function
- `cfg.method = 'montecarlo'` - Permutation method
- `cfg.statistic = 'depsamplesT'` (paired) or `'indepsamplesT'` (independent)
- `cfg.clusterstatistic = 'maxsum'` - Cluster statistic method
- `cfg.clusteralpha = 0.05` - Cluster-forming threshold
- `cfg.clustercritval` - Cluster threshold

**Important FieldTrip behaviors:**
- Uses two-tailed test by default
- Handles positive and negative clusters separately
- Returns cluster p-values
- Provides cluster locations and statistics

## Questions for FieldTrip Code Review

1. Exact thresholding method (p-value vs t-value)?
2. How are positive/negative clusters handled?
3. Cluster statistic computation details?
4. Permutation strategy for paired designs?
5. How is the permutation distribution used?
6. Output structure and fields?

## Next Steps

1. **Review FieldTrip code** (if available) to match exact approach
2. **Implement Phase 1**: Basic t-test computation
3. **Test with simple data**: Known differences, verify statistics
4. **Implement Phase 2**: Clustering algorithm
5. **Implement Phase 3**: Permutation framework
6. **Implement Phase 4**: Inference and results
7. **Implement Phase 5**: Visualization functions
8. **Validate**: Compare results with FieldTrip on same data

