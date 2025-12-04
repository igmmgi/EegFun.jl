# Cluster-Based Permutation Test Implementation Plan

## Current Status
✅ `prepare_permutation_data` - Working
✅ `paired_ttest` and `independent_ttest` - Implemented
✅ `Layout.neighbours` - Available for spatial clustering
✅ Data structure: `PermutationTestData` - Ready

## Implementation Strategy

### Phase 0: Infrastructure and Validation (Start Here)
**Goal**: Set up foundation and validation

1. **Input Validation Function**
   - Check `Layout.neighbours` exists (compute if missing with default criterion)
   - Validate `n_permutations > 0`
   - Validate `threshold` in valid range (0 < threshold < 1)
   - Check data dimensions match
   - Validate design type (`:paired` or `:independent`)
   - Check minimum participants (at least 2 per group/condition)

2. **Result Structure Definition**
   ```julia
   struct Cluster
       id::Int
       electrodes::Vector{Symbol}
       time_indices::Vector{Int}
       time_range::Tuple{Float64, Float64}
       cluster_stat::Float64
       p_value::Float64
       is_significant::Bool
       polarity::Symbol  # :positive or :negative
   end
   
   struct ClusterPermutationResult
       # Test information
       test_type::Symbol  # :paired or :independent
       design::Symbol
       
       # Observed statistics
       t_matrix::Array{Float64, 2}  # [electrodes × time]
       df_matrix::Array{Float64, 2}  # Degrees of freedom per point
       significant_mask_positive::BitArray{2}  # [electrodes × time]
       significant_mask_negative::BitArray{2}  # [electrodes × time]
       
       # Clusters
       positive_clusters::Vector{Cluster}
       negative_clusters::Vector{Cluster}
       cluster_stats_positive::Vector{Float64}
       cluster_stats_negative::Vector{Float64}
       
       # Permutation info
       n_permutations::Int
       permutation_max_positive::Vector{Float64}  # Distribution
       permutation_max_negative::Vector{Float64}
       random_seed::Union{Int, Nothing}
       
       # Metadata
       electrodes::Vector{Symbol}
       time_points::Vector{Float64}
       threshold::Float64
       threshold_method::Symbol  # :parametric, :nonparametric_individual, :nonparametric_common
       cluster_type::Symbol  # :spatial, :temporal, :spatiotemporal
       cluster_statistic::Symbol  # :sum, :max, :size, :wcm
       min_cluster_size::Int
       critical_t_values::Array{Float64, 2}  # Critical t-values used
   end
   ```

### Phase 1: Core Statistics Computation
**Goal**: Compute t-statistics for all electrode × time points

1. **`compute_t_matrix` function**
   - Input: `PermutationTestData`
   - Output: `Array{Float64, 2}` [electrodes × time] of t-values, `Array{Float64, 2}` of df
   - Use existing `paired_ttest` or `independent_ttest` functions
   - Loop over electrodes and time points
   - Vectorize where possible for performance
   - Handle NaN/Inf edge cases

**Implementation approach:**
```julia
function compute_t_matrix(prepared::PermutationTestData)
    n_electrodes = length(prepared.electrodes)
    n_time = length(prepared.time_points)
    t_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    df_matrix = Array{Float64, 2}(undef, n_electrodes, n_time)
    
    for e_idx in 1:n_electrodes
        for t_idx in 1:n_time
            data_A = prepared.data_A[:, e_idx, t_idx]
            data_B = prepared.data_B[:, e_idx, t_idx]
            
            if prepared.design == :paired
                result = paired_ttest(data_A, data_B)
            else
                result = independent_ttest(data_A, data_B)
            end
            t_matrix[e_idx, t_idx] = result.t
            df_matrix[e_idx, t_idx] = result.df
        end
    end
    return t_matrix, df_matrix
end
```

2. **`compute_critical_t_values` function**
   - Input: df_matrix, alpha (default 0.05), tail
   - Output: `Array{Float64, 2}` [electrodes × time] of critical t-values
   - Use TDist to compute critical values
   - Handle two-tailed (split alpha/2) vs one-tailed

### Phase 2: Thresholding
**Goal**: Create binary mask of significant points

1. **`threshold_t_matrix_parametric` function**
   - Input: t_matrix, df_matrix, critical_t_values, tail
   - Output: `BitArray{2}` [electrodes × time] masks (separate for positive/negative)
   - Convert t-values to p-values using t-distribution
   - Threshold based on critical t-values
   - Handle two-tailed tests (both positive and negative separately)

2. **`threshold_t_matrix_nonparametric` function** (Optional, Phase 2b)
   - Input: t_matrix, permutation_t_matrices, alpha, tail
   - Output: `BitArray{2}` masks
   - Two variants:
     - `nonparametric_individual`: Each electrode/time gets its own threshold
     - `nonparametric_common`: Single threshold across all points
   - Compute percentiles from permutation distribution

**Options:**
- Parametric threshold: Use t-distribution to convert t to p, then threshold (default)
- Non-parametric threshold: Use permutation distribution (more robust, slower)

### Phase 3: Connectivity Matrix Construction
**Goal**: Build adjacency matrix for clustering

1. **`build_connectivity_matrix` function**
   - Input: electrodes, layout, cluster_type
   - Output: `SparseMatrixCSC{Bool}` connectivity matrix
   - Convert `Layout.neighbours` to adjacency matrix
   - For spatiotemporal: extend to [electrodes × time] space
   - Handle spatial, temporal, and spatiotemporal adjacency

**Implementation approach:**
```julia
function build_connectivity_matrix(
    electrodes::Vector{Symbol},
    layout::Layout,
    cluster_type::Symbol
) -> SparseMatrixCSC{Bool}
    # Build spatial connectivity from Layout.neighbours
    # For temporal: add connections between consecutive time points
    # For spatiotemporal: combine both
    # Return sparse matrix for efficiency
end
```

### Phase 4: Cluster Finding
**Goal**: Find connected clusters in thresholded data

1. **`find_clusters_connected_components` function** (Core algorithm)
   - Input: thresholded mask, connectivity matrix, min_cluster_size
   - Output: Vector of Cluster structs (without statistics)
   - Use BFS/DFS graph traversal to find connected components
   - Handle spatial, temporal, and spatiotemporal adjacency
   - Filter clusters by minimum size

**Algorithm details:**
- Start with first significant point
- Use BFS to find all connected significant points
- Label cluster with unique ID
- Move to next unlabeled significant point
- Repeat until all points labeled
- Filter out clusters smaller than min_cluster_size

**Spatial adjacency:**
- Use connectivity matrix built from `Layout.neighbours`
- Two electrodes are adjacent if connected in matrix

**Temporal adjacency:**
- Consecutive time points (i and i+1)
- Handled in connectivity matrix construction

**Spatiotemporal:**
- Point (electrode_i, time_j) is adjacent to:
  - All spatial neighbors of electrode_i at time_j
  - electrode_i at time_j±1
  - All handled via connectivity matrix

2. **`find_clusters` wrapper function**
   - Input: thresholded masks (positive/negative), connectivity, min_size, time_points
   - Output: Separate vectors of positive and negative clusters
   - Calls `find_clusters_connected_components` for each polarity
   - Adds time_range information to clusters

### Phase 5: Cluster Statistics
**Goal**: Compute cluster-level statistics

1. **`compute_cluster_statistics` function**
   - Input: clusters, t_matrix, statistic_type
   - Output: Updated clusters with cluster_stat
   - Support multiple statistic types:
     - `:sum` - Sum of t-values within cluster (default, FieldTrip `maxsum`)
     - `:max` - Maximum t-value in cluster
     - `:size` - Number of points in cluster (FieldTrip `maxsize`)
     - `:wcm` - Weighted cluster mass (FieldTrip `wcm`)
   - Handle positive and negative clusters separately

2. **`compute_cluster_stat_wcm` function** (Weighted Cluster Mass)
   - Input: cluster, t_matrix, weight (default 1.0)
   - Output: WCM statistic
   - Combines cluster size and intensity
   - Formula: WCM = sum(|t|^weight) or similar (see FieldTrip implementation)

### Phase 6: Permutation Loop
**Goal**: Generate permutation distribution

1. **`shuffle_labels` function**
   - Input: prepared data, design type, random number generator
   - Output: Shuffled PermutationTestData (or indices for shuffling)
   - **Paired design**: Randomly swap A/B labels within each participant
     - For each participant: randomly decide whether to swap conditions
     - Maintains within-subject structure
   - **Independent design**: Randomly shuffle participants between groups
     - Randomly assign participants to groups
     - Maintains group sizes (n_A and n_B)

2. **`run_permutations` function**
   - Input: prepared data, n_permutations, all threshold/cluster parameters, progress callback
   - Output: Vectors of max cluster statistics (positive/negative) from each permutation
   - **Progress reporting**: Use `ProgressMeter.jl` to show progress bar, ETA, current permutation
   - **Random seed**: Set seed at start if provided for reproducibility
   - Loop N times (default 1000):
     - Shuffle labels (using `shuffle_labels`)
     - Recompute t_matrix (using `compute_t_matrix`)
     - Re-threshold (using `threshold_t_matrix_parametric`)
     - Find clusters (using `find_clusters`)
     - Compute cluster statistics (using `compute_cluster_statistics`)
     - Store maximum cluster statistic (separate for positive/negative)
     - Update progress bar
   - **Edge case handling**: 
     - If no clusters found in permutation, store 0.0
     - Handle NaN/Inf values gracefully
   - **Performance**: Consider parallelization with `Threads.@threads` for large N

### Phase 7: Inference
**Goal**: Compute p-values and determine significance

1. **`compute_cluster_pvalues` function**
   - Input: observed clusters, permutation distributions (positive/negative)
   - Output: Updated clusters with p-values and is_significant flags
   - **For each observed cluster**:
     - Compare cluster statistic to permutation distribution
     - Handle positive and negative clusters separately
     - p = (number of permutations with max_cluster_stat >= observed) / N
     - Add 1 to numerator and denominator (standard correction)
     - Mark as significant if p < 0.05 (or specified alpha)
   - **Edge case**: If no clusters found in observed data, return empty result

2. **`assemble_results` function**
   - Input: All computed components
   - Output: Complete `ClusterPermutationResult` struct
   - Populate all fields
   - Handle edge cases (no clusters, empty results)

### Phase 8: Main Function
**Goal**: Put it all together

```julia
function cluster_permutation_test(
    prepared::PermutationTestData;
    n_permutations::Int = 1000,
    threshold::Float64 = 0.05,  # p-value threshold
    threshold_method::Symbol = :parametric,  # :parametric, :nonparametric_individual, :nonparametric_common
    cluster_type::Symbol = :spatiotemporal,  # :spatial, :temporal, :spatiotemporal
    cluster_statistic::Symbol = :sum,  # :sum, :max, :size, :wcm
    min_cluster_size::Int = 0,  # Minimum number of points in cluster
    tail::Symbol = :both,  # :both, :left, :right
    random_seed::Union{Int, Nothing} = nothing,  # For reproducibility
    show_progress::Bool = true  # Show progress bar
) -> ClusterPermutationResult
```

**Function flow:**
1. Validate inputs (Phase 0)
2. Ensure Layout.neighbours exists (compute if needed)
3. Build connectivity matrix (Phase 3)
4. Compute observed t-matrix and df-matrix (Phase 1)
5. Compute critical t-values (Phase 1)
6. Threshold observed data (Phase 2)
7. Find observed clusters (Phase 4)
8. Compute observed cluster statistics (Phase 5)
9. Run permutations (Phase 6)
10. Compute p-values (Phase 7)
11. Assemble and return results (Phase 7)

### Phase 9: Visualization (Optional, Lower Priority)
**Goal**: Provide plotting functions for results

1. **`plot_cluster_topography` function**
   - Plot significant clusters on electrode layout
   - Color-code by cluster
   - Show time window for each cluster
   - Overlay on head plot (using Makie)

2. **`plot_cluster_timecourse` function**
   - Plot ERP waveforms for each condition/group
   - Highlight significant time windows
   - Mark cluster boundaries
   - Show t-values over time

3. **`plot_spatiotemporal_map` function**
   - 2D plot: electrodes × time
   - Color-code by t-value
   - Outline significant clusters
   - Show cluster IDs and p-values

4. **`plot_permutation_distribution` function**
   - Histogram of permutation distribution
   - Mark observed cluster statistics
   - Show significance threshold

5. **`print_cluster_summary` function**
   - Pretty table of all clusters with:
     - Cluster ID
     - Electrodes
     - Time window
     - Cluster statistic
     - p-value
     - Significant (yes/no)

### Phase 10: Result Export (Optional, Lower Priority)
**Goal**: Save results for later analysis

1. **`save_permutation_results` function**
   - Save `ClusterPermutationResult` to JLD2 file
   - Include all necessary data for re-analysis
   - Save metadata for reproducibility (random_seed, parameters, etc.)

## Suggested Implementation Order

### Phase 0: Foundation (Week 1)
1. **Define result structures** - Cluster and ClusterPermutationResult
2. **Implement input validation** - Catch errors early
3. **Test with prepared data** - Ensure data structure works

### Phase 1-2: Core Statistics (Week 1-2)
4. **Implement `compute_t_matrix`** - Test with your prepared data
5. **Implement `compute_critical_t_values`** - For parametric thresholding
6. **Implement `threshold_t_matrix_parametric`** - Verify thresholding works
7. **Test with simple cases** - Verify t-values and thresholding are correct

### Phase 3-4: Clustering (Week 2-3)
8. **Implement `build_connectivity_matrix`** - Start with temporal only (simpler)
9. **Implement `find_clusters_connected_components`** - Core BFS/DFS algorithm
10. **Test temporal clustering** - Verify clusters are found correctly
11. **Add spatial connectivity** - Extend to spatial clustering
12. **Add spatiotemporal** - Combine both
13. **Test with different cluster types** - Verify all work

### Phase 5: Cluster Statistics (Week 3)
14. **Implement `compute_cluster_statistics`** - Simple sum first
15. **Add other statistic types** - max, size, wcm
16. **Test with different statistics** - Verify calculations

### Phase 6: Permutation (Week 3-4)
17. **Implement `shuffle_labels`** - Test with small N
18. **Implement `run_permutations`** - Start with 10-100 permutations for testing
19. **Add progress reporting** - Use ProgressMeter
20. **Add random seed control** - For reproducibility
21. **Test permutation loop** - Verify shuffling and statistics

### Phase 7: Inference (Week 4)
22. **Implement `compute_cluster_pvalues`** - Final inference step
23. **Implement `assemble_results`** - Put everything together
24. **Test complete pipeline** - End-to-end test

### Phase 8: Main Function (Week 4)
25. **Create main `cluster_permutation_test` function** - Integrate everything
26. **Add comprehensive error handling** - Edge cases
27. **Test with real data** - Your prepared data

### Phase 9-10: Polish (Week 5+)
28. **Add visualization functions** - As needed
29. **Add result export** - Save to JLD2
30. **Documentation and examples** - Comprehensive docs

## Key Design Decisions

1. **Threshold**: Use parametric threshold by default (p < 0.05) - convert t to p using t-distribution
   - Support non-parametric as optional (slower but more robust)
2. **Cluster statistic**: Sum of t-values (FieldTrip default `maxsum`)
   - Support alternatives: max, size, weighted cluster mass
3. **Permutations**: Default 1000, make configurable
   - Show progress for user feedback
4. **Positive/negative**: Handle separately (two-tailed test)
   - Separate permutation distributions
   - Separate p-value computation
5. **Spatial clustering**: Check Layout.neighbours exists, compute if missing
   - Use default distance criterion if not specified
6. **Minimum cluster size**: Default 0 (no filtering)
   - Make configurable via `min_cluster_size` parameter
7. **Random seed**: Optional, for reproducibility
   - Set at start of permutation loop
8. **Progress reporting**: Always show for permutation loop
   - Use ProgressMeter.jl (already in dependencies)
9. **Performance**: 
   - Use sparse matrices for connectivity
   - Consider parallelization for large permutation counts
   - Vectorize where possible

## Testing Strategy

1. **Unit tests**: Test each component independently
   - `compute_t_matrix` with known data
   - `threshold_t_matrix` with simple cases
   - `find_clusters` with manual test cases
2. **Integration tests**: Test components together
   - Thresholding + clustering
   - Permutation loop end-to-end
3. **Validation tests**: Compare to FieldTrip output if possible
   - Use same data and parameters
   - Verify similar results
4. **Edge case tests**:
   - No clusters found
   - All points significant
   - Single participant/electrode
   - NaN/Inf values in data
   - Empty clusters after filtering
5. **Performance tests**: 
   - Test with small datasets first
   - Scale up gradually
   - Measure timing for optimization

## Error Handling

1. **Input validation**: Check all parameters at start
2. **Data validation**: Check data dimensions, missing values
3. **Layout validation**: Ensure neighbours exist or can be computed
4. **Edge cases**: Handle gracefully with informative messages
5. **Progress reporting**: Show errors clearly in progress output

## Next Steps

1. **Start with Phase 0**: Define structures and validation
2. **Implement Phase 1**: `compute_t_matrix` and critical values
3. **Test with prepared data**: Use `test_permutation_data.jl`
4. **Iterate**: Build and test each component incrementally
5. **Validate**: Compare with FieldTrip when possible

