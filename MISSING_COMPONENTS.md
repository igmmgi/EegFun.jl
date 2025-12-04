# Missing Components in Permutation Test Implementation Plan

## Critical Missing Components

### 1. **Connectivity Matrix Construction**
**Status**: Not mentioned in plan
**Why needed**: FieldTrip uses a connectivity matrix (sparse or dense) for clustering, not just the neighbours dict directly
**Implementation needed**:
```julia
function build_connectivity_matrix(electrodes::Vector{Symbol}, layout::Layout) -> SparseMatrixCSC
    # Convert Layout.neighbours to adjacency matrix
    # For spatiotemporal: need to extend to [electrodes × time] space
end
```

### 2. **Graph Traversal Algorithm (findcluster equivalent)**
**Status**: Mentioned but not detailed
**Why needed**: Core algorithm for finding connected components
**Implementation needed**:
```julia
function find_clusters_connected_components(
    mask::BitArray{2},  # [electrodes × time]
    connectivity::SparseMatrixCSC,
    min_size::Int = 0
) -> Vector{Cluster}
    # BFS/DFS to find connected components
    # Handle spatial, temporal, and spatiotemporal adjacency
end
```

### 3. **Degrees of Freedom Handling for Thresholding**
**Status**: Not addressed
**Why needed**: For parametric thresholding, need critical t-values which depend on df
- Paired design: df = n - 1 (same for all points)
- Independent design: df = n_A + n_B - 2 (same for all points)
- But if missing data: df can vary per electrode/time
**Implementation needed**:
```julia
function compute_critical_t_values(
    prepared::PermutationTestData,
    alpha::Float64 = 0.05
) -> Array{Float64, 2}  # [electrodes × time] of critical t-values
```

### 4. **Minimum Cluster Size Parameter**
**Status**: Not mentioned
**Why needed**: FieldTrip has `minnbchan` - filters out clusters smaller than threshold
**Implementation needed**:
- Add `min_cluster_size::Int = 0` parameter
- Filter clusters after finding them

### 5. **Non-Parametric Thresholding Option**
**Status**: Not mentioned
**Why needed**: FieldTrip supports `nonparametric_individual` and `nonparametric_common` thresholding
**Implementation needed**:
```julia
function threshold_nonparametric(
    t_matrix::Array{Float64, 2},
    permutation_t_matrices::Array{Float64, 3},  # [electrodes × time × permutations]
    alpha::Float64 = 0.05
) -> BitArray{2}
```

### 6. **Progress Reporting**
**Status**: Not mentioned
**Why needed**: Permutation loops can take hours - users need feedback
**Implementation needed**:
- Use `ProgressMeter.jl` (already in dependencies)
- Show progress bar, ETA, current permutation number

### 7. **Random Seed Control**
**Status**: Not mentioned
**Why needed**: Reproducibility for scientific results
**Implementation needed**:
- Add `random_seed::Union{Int, Nothing} = nothing` parameter
- Set seed at start of permutation loop

### 8. **Edge Case Handling**
**Status**: Not detailed
**Why needed**: Robust implementation needs to handle:
- No clusters found in observed data
- All NaN values in t_matrix
- Empty clusters after filtering
- Single participant/electrode edge cases
**Implementation needed**:
- Check for empty results
- Return appropriate empty structures
- Warn user appropriately

### 9. **Complete Result Structure Definition**
**Status**: Partially defined
**Why needed**: Need full structure with all fields
**Missing fields**:
```julia
struct ClusterPermutationResult
    # ... existing fields ...
    
    # Missing:
    df_matrix::Array{Float64, 2}  # Degrees of freedom per point
    critical_t_values::Array{Float64, 2}  # Critical t-values used
    positive_clusters::Vector{Cluster}
    negative_clusters::Vector{Cluster}
    max_cluster_stat_positive::Float64  # From observed data
    max_cluster_stat_negative::Float64
    permutation_max_positive::Vector{Float64}  # Distribution
    permutation_max_negative::Vector{Float64}
    threshold_method::Symbol  # :parametric, :nonparametric_individual, etc.
    cluster_type::Symbol  # :spatial, :temporal, :spatiotemporal
    min_cluster_size::Int
    random_seed::Union{Int, Nothing}
end
```

### 10. **Input Validation**
**Status**: Not mentioned
**Why needed**: Catch errors early with helpful messages
**Validation needed**:
- Check Layout.neighbours exists (or compute if missing)
- Validate n_permutations > 0
- Validate threshold in valid range (0 < threshold < 1)
- Check data dimensions match
- Validate design type

### 11. **Weighted Cluster Mass (WCM) Option**
**Status**: Not mentioned
**Why needed**: FieldTrip supports this alternative cluster statistic
**Implementation needed**:
```julia
function compute_cluster_stat_wcm(
    cluster::Cluster,
    t_matrix::Array{Float64, 2},
    weight::Float64 = 1.0
) -> Float64
    # Weighted combination of cluster size and intensity
end
```

### 12. **Separate Positive/Negative Cluster Handling**
**Status**: Mentioned but not detailed
**Why needed**: Two-tailed tests need separate handling
**Implementation needed**:
- Find positive and negative clusters separately
- Compute statistics separately
- Compare to separate permutation distributions
- Return both in result structure

### 13. **Performance Optimizations**
**Status**: Not mentioned
**Why needed**: Permutation tests are computationally intensive
**Considerations**:
- Vectorize t-matrix computation where possible
- Use sparse matrices for connectivity
- Parallelize permutation loop (Threads.@threads)
- Pre-allocate arrays
- Consider GPU acceleration for large datasets (future)

### 14. **Visualization Functions**
**Status**: Mentioned in docs but not in plan
**Why needed**: Users need to visualize results
**Functions needed**:
```julia
function plot_cluster_topography(result::ClusterPermutationResult, cluster_id::Int)
function plot_cluster_timecourse(result::ClusterPermutationResult, cluster_id::Int)
function plot_spatiotemporal_map(result::ClusterPermutationResult)
function plot_permutation_distribution(result::ClusterPermutationResult)
function print_cluster_summary(result::ClusterPermutationResult)
```

### 15. **Result Export/Saving**
**Status**: Not mentioned
**Why needed**: Users may want to save results for later analysis
**Implementation needed**:
- Save to JLD2 format
- Include all necessary data for re-analysis
- Metadata for reproducibility

### 16. **Documentation and Examples**
**Status**: Not mentioned
**Why needed**: Users need to understand how to use the function
**Needed**:
- Comprehensive docstrings
- Usage examples
- Comparison with FieldTrip workflow
- Troubleshooting guide

## Implementation Priority

### High Priority (Core Functionality)
1. Connectivity matrix construction
2. Graph traversal algorithm (find_clusters)
3. Degrees of freedom handling
4. Complete result structure
5. Input validation
6. Edge case handling

### Medium Priority (Important Features)
7. Progress reporting
8. Random seed control
9. Minimum cluster size
10. Separate positive/negative handling
11. Performance optimizations

### Lower Priority (Nice to Have)
12. Non-parametric thresholding
13. Weighted cluster mass
14. Visualization functions
15. Result export
16. Advanced documentation

## Recommended Next Steps

1. **Start with connectivity matrix** - Foundation for clustering
2. **Implement graph traversal** - Core clustering algorithm
3. **Add df handling** - Needed for proper thresholding
4. **Complete result structure** - Define all fields
5. **Add validation** - Catch errors early
6. **Test with simple cases** - Temporal-only first, then add spatial

