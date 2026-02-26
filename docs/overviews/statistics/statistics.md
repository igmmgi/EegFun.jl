This demo shows permutation-based statistical testing for ERP data, with multiple thresholding and cluster correction approaches.

### Key Functions

| Function | Purpose |
| --- | --- |
| `prepare_stats` | Load multi-participant ERPs and prepare for statistical comparison |
| `analytic_test` | Parametric t-test (uncorrected or Bonferroni) |
| `permutation_test` | Cluster-based permutation test with multiple thresholding options |
| `plot_erp_stats` | Visualise results with significance shading and critical t-values |

### Thresholding Methods

| Method | Description |
| --- | --- |
| `:parametric` | T-distribution threshold (fastest, default) |
| `:nonparametric_common` | Single threshold from pooled permutation distribution |
| `:nonparametric_individual` | Point-specific thresholds from permutation distributions |

### Cluster-Based Permutation Testing

Addresses the multiple comparisons problem in ERP analysis:

1. Compute test statistic at each time point × channel
2. Identify clusters of contiguous spatiotemporal significant effects
3. Permute condition labels and repeat (e.g., 1000 times)
4. Compare observed cluster mass to permutation distribution

**Key advantages:**
- Controls family-wise error rate without being overly conservative
- Sensitive to spatially and temporally distributed effects

## Workflow Summary

### Prepare Data

- Load group ERP data with `prepare_stats` specifying conditions, channels, and time intervals
- Supports `:paired` (within-subject) and `:independent` (between-subject) designs

### Analytic Tests

- Quick uncorrected or Bonferroni-corrected t-tests for exploration

### Permutation Tests

- Choose thresholding method based on assumptions and computational budget
- Set `cluster_type = :spatiotemporal` and `min_num_neighbors` for spatial constraints
- Visualise results with `plot_erp_stats`
