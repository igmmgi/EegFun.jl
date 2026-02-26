This demo shows statistical testing and topography visualization for time-frequency data, including analytic tests and cluster-based permutation tests.

### Key Functions

| Function | Purpose |
| --- | --- |
| `tf_morlet` | Compute TF decomposition (batch, from epoch files) |
| `prepare_stats` | Prepare multi-participant TF data for statistical comparison |
| `analytic_test` | Parametric test (fast, no permutations) |
| `permutation_test` | Cluster-based permutation test (robust) |
| `plot_tf` | Plot TF power for a single channel |
| `plot_tf_stats` | Per-channel TF heatmap with significance |
| `plot_topography` | TF topography for a frequency/time window |
| `plot_topo_stats` | Topography with statistical significance overlay |

## Workflow Summary

### TF Decomposition

- Batch compute Morlet wavelets from epoch data with `tf_morlet`

### Single-Participant Visualization

- `plot_tf` for per-channel time-frequency plots with baseline correction
- `plot_topography` for scalp maps at specific frequency bands and time windows

### Group Statistics

- `prepare_stats` loads all participants and prepares for comparison (paired or independent)
- Supports frequency selection, time intervals, and baseline correction
- `analytic_test` for quick exploration
- `permutation_test` with spatiotemporal clustering for robust inference

### Statistical Visualization

- `plot_tf_stats` shows per-channel significance as contours
- `plot_topo_stats` shows scalp topography with significance highlighting for specific frequency bands
