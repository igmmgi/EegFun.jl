# TF Statistics

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
| `plot_topography_stats` | Topography with statistical significance overlay |

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
- `plot_topography_stats` shows scalp topography with significance highlighting for specific frequency bands


## Code Examples

::: details Show Code

```julia
# Demo: Time-Frequency Statistics & Topography Plots
# Shows TF topography plotting, statistical tests, and
# visualization using plot_topography_stats for TF data.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# ==============================================================================
# Regenerate TF morlet data from epoch data 
# ==============================================================================
epoch_dir = joinpath(homedir(), "Documents", "Julia", "TestDataSets", "AttentionExp", "output_data")  # ← Update to your data path
println("Regenerating TF morlet data from epochs...")

EegFun.tf_morlet(
    "epochs_final";
    input_dir = epoch_dir,
    condition_selection = EegFun.conditions([1, 2]),
    frequencies = 1:2:40,
    cycles = 7,
    pad = :both,
)

# ==============================================================================
# Setup
# ==============================================================================
input_dir = joinpath(epoch_dir, "tf_morlet_epochs_final")
file_pattern = "epochs_final"

# ==============================================================================
# Load a single participant to test plot_topography for TF data
# ==============================================================================
println("Loading single participant TF data...")
single_data = EegFun.read_data(joinpath(input_dir, "example1_epochs_final.jld2"))

# Inspect conditions
println("Number of datasets: ", length(single_data))
for d in single_data
    println("  Condition $(d.condition): $(d.condition_name), file=$(d.file)")
end

# Pick the first condition for topo testing
tf_cond1 = single_data[1]
println("\nFrequencies: ", sort(unique(tf_cond1.data_power.freq)))
println("Time range: ", extrema(tf_cond1.data_power.time))
println("Channels: ", EegFun.channel_labels(tf_cond1))

# ==============================================================================
# Test plot_topography for single TimeFreqData
# ==============================================================================
println("\n--- Testing plot_topography (single TF) ---")

EegFun.plot_tf(tf_cond1, baseline_interval = (-0.5, -0.3), baseline_method = :db, channel_selection = EegFun.channels(:Fp1))

# Alpha band (8-12 Hz) in the 200-400 ms window, with dB baseline
EegFun.plot_topography(
    tf_cond1,
    freq_range = (8.0, 20.0),
    interval_selection = EegFun.times(0.2, 0.4),
    baseline_interval = (-0.5, -0.3),
    baseline_method = :percent,
)

# ==============================================================================
# Test plot_topography for Vector{TimeFreqData} (multiple conditions)
# ==============================================================================
println("\n--- Testing plot_topography (multi-condition) ---")

EegFun.plot_topography(single_data, freq_range = (8.0, 12.0), interval_selection = EegFun.times(0.2, 0.4), baseline_interval = (-0.3, 0.0))

# ==============================================================================
# Prepare TF statistics (all participants, 2 conditions)
# ==============================================================================
println("\n--- Preparing TF statistics ---")

tf_stat_data = EegFun.prepare_stats(
    file_pattern,
    :paired;
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    frequency_selection = (2.0, 40.0),          # 2-40 Hz
    interval_selection = (0.0, 1.0),            # 0-1000 ms
    baseline_interval = (-0.5, -0.3),           # Baseline: -500 to -300 ms
    baseline_method = :db,
)
println(tf_stat_data)

# ==============================================================================
# Analytic test (fast, no permutations)
# ==============================================================================
println("\n--- Running analytic test ---")
result_analytic = EegFun.analytic_test(tf_stat_data)
println(result_analytic)

# ==============================================================================
# Plot TF stats (existing per-channel heatmap)
# ==============================================================================
println("\n--- Plotting TF stats (heatmap) ---")

EegFun.plot_tf_stats(result_analytic, channel_selection = EegFun.channels(:Cz), significance = :contour)

# ==============================================================================
# Plot TF topo stats (NEW — topography with significance)
# ==============================================================================
println("\n--- Plotting TF topo stats (alpha band) ---")

EegFun.plot_topography_stats(result_analytic, freq_range = (8.0, 12.0), n_topos = 10, highlight_threshold = 0.5)

# Theta band
println("\n--- Plotting TF topo stats (theta band) ---")

EegFun.plot_topography_stats(result_analytic, freq_range = (4.0, 7.0), n_topos = 8, highlight_threshold = 0.3, topo_data = :difference)

# ==============================================================================
# Cluster permutation test (slower but more robust)
# ==============================================================================
println("\n--- Running cluster permutation test ---")

result_perm = EegFun.permutation_test(
    tf_stat_data,
    n_permutations = 500,
    threshold_method = :parametric,
    cluster_type = :spatiotemporal,
    min_num_neighbors = 3,
    show_progress = true,
)
println(result_perm)

# ==============================================================================
# Plot cluster permutation TF topo stats
# ==============================================================================
println("\n--- Plotting permutation TF topo stats ---")

EegFun.plot_topography_stats(result_perm, freq_range = (8.0, 12.0), n_topos = 10, highlight_threshold = 0.5)

println("\n✓ All TF stats tests completed!")
```

:::

## See Also

- [API Reference](../../reference/index.md)
