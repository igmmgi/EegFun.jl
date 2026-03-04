# Plot Statistics

This demo shows how to visualise ERP statistical test results, including condition comparisons, significance markers, and t-values.

### Why Plot Statistics?

- **Identify significant differences** — see where conditions diverge with corrected significance bars
- **Inspect effect size** — overlay difference waves and t-statistics
- **Compare methods** — works with both analytic (t-test) and permutation test results

### Key Functions

| Function | Purpose |
| --- | --- |
| `plot_erp_stats(result, channel_selection=channels(:Cz))` | Plot analytic test results |
| `plot_erp_stats(perm_result, channel_selection=channels(:Cz))` | Plot permutation test results |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `layout` | `:single` | Layout: `:single`, `:grid`, or `:topo` |
| `channel_selection` | `channels()` | Channel predicate (e.g., `channels(:Cz)`) |
| `plot_erp` | `true` | Show condition ERP waveforms |
| `plot_difference` | `false` | Show difference wave (A − B) |
| `plot_tvalues` | `false` | Overlay t-statistic curve |
| `plot_significance` | `false` | Add significance bars |
| `plot_critical_t` | `false` | Show critical t-value lines |
| `difference_offset` | `0.0` | Vertical offset for difference wave (μV) |
| `significance_position` | `:auto` | `:auto`, `:bottom`, `:zero`, or a `Float64` |
| `significance_color` | `(:gray, 0.6)` | Colour for significance bars |

### What You'll Learn

1. Plotting ERP waveforms with significance markers
2. Adding difference waves and t-value overlays
3. Controlling significance bar position and colour
4. Working with both analytic and permutation test results
5. Using layout options for multi-channel views


## Code Examples

::: details Show Code

```julia
# Demo: Plotting Statistical Test Results
# Shows how to visualise analytic and permutation test results with
# ERP waveforms, difference waves, t-values, and significance markers.

using EegFun

#######################################################################
# PREPARE DATA (requires pre-saved per-participant ERPs)
#######################################################################

input_dir = "./resources/data/julia/erps"
file_pattern = "erps_good"

prepared = EegFun.prepare_stats(
    file_pattern,
    :paired;
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    baseline_interval = EegFun.times((-0.2, 0.0)),
    analysis_interval = EegFun.times((0.0, 2.0)),
)

#######################################################################
# RUN ANALYTIC TEST
#######################################################################

result = EegFun.analytic_test(prepared, correction_method = :no)

#######################################################################
# 1. BASIC — ERP waveforms (single channel)
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:PO7))

#######################################################################
# 2. ADD DIFFERENCE WAVE
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_difference = true)

#######################################################################
# 3. ADD SIGNIFICANCE MARKERS
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:PO7), plot_difference = true, plot_significance = true)

#######################################################################
# 4. T-VALUES ONLY (no ERPs)
#######################################################################

EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels(:PO8),
    plot_erp = false,
    plot_tvalues = true,
    plot_critical_t = true,
    plot_significance = true,
)

#######################################################################
# 5. T-VALUES + ERPs TOGETHER
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Pz), plot_tvalues = true, plot_critical_t = true)

#######################################################################
# 6. SHIFT DIFFERENCE WAVE (for visibility)
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_difference = true, difference_offset = 3.0)

#######################################################################
# 7. SIGNIFICANCE BAR POSITIONING
#######################################################################

# :bottom
EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels(:Cz),
    plot_significance = true,
    significance_position = :bottom,
    significance_color = (:red, 0.5),
)

# :zero
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_significance = true, significance_position = :zero)

# custom y-position
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_significance = true, significance_position = -5.0)

#######################################################################
# 8. CUSTOM FIGURE TITLE
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_significance = true, figure_title = "N2pc Analytic Test")

#######################################################################
# 9. GRID LAYOUT — multiple channels
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels([:Cz, :Pz, :Oz, :Fz]), layout = :grid, plot_significance = true)

#######################################################################
# 10. TOPO LAYOUT — channels at electrode positions
#######################################################################

EegFun.plot_erp_stats(result, layout = :topo, plot_significance = true)

#######################################################################
# 11. CAPTURE RETURN VALUE for further customisation
#######################################################################

(; fig, axes) = EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels(:Cz),
    plot_difference = true,
    plot_significance = true,
    display_plot = false,
)
# axes[1] is a regular Makie Axis — add annotations, save, etc.

#######################################################################
# 12. SEM BANDS (±1 SEM around waveforms)
#######################################################################

# ±1 SEM
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:PO7), plot_se = true)

# SEM + difference wave
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:PO7), plot_se = true, plot_difference = true)

# SEM + significance markers + t-values
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:PO7), plot_se = true, plot_significance = true, plot_tvalues = true)

#######################################################################
# RUN PERMUTATION TEST
#######################################################################

result_perm = EegFun.permutation_test(prepared, n_permutations = 1000)

#######################################################################
# 12. PERMUTATION RESULTS — ERPs + significance
#######################################################################

EegFun.plot_erp_stats(result_perm, channel_selection = EegFun.channels(:Cz), plot_significance = true)

#######################################################################
# 13. PERMUTATION RESULTS — t-values + critical t
#######################################################################

EegFun.plot_erp_stats(result_perm, channel_selection = EegFun.channels(:Cz), plot_erp = false, plot_tvalues = true, plot_critical_t = true)
```

:::

## See Also

- [API Reference](../../reference/index.md)
