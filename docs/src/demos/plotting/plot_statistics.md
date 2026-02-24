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
# 1. LOAD DATA AND CREATE ERPS
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

#######################################################################
# 2. RUN STATISTICAL TEST
#######################################################################

# prepare for comparison and run analytic test
prepared = EegFun.prepare_statistics(epochs)
result = EegFun.analytic_test(prepared, correction_method = :bonferroni)

#######################################################################
# 3. BASIC PLOT — ERP WAVEFORMS (single channel)
#######################################################################

# plot ERP waveforms for both conditions at a channel
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true)

#######################################################################
# 4. ADD DIFFERENCE WAVE
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true)

#######################################################################
# 5. ADD SIGNIFICANCE MARKERS
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true, plot_significance = true)

#######################################################################
# 6. SHOW T-VALUES
#######################################################################

EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Pz), plot_tvalues = true, plot_critical_t = true)

#######################################################################
# 7. SHIFT DIFFERENCE WAVE
#######################################################################

# offset the difference wave for visibility (similar to MATLAB)
EegFun.plot_erp_stats(result, channel_selection = EegFun.channels(:Cz), plot_erp = true, plot_difference = true, difference_offset = 3.0)

#######################################################################
# 8. SIGNIFICANCE BAR POSITIONING
#######################################################################

# position significance bars at the bottom, at zero, or at a custom y
EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels(:Cz),
    plot_erp = true,
    plot_significance = true,
    significance_position = :bottom,
    significance_color = (:red, 0.5),
)

#######################################################################
# 9. GRID LAYOUT — multiple channels
#######################################################################

EegFun.plot_erp_stats(
    result,
    channel_selection = EegFun.channels([:Cz, :Pz, :Oz, :Fz]),
    layout = :grid,
    plot_erp = true,
    plot_significance = true,
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
