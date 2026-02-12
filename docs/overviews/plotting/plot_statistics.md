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
