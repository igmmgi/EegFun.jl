This demo shows how to visualise ERP statistical test results, including condition comparisons, significance markers, and t-values.

### Why Plot Statistics?

- **Identify significant differences** — see where conditions diverge with corrected significance bars
- **Inspect effect size** — overlay difference waves and t-statistics
- **Compare methods** — works with both analytic (t-test) and permutation test results

### Key Functions

| Function | Purpose |
| --- | --- |
| `plot_analytic_test(result, channel=:Cz)` | Plot analytic test results |
| `plot_analytic_test(perm_result, channel=:Cz)` | Plot permutation test results |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `plot_erp` | `false` | Show condition ERP waveforms |
| `plot_difference` | `false` | Show difference wave (A − B) |
| `plot_tvalues` | `false` | Overlay t-statistic curve |
| `show_significance` | `false` | Add significance bars |
| `show_critical_t` | `false` | Show critical t-value lines |
| `shift_difference` | `false` | Offset difference wave for visibility |
| `sig_bar_position` | `:auto` | `:auto`, `:bottom`, `:zero`, or a `Float64` |
| `sig_bar_color` | `(:gray, 0.6)` | Colour for significance bars |

### What You'll Learn

1. Plotting ERP waveforms with significance markers
2. Adding difference waves and t-value overlays
3. Controlling significance bar position and colour
4. Working with both analytic and permutation test results
