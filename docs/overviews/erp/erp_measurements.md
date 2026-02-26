This demo shows how to extract quantitative measurements from ERPs and analyse them with traditional statistics using AnovaFun.

### Key Functions

| Function | Purpose |
| --- | --- |
| `erp_measurements` | Batch extract amplitude/latency/area from all participant ERPs → CSV |
| `plot_erp_measurement_gui` | Interactive GUI to explore measurements before batch extraction |
| `paired_ttest` | AnovaFun: paired t-test between two conditions |
| `independent_ttest` | AnovaFun: independent samples t-test |
| `anova` | AnovaFun: repeated-measures or between-subject ANOVA |
| `emmeans` / `pairwise` | AnovaFun: estimated marginal means and pairwise contrasts |

### Available Measurement Types

| Category | Types |
| --- | --- |
| **Amplitude** | `mean_amplitude`, `max_peak_amplitude`, `min_peak_amplitude`, `peak_to_peak_amplitude` |
| **Latency** | `max_peak_latency`, `min_peak_latency`, `peak_to_peak_latency`, `fractional_area_latency`, `fractional_peak_latency` |
| **Area** | `rectified_area`, `integral`, `positive_area`, `negative_area` |

## Workflow Summary

### Extract Measurements

- Use `plot_erp_measurement_gui` to explore and choose parameters interactively
- Run `erp_measurements` to batch-extract values across all participants
- Output is a DataFrame and a saved CSV with columns: participant, condition, channel, measurement

### Traditional Statistics with AnovaFun

- Load the CSV or use the in-memory DataFrame
- Use `paired_ttest` for condition comparisons at specific channels
- Use `anova` for multi-factor within/between-subject designs
- Follow up with `emmeans` and `pairwise` for post-hoc contrasts
