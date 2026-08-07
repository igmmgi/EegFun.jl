# ERP Measurements

This demo shows how to extract quantitative measurements from ERPs and analyse them with traditional statistics using AnovaFun.

## Key Functions

| Function | Purpose |
| --- | --- |
| `erp_measurements` | Batch extract amplitude/latency/area from all participant ERPs → CSV |
| `plot_erp_measurement_gui` | Interactive GUI to explore measurements before batch extraction |
| `paired_ttest` | AnovaFun: paired t-test between two conditions |
| `independent_ttest` | AnovaFun: independent samples t-test |
| `anova` | AnovaFun: repeated-measures or between-subject ANOVA |
| `emmeans` / `pairwise` | AnovaFun: estimated marginal means and pairwise contrasts |

## Available Measurement Types

| Category | Types |
| --- | --- |
| **Amplitude** | `mean_amplitude`, `max_peak_amplitude`, `min_peak_amplitude`, `peak_to_peak_amplitude` |
| **Latency** | `max_peak_latency`, `min_peak_latency`, `peak_to_peak_latency`, `fractional_area_latency`, `fractional_peak_latency` |
| **Area** | `rectified_area`, `integral`, `positive_area`, `negative_area` |

## Workflow Summary

## Extract Measurements

- Use `plot_erp_measurement_gui` to explore and choose parameters interactively
- Run `erp_measurements` to batch-extract values across all participants
- Output is a DataFrame and a saved CSV with columns: participant, condition, channel, measurement

## Traditional Statistics with AnovaFun

- Load the CSV or use the in-memory DataFrame
- Use `paired_ttest` for condition comparisons at specific channels
- Use `anova` for multi-factor within/between-subject designs
- Follow up with `emmeans` and `pairwise` for post-hoc contrasts


## Code Examples

::: details Show Code

```julia
# Demo: ERP Measurements and Traditional Statistics
# Shows how to extract ERP measurements (amplitudes, latencies) and 
# analyse them using AnovaFun for traditional t-tests and ANOVA.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie
using AnovaFun

#######################################################################
# EXTRACT ERP MEASUREMENTS
#######################################################################

# This batch function loads all participant ERP files, extracts 
# measurements, and saves results to a CSV file.

input_dir = EegFun.example_path("data/julia/erps")
file_pattern = "erps_final"

# Mean amplitude in the P300 window
mean_amp = EegFun.erp_measurements(
    file_pattern,
    "mean_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    channel_selection = EegFun.channels([:Pz, :Cz, :Fz]),
    analysis_interval = (0.3, 0.5),
    baseline_interval = (-0.2, 0.0),
)

# Result is an ErpMeasurementsResult containing a DataFrame
# It's also saved as a CSV to the output directory
mean_amp.data  # DataFrame with: participant, condition, channel, measurement


#######################################################################
# AVAILABLE MEASUREMENT TYPES
#######################################################################

# Amplitude measurements
# "mean_amplitude"          — mean voltage in interval
# "max_peak_amplitude"      — maximum peak value (robust detection)
# "min_peak_amplitude"      — minimum peak value (robust detection)
# "peak_to_peak_amplitude"  — difference between max and min peaks

# Latency measurements
# "max_peak_latency"        — time of maximum peak
# "min_peak_latency"        — time of minimum peak
# "peak_to_peak_latency"    — time difference between max and min peaks
# "fractional_area_latency" — time point dividing area into fraction
# "fractional_peak_latency" — time where amplitude is fraction of peak

# Area measurements
# "rectified_area"          — sum of absolute voltage values
# "integral"                — signed area (positive minus negative)
# "positive_area"           — area of positive deflections only
# "negative_area"           — area of negative deflections only


#######################################################################
# EXPLORE MEASUREMENTS INTERACTIVELY
#######################################################################

# GUI for exploring measurements interactively before batch extraction
dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))
EegFun.plot_erp_measurement_gui(dat)     # all conditions
# EegFun.plot_erp_measurement_gui(dat[1]) # first condition only


#######################################################################
# PEAK AMPLITUDE AND LATENCY
#######################################################################

# Max peak with robust detection (requires peak to be larger than neighbors)
max_peak = EegFun.erp_measurements(
    file_pattern,
    "max_peak_amplitude",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = (0.3, 0.6),
    baseline_interval = (-0.2, 0.0),
    local_interval = 3,  # peak must be larger than 3 neighbors on each side
)

# Corresponding latency
max_latency = EegFun.erp_measurements(
    file_pattern,
    "max_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    channel_selection = EegFun.channels([:Pz]),
    analysis_interval = (0.3, 0.6),
    baseline_interval = (-0.2, 0.0),
)


#######################################################################
# TRADITIONAL STATISTICS WITH ANOVAFUN
#######################################################################

# The CSV output from erp_measurements can be used directly with AnovaFun.
# Here we demonstrate the workflow in-memory using the DataFrame.

# --- Paired t-test: condition 1 vs condition 2 at Pz ---
df = mean_amp.data

# Channels are in wide format (each channel is a column, e.g. :Pz, :Cz, :Fz)
cond1_pz = df[df.condition.==1, :Pz]
cond2_pz = df[df.condition.==2, :Pz]
result = paired_ttest(cond1_pz, cond2_pz)
result.t   # t-statistic
result.p   # p-value

# --- Independent t-test ---
# result = independent_ttest(group1_data, group2_data)

# --- Repeated-measures ANOVA ---
# For a 2-factor within-subject ANOVA (condition × channel), reshape your 
# data into the format AnovaFun expects and call:
# result = anova(data, within = [:condition, :channel])
# anova_table(result)
# emmeans(result, :condition)
# pairwise(result, :condition)
```

:::

## See Also

- [API Reference](../../reference/index.md)
