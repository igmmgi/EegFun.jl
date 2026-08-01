# Jackknife Average

This demo shows how to create jackknife (leave-one-out) averages for robust statistical testing of ERP latency measures.

## What is Jackknife Averaging?

The jackknife technique creates N grand averages from N participants, where each average excludes one participant. This produces smoother waveforms with clearly defined peaks, enabling reliable measurement of onset or peak latencies that would be noisy in individual-participant data.

## Why Use It?

Standard ERP peak latency measures are unreliable in noisy single-participant waveforms. The jackknife method:

- Produces cleaner waveforms by averaging many participants
- Enables latency measurement on smooth grand-average-like waveforms
- Requires a simple correction to the t-statistic: t_corrected = t / (n − 1)

## Key Functions

| Function | Purpose |
| --- | --- |
| `jackknife_average(erps)` | Single-condition jackknife from a vector of ERPs |
| `jackknife_average(file_pattern)` | Batch jackknife across participant files |

## Reference

Miller, Patterson, & Ulrich (1998). Jackknife-based method for measuring LRP onset latency differences. *Psychophysiology*, 35, 99–115.

## Workflow Summary

## Single-Condition Jackknife

- Create leave-one-out averages from a vector of participant ERPs

## Batch Processing

- Process all participant files and save jackknife averages

## Typical Pipeline

- Calculate LRP → Jackknife average → Measure onset latency → Apply correction to t-test


## Code Examples

::: details Show Code

```julia
# Demo: Jackknife Averaging
# Shows how to create leave-one-out (jackknife) averages for statistical
# testing of ERP latency measures.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# LOAD MULTI-PARTICIPANT ERP DATA
#######################################################################

# We need multiple participants for jackknifing
erps_p1 = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_final.jld2"))
erps_p2 = EegFun.read_data(EegFun.example_path("data/julia/erps/example2_erps_final.jld2"))
erps_p3 = EegFun.read_data(EegFun.example_path("data/julia/erps/example3_erps_final.jld2"))

# Create a list of one condition across participants (e.g., Condition 1)
cond_1_erps = [erps_p1[1], erps_p2[1], erps_p3[1]]

#######################################################################
# JACKKNIFE AVERAGING (IN-MEMORY)
#######################################################################

# Calculate jackknife averages: For N participants, this creates N averages 
# where the i-th average includes all participants EXCEPT participant i.
jackknife_erps = EegFun.jackknife_average(cond_1_erps)

# Plot the jackknife averages
EegFun.plot_erp(jackknife_erps, channel_selection = EegFun.channels(:Pz))


#######################################################################
# BATCH JACKKNIFE AVERAGING (FROM DISK)
#######################################################################

# Averages all JLD2 files whose name contains "erps_final" in the given directory.
# This automatically handles multiple conditions and saves the leave-one-out 
# datasets for each participant.
EegFun.jackknife_average("erps_final", input_dir = EegFun.example_path("data/julia/erps/"), output_dir = DEMO_OUTPUT)
```

:::

## See Also

- [API Reference](../../reference/index.md)
