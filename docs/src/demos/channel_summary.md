# Channel Summary

This demo demonstrates how to generate summary statistics across all channels for data quality assessment.

This demo demonstrates how to generate summary statistics across all channels for data quality assessment.

### What is Channel Summary?

Channel summary provides aggregate statistics across all EEG channels, offering a quick overview of data characteristics. This complements channel-specific metrics by revealing patterns across the entire montage.

**Statistical measures**:

- Mean amplitude per channel
- Variance/standard deviation
- Minimum and maximum values
- Sample counts

### Use Cases

**Data quality overview**:

- Quick assessment of recording quality across all channels
- Identify systematic patterns or issues across the montage
- Compare data quality before and after preprocessing

**Channel comparison**:

- Identify channels with systematically different properties
- Detect asymmetries or systematic biases in the recording
- Verify that preprocessing affected channels as expected

## Workflow Summary

This demo shows channel summary analysis:

### 1. Generate Initial Summary

- Load and preprocess data (average reference, high-pass filter)
- Calculate summary statistics across all channels
- Display results using formatted table

### 2. Summary with Sample Selection

- Mark extreme values for exclusion
- Recalculate summary excluding artifacts
- Compare clean vs. raw summary statistics


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# summary statistics across all channels
summary = EegFun.channel_summary(dat)
EegFun.log_pretty_table(summary; title = "Initial Channel Summary")

# summary statistics across all channels excluding v. extreme values
EegFun.is_extreme_value!(dat, 200);
summary = EegFun.channel_summary(dat, sample_selection = EegFun.samples_not(:is_extreme_value_200))
EegFun.log_pretty_table(summary; title = "Channel Summary (excluding extreme values)")

# summary statistics across all Midline channels via predicate selection
summary = EegFun.channel_summary(dat, channel_selection = EegFun.channels(x -> endswith.(string.(x), "z")))
EegFun.log_pretty_table(summary; title = "Channel Summary (Midline)")
```

:::

## See Also

- [API Reference](../reference/index.md)
