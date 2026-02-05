# Plot Joint Probability

## Overview

## Overview

This demo shows how to visualize channel joint probability distributions for quality assessment.

### Joint Probability Analysis

Joint probability detects outlier channels based on multi-dimensional distributions:
- Compares each channel to all others simultaneously
- More sensitive than univariate metrics
- Identifies channels with unusual data characteristics

### Visualization

Plot shows how each channel relates to the overall distribution:
- **Low probability**: Potential bad channels
- **Consistent probabilities**: Clean data
- **Spatial patterns**: Systematic issues with recording

### Use Cases

- Automated bad channel detection
- Complement other quality metrics
- Pre-processing quality control


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# test plot_channel_summary
jp = EegFun.channel_joint_probability(dat)
EegFun.plot_joint_probability(jp)
```

:::

## See Also

- [API Reference](../reference/index.md)
