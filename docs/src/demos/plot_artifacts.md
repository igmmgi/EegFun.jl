# Plot Artifacts

## Overview

## Overview

This demo shows how to visualize detected artifacts in epoched data.

### Artifact Visualization

Visualize which channels and epochs contain artifacts:
- **Channel-wise**: See which electrodes are most affected
- **Epoch-wise**: Identify problematic trials
- **Temporal patterns**: When artifacts occur within epochs

### Plot Components

- Individual epoch waveforms with artifact markers
- Bad channel/epoch summaries
- Rejection statistics across conditions


## Code Examples

::: details Show Code

```julia
using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs[1], artifacts[1])
```

:::

## See Also

- [API Reference](../reference/index.md)
