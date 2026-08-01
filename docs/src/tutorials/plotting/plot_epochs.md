# Plot Epochs

This demo demonstrates visualizing individual epoch waveforms with flexible layout options.

## What is Epoch Plotting?

Epoch plots display individual trial timecourses, showing:

- **Single-trial waveforms**: Raw epoch data (not averaged)
- **Time-locked activity**: Aligned to stimulus or response events
- **Trial variability**: Individual differences across repetitions
- **ERP overlay**: Optional averaged waveform on top

## Layout Options

**`:single`** (default):

- All channels overlaid on one plot
- Best for comparing few channels
- Clear view of channel differences

**`:grid`**:

- Each channel in separate subplot
- Auto-calculated grid dimensions
- Best for many channels
- Synchronized axes

**`:topo`**:

- Channels arranged by scalp position
- Preserves spatial relationships
- Intuitive interpretation
- Optional scale axis

## Use Cases

**Trial-level inspection**:

- View variability across trials
- Identify outlier epochs
- Assess signal quality

**Verify epoching**:

- Check time interval selection
- Confirm event alignment
- Validate baseline correction

**Artifact detection**:

- Spot trial-specific artifacts
- Find contaminated epochs
- Guide rejection decisions

**Condition comparison**:

- Overlay multiple conditions
- Compare trial-level activity
- Assess consistency

## Customization Options

**Channel selection**:

- Select specific channels or regions
- Focus on electrodes of interest

**Visual styling**:

- Individual trial transparency
- Average line width
- Color schemes
- Axis limits

**Layout parameters**:

- Grid dimensions
- Spacing and gaps
- Scale axis (topo layout)

## Interactive Features

When `interactive=true` (default):

- **Control panel** (press 'c'): Toggle conditions, adjust baseline
- **Keyboard shortcuts**: Scaling, help
- **Linked axes**: Synchronized zoom/pan across subplots

## Workflow Summary

This demo shows epoch visualization workflows:

## Load and Prepare Data

- Read raw data
- Apply preprocessing
- Extract epochs

## Basic Plotting

- Plot single condition
- Select specific channels
- Compare conditions with overlay

## Layout Options

- Single plot: All channels together
- Grid layout: Individual subplots
- Topo layout: Spatial arrangement

## Channel Selection

- Focus on specific channels
- Combine with different layouts
- Maintain consistent scaling


## Code Examples

::: details Show Code

```julia
# Demo: Epoch Plotting
# Shows epoch visualization with channel selection and overlay options.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# Basic plots for Epochs
EegFun.plot_epochs(epochs) # not that useful as crowded!
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:PO7, :PO8]))

EegFun.plot_epochs(epochs, layout = :grid)
EegFun.plot_epochs(epochs, layout = :grid, channel_selection = EegFun.channels([:Fp1, :Fp2, :PO7, :PO8]))

EegFun.plot_epochs(epochs, layout = :topo)
EegFun.plot_epochs(epochs, layout = :topo, add_xy_origin = false, theme_fontsize = 10, layout_topo_show_scale = true)

EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)






```

:::

## See Also

- [API Reference](../../reference/index.md)
