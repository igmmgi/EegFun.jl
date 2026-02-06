# Plot Epochs

This demo demonstrates visualizing individual epoch waveforms with flexible layout options.

This demo demonstrates visualizing individual epoch waveforms with flexible layout options.

### What is Epoch Plotting?

Epoch plots display individual trial timecourses, showing:

- **Single-trial waveforms**: Raw epoch data (not averaged)
- **Time-locked activity**: Aligned to stimulus or response events
- **Trial variability**: Individual differences across repetitions
- **ERP overlay**: Optional averaged waveform on top

### Layout Options

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

### Use Cases

**Trial-level inspection**:

- View variability across trials
- Identify outlier epochs
- Assess signal quality

**Verify epoching**:

- Check time window selection
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

### Customization Options

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

### Interactive Features

When `interactive=true` (default):

- **Control panel** (press 'c'): Toggle conditions, adjust baseline
- **Keyboard shortcuts**: Scaling, help
- **Linked axes**: Synchronized zoom/pan across subplots

## Workflow Summary

This demo shows epoch visualization workflows:

### 1. Load and Prepare Data

- Read raw data
- Apply preprocessing
- Extract epochs

### 2. Basic Plotting

- Plot single condition
- Select specific channels
- Compare conditions with overlay

### 3. Layout Options

- Single plot: All channels together
- Grid layout: Individual subplots
- Topo layout: Spatial arrangement

### 4. Channel Selection

- Focus on specific channels
- Combine with different layouts
- Maintain consistent scaling


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# Basic plots for Epochs
EegFun.plot_epochs(epochs[1]) # not that useful as crowded!
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:PO7, :PO8]))

EegFun.plot_epochs(epochs[1], layout = :grid)
EegFun.plot_epochs(epochs[1], layout = :grid, channel_selection = EegFun.channels([:Fp1, :Fp2, :PO7, :PO8]))

EegFun.plot_epochs(epochs[1], layout = :topo)
EegFun.plot_epochs(epochs[1], layout = :topo, add_xy_origin = false, theme_fontsize = 10, layout_topo_show_scale = true)

EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]); layout = :topo)

```

:::

## See Also

- [API Reference](../reference/index.md)
