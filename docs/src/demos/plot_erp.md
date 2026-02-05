# Plot ERP

## Overview

## Overview

This demo demonstrates plotting event-related potentials with various customization options.

### ERP Plotting

Visualize averaged waveforms across conditions:
- **Multiple conditions**: Overlay for comparison
- **Channel selection**: Plot specific electrodes or regions
- **Customization**: Colors, line styles, confidence intervals

### Plot Features

- Time on x-axis, amplitude on y-axis
- Condition overlays with labels
- Shaded confidence intervals (optional)
- Time markers for stimulus onset, analysis windows
- Topographic insets at specific time points

### Applications

- Publication-quality figures
- Compare experimental conditions
- Identify and measure ERP components


## Code Examples

::: details Show Code

```julia
using EegFun
using JLD2

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# EPOCHS -> ERPs
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))
erps = EegFun.average_epochs(epochs)

out = EegFun.plot_erp(erps, layout = :single)
EegFun.plot_erp(erps, layout = :grid)
EegFun.plot_erp(erps[1], layout = :topo)


EegFun.plot_erp(erps, layout = :grid, legend_channel = [:Fp1, :M2], yreversed = true)

EegFun.plot_erp(
    erps,
    channel_selection = EegFun.channels([:F3, :Cz, :PO7, :PO8, :Fp1, :Fp2]),
    layout = :grid,
    layout_grid_dims = (3, 2),
    layout_grid_skip_positions = [(2, 1)],
)

EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8, :Fp1, :Fp2, :F3]), layout = :grid, layout_grid_dims = (2, 3))

EegFun.plot_erp(
    erps,
    channel_selection = EegFun.channels([:Cz, :PO7, :PO8, :Fp1, :Fp2, :F3, :T8, :F4]),
    layout = :grid,
    layout_grid_dims = (3, 4),
    layout_grid_skip_positions = [(2, 1), (2, 3)],
)

EegFun.plot_erp(
    erps,
    channel_selection = EegFun.channels([:Cz, :PO7, :PO8, :Fp1, :Fp2, :F3]),
    layout = :grid,
    layout_grid_dims = (2, 4),
    layout_grid_skip_positions = [(2, 1), (2, 3)],
    layout_grid_rowgap = 0,
    layout_grid_colgap = 0,
    figure_padding = (150, 150, 150, 150),
)

# ERP Plots (:single)
EegFun.plot_erp(erps, average_channels = false, colormap = :viridis, legend_nbanks = 12)
EegFun.plot_erp(erps[1], average_channels = false, colormap = :viridis, legend_nbanks = 12)

EegFun.plot_erp(erps, average_channels = true)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8]), average_channels = true)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8]), average_channels = false)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8]), average_channels = false, legend_nbanks = 3)
EegFun.plot_erp([erps[1], erps[1]], channel_selection = EegFun.channels([:PO8]))

# ERP Plots (:grid)
EegFun.plot_erp(erps, layout = :grid)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz, :PO7, :PO8, :Fp1, :Fp2, :F3]), layout = :grid)

# EPR Plots (:topo)
EegFun.plot_erp(erps, layout = :topo)
EegFun.plot_erp([erps[1], erps[1]], layout = :topo)
EegFun.plot_erp(erps, layout = :topo, channel_selection = EegFun.channels([:Fp1, :Fp2, :PO8]))

# Combined plots
using GLMakie
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 1], width = Relative(0.2), height = Relative(0.2), halign = 0, valign = 0)
EegFun.plot_erp!(fig, ax1, erps, average_channels = true)
EegFun.plot_topography!(
    fig,
    ax2,
    erps[1];
    point_plot = false,
    label_plot = false,
    colorbar_plot = true,
    colorbar_width = Relative(0.03),
    colorbar_height = Relative(0.2),
    colorbar_tellheight = false,
    colorbar_tellwidth = false,
    colorbar_position = (1, 1),
    colorbar_halign = 0.25,
    colorbar_valign = 0,
    colorbar_flipaxis = true,
)
fig

GLMakie.closeall()
```

:::

## See Also

- [API Reference](../reference/index.md)
