# Plot ERP

This demo demonstrates plotting event-related potentials (ERPs) with comprehensive customization options.

This demo demonstrates plotting event-related potentials (ERPs) with comprehensive customization options.

### What is an ERP Plot?

ERP plots visualize averaged brain responses time-locked to events:

- **Time course**: Amplitude changes over time
- **Waveforms**: Characteristic positive and negative deflections
- **Condition comparison**: Overlay multiple experimental conditions
- **Channel-specific**: View individual electrodes or averages

### Layout Options

The demo shows three layout modes:

| Layout | Description |
|--------|-------------|
| **:single** | One plot with overlaid conditions |
| **:grid** | Multiple subplots (one per channel) |
| **:topo** | Channels arranged by scalp location |

### Single Layout

**Average across channels**:
```julia
plot_erp(erps, average_channels = true)
```
Shows grand average waveform across all selected channels.

**Individual channels**:
```julia
plot_erp(erps, average_channels = false, colormap = :viridis)
```
Overlays all channels with color-coding.

**Selected channels**:
```julia
plot_erp(erps, 
    channel_selection = channels([:Cz, :PO7, :PO8]),
    average_channels = false
)
```
Shows only specified channels.

### Grid Layout

Displays multiple channels as subplots:

```julia
plot_erp(erps, layout = :grid)
```

**Custom grid dimensions**:
```julia
plot_erp(erps,
    channel_selection = channels([:F3, :Cz, :PO7, :PO8, :Fp1, :Fp2]),
    layout = :grid,
    layout_grid_dims = (3, 2)  # 3 rows × 2 columns
)
```

**Skip positions**:
```julia
plot_erp(erps,
    layout_grid_dims = (3, 4),
    layout_grid_skip_positions = [(2, 1), (2, 3)]  # Leave empty
)
```
Creates custom layouts with empty spaces.

**Adjust spacing**:
```julia
plot_erp(erps,
    layout = :grid,
    layout_grid_rowgap = 0,  # No vertical gap
    layout_grid_colgap = 0   # No horizontal gap
)
```

### Topographic Layout

Arranges channels by scalp position:

```julia
plot_erp(erps, layout = :topo)
```

Each channel plotted at its actual spatial location for intuitive interpretation.

### Customization Options

**Y-axis orientation**:
```julia
plot_erp(erps, yreversed = true)  # Negative up (common convention)
```

**Legend placement**:
```julia
plot_erp(erps,
    legend_channel = [:Fp1, :M2],  # Channels for legend
    legend_nbanks = 3              # Number of legend columns
)
```

**Figure padding**:
```julia
plot_erp(erps,
    figure_padding = (150, 150, 150, 150)  # left, right, bottom, top
)
```

### Combining with Topography

Create publication-quality figures with embedded topographies:

```julia
using GLMakie
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 1], width = Relative(0.2), height = Relative(0.2))

plot_erp!(fig, ax1, erps, average_channels = true)
plot_topography!(fig, ax2, erps[1],
    point_plot = false,
    label_plot = false,
    colorbar_plot = true
)
```

Shows ERP waveform with scalp distribution at a specific time point.

### Common Use Cases

**Condition comparison**:
- Overlay multiple experimental conditions
- Identify differences in amplitude or latency
- Statistical windows highlighted

**Component identification**:
- Classic ERP components (N1, P1, N170, P3, etc.)
- Measure peak amplitudes and latencies
- Compare across channels

**Publication figures**:
- High-quality vector graphics
- Customizable colors and styles
- Grid layouts for multiple channels

### Interpretation

**Positive/negative deflections**:
- **P1, P2, P3**: Positive peaks (often plotted downward with yreversed = true)
- **N1, N2, N4**: Negative peaks (plotted upward)

**Typical components**:
- **P1/N1** (~100-200 ms): Early sensory processing
- **N170** (~170 ms): Face perception (occipito-temporal)
- **P3** (~300-600 ms): Attention, memory updating
- **N400** (~400 ms): Semantic processing

### Workflow Summary

This demo shows:

1. **Basic plotting**: All three layouts
2. **Channel averaging**: Grand average vs individual channels
3. **Custom grids**: Flexible subplot arrangements
4. **Grid customization**: Gaps, skip positions, padding
5. **Legend control**: Placement and formatting
6. **Combined plots**: ERP + topography insets

ERP plots are the foundation of event-related brain potential analysis!


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
