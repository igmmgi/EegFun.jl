# Plot Topography

This demo demonstrates creating topographic scalp maps to visualize the spatial distribution of EEG activity at specific time points.

This demo demonstrates creating topographic scalp maps to visualize the spatial distribution of EEG activity at specific time points.

### What are Topographic Maps?

Topographic maps (topoplots) show the spatial distribution of electrical activity across the scalp:

- **Color-coded amplitudes**: Activity levels represented by color
- **Interpolated surfaces**: Smooth maps between discrete electrode positions
- **Time-specific snapshots**: Activity at particular latencies or time windows

### Interpolation Methods

**`:thin_plate`**: Thin-plate spline (default; smooth, natural-looking)

**`:nearest`**: Nearest neighbor (no interpolation, fastest)

**`:shepard`**: Inverse distance weighting

**`:multiquadratic`**: Radial basis function

**`:spherical_spline`**: Spherical spline (accounts for head curvature)

> **Note**: All interpolation methods except `:spherical_spline` are implemented using the [ScatteredInterpolation.jl](https://github.com/eljungsk/ScatteredInterpolation.jl) package. The spherical spline method uses a custom implementation based on Perrin et al. (1989).

### Interpretation

**Focal activity**:

- Localized color patches suggest specific neural sources
- Sharp gradients indicate nearby sources

**Widespread activity**:

- Diffuse patterns suggest distributed processing
- Gradual transitions indicate distant or multiple sources

**Polarity conventions**:

- Warm colors (red/yellow): Positive voltage
- Cool colors (blue): Negative voltage
- Convention may vary by field

### Use Cases

**Visualize ERP components**:

- Show spatial distribution of P1, N170, P300, etc.
- Identify component topographies

**Compare conditions**:

- Side-by-side condition comparisons
- Difference topographies (condition A - B)

**Publication figures**:

- High-quality scalp maps
- Customizable appearance
- Multiple time points or conditions

### Customization Options

**Interpolation**:


- Method selection
- Grid resolution (gridscale)

**Appearance**:

- Colormap selection
- Color limits (ylim)
- Head outline radius

**Labels and markers**:

- Electrode positions
- Channel labels
- Font sizes and colors

**Colorbar**:

- Position and orientation
- Size and tick labels
- Show/hide per plot

### Working with Different Data Types

**Continuous data**: Average over time window

**Epoched data**: Specify epoch number and time window

**ERP data**: Average directly (already averaged)

## Workflow Summary

This demo shows topographic visualization workflows:

### 1. Basic Continuous Data

- Load and preprocess data
- Create topographic maps with different methods
- Customize time windows and appearance

### 2. Epoched Data

- Extract epochs from continuous data
- Plot topographies for specific epochs
- Customize interpolation and display

### 3. ERP Data

- Average epochs into ERPs
- Create condition-specific topographies
- Control colorbar placement

### 4. Multi-Panel Figures

- Create custom figure layouts
- Combine multiple topographies
- Control colorbar positions for each plot


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

# visually selected blink like artifact
EegFun.plot_topography(
    dat,
    sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02,
    gridscale = 75,
    ylim = (-200, 200),
    head_radius = 1.0,
)

# different methods
EegFun.plot_topography(
    dat,
    sample_selection = x -> x.time .>= 5.973 .&& x.time .<= 6.02,
    gridscale = 75,
    ylim = (-200, 200),
    head_radius = 1.0,
    method = :multiquadratic, # :nearest, :shepard, :spherical_spline, :thin_plate, :multiquadratic
)

# Various combinations
EegFun.plot_topography(dat, colorbar_plot = false, head_radius = 1.25)
EegFun.plot_topography(dat, gridscale = 250)
EegFun.plot_topography(dat, colormap = :inferno)
EegFun.plot_topography(dat, title = "Custom Title", title_fontsize = 30)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)
EegFun.plot_topography(dat, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, method = :spherical_spline)
EegFun.plot_topography(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :Cz]))
EegFun.plot_topography(dat, colorbar_size = 20, colorbar_position = (2, 1), colorbar_vertical = false)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Separate plots
EegFun.plot_topography(epochs[1], 1) # epoch 1
EegFun.plot_topography(epochs[2], 1) # epoch 2
EegFun.plot_topography(epochs) # TODO: aspect ration?; global scale?
EegFun.plot_topography(epochs, ylim = (-0.1, 0.1)) # TODO: aspect ration?; global scale?
EegFun.plot_topography(epochs, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6)



#################################
# ERP like data
#################################
erps = EegFun.average_epochs(epochs)

EegFun.plot_topography(erps, sample_selection = x -> x.time .>= 0.4 .&& x.time .<= 0.6, ylim = (-2, 2))
EegFun.plot_topography(erps, interval_selection = (0.4, 0.6), ylim = (-2, 2))

EegFun.plot_topography(erps, interval_selection = (0.4, 0.6), ylim = (-2, 2), colorbar_plot_numbers = [2], dims = (1, 2))

EegFun.plot_topography(
    erps,
    interval_selection = (0.4, 0.6),
    ylim = (-2, 2),
    colorbar_plot = true,
    colorbar_position = (2, 1),
    colorbar_vertical = false,
)

EegFun.plot_topography(erps)
EegFun.plot_topography(erps, ylim = (-2, 2))
EegFun.plot_topography(erps[2])
EegFun.plot_topography(erps[1], gridscale = 50)
EegFun.plot_topography(erps[2], gridscale = 1000)
EegFun.plot_topography(erps[1], colormap = :inferno)
EegFun.plot_topography(erps[2], title = "Custom Title", title_fontsize = 30)

```

:::

## See Also

- [API Reference](../reference/index.md)
