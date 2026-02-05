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
