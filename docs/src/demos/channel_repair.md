# Channel Repair

This demo demonstrates channel interpolation methods for repairing bad electrodes using spatial interpolation techniques.

This demo demonstrates channel interpolation methods for repairing bad electrodes using spatial interpolation techniques.

### What is Channel Repair?

Channel repair (interpolation) estimates the signal at bad electrodes using data from neighboring good channels. This preserves data quality while maintaining electrode count for spatial analyses.

### When to Repair Channels

Channel repair is appropriate when:

- Individual electrodes have poor contact or high impedance
- Isolated channels show excessive noise or artifacts
- A small number of channels are affected while most data is clean
- You need to maintain electrode count for spatial analyses (e.g., source localization)

**Do not repair** when:

- Too many channels are bad (>10-15% of total)
- The entire dataset is noisy
- Bad channels cluster together spatially
- Channel shows intermittent problems (repair won't fix instability)

### Interpolation Methods

**Neighbor Interpolation**:

- Weighted average of spatially nearby electrodes
- Fast and computationally efficient
- Good for isolated bad channels
- Requires neighbor calculation based on distance threshold

**Spherical Spline**:

- Uses spherical spline functions to model scalp potential distribution
- More accurate spatial interpolation
- Preserves topographic patterns better
- Recommended for publication-quality analyses

### Best Practices

**Timing**:

- Identify bad channels using quality metrics first
- Repair before re-referencing (reference calculation needs all channels)
- Repair before averaging or statistical analysis

**Limits**:

- Avoid interpolating too many channels (>10% is problematic)
- Document exact number and labels of repaired channels in methods
- Consider rejecting datasets with excessive bad channels

**Validation**:

- Visually verify repair quality (before/after comparison)
- Check that interpolated channels match neighbors
- Review repair plots to confirm successful interpolation

## Workflow Summary

This demo shows channel repair workflows:

### 1. Identify Channels to Repair

- Load and preprocess data
- Select specific channels for demonstration
- Verify channels exist in the dataset

### 2. Calculate Neighbor Relationships

- Compute spatial neighbors based on 3D electrode positions
- Use distance threshold to define neighborhood
- Required for neighbor interpolation method

### 3. Apply Interpolation

- Test neighbor interpolation method
- Test spherical spline method
- Store original data for comparison

### 4. Validate Repair Quality

- Compare interpolated vs. original data
- Verify that interpolation changed the data as expected
- Use visual inspection to confirm successful repair


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

# Select a few channels to repair
test_channels = [:Fp1]
available_channels = dat.layout.data.label
channels_to_repair = filter(ch -> ch in available_channels, test_channels)

# Calculate neighbors first
EegFun.get_neighbours_xyz!(dat.layout, 0.5)

# Store original data for comparison
original_data = copy(dat.data)

# Try neighbor interpolation
EegFun.repair_channels!(dat, channels_to_repair, method = :neighbor_interpolation)

# Try spherical spline
EegFun.repair_channels!(dat, channels_to_repair, method = :spherical_spline)

# Check if data changed (using isapprox to handle floating point precision)
data_changed = any(any(.!isapprox.(dat.data[!, ch], original_data[!, ch], rtol = 1e-10)) for ch in channels_to_repair)
```

:::

## See Also

- [API Reference](../reference/index.md)
