# Channel Repair

This demo demonstrates channel interpolation methods for repairing bad electrodes using spatial interpolation techniques.

This demo demonstrates channel interpolation methods for repairing bad electrodes using spatial interpolation techniques.

### What is Channel Repair?

Channel repair (interpolation) estimates the signal at bad electrodes using data from neighboring good channels. This preserves data quality while maintaining electrode count for spatial analyses.

### When to Repair Channels

Channel repair is appropriate when:

- Individual electrodes have poor contact/noisy signal
- Isolated channels show excessive noise or artifacts
- A small number of channels are affected while most data is clean
- You need to maintain electrode count for spatial analyses (e.g., source localization)

**Do not repair** when:

- Too many channels are bad
- The entire dataset is noisy
- Bad channels cluster together spatially

### Interpolation Methods

**Neighbor Interpolation**:

- Weighted average of spatially nearby electrodes
- Fast and computationally efficient
- Good for isolated bad channels
- Requires neighbor calculation based on distance threshold

**Spherical Spline**:

- Uses spherical spline functions to model scalp potential distribution

### Best Practices

**Timing**:

- Identify bad channels using quality metrics first
- Repair before re-referencing (reference calculation needs all channels)
- Repair before averaging or statistical analysis

**Limits**:

- Avoid interpolating too many channels
- Consider rejecting datasets with an excessive number of bad channels

**Validation**:

- Visually verify repair quality (before/after comparison)
- Check that interpolated channels match neighbors

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

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eegfun_data(dat, layout_file)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

EegFun.plot_databrowser(dat)

# Select a channel to repair and make this channel noisy!
channel_to_repair = :Cz
dat.data[!, channel_to_repair] .+= randn(size(dat.data[:, channel_to_repair])) * 200 # v. noisy!

# We can now see this noisy channel in the databrowser
# NB. we can actually press "R" and select Cz and apply the repair in the browser
EegFun.plot_databrowser(dat)

# Try neighbor interpolation
EegFun.repair_channels!(dat, [channel_to_repair], method = :neighbor_interpolation)

# Cz is now repaired
EegFun.plot_databrowser(dat)

# Try neighbor interpolation
EegFun.repair_channels!(dat, [channel_to_repair], method = :spherical_spline)

# Cz is now repaired
EegFun.plot_databrowser(dat)


```

:::

## See Also

- [API Reference](../reference/index.md)
