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
