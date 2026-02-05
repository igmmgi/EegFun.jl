# Plot Correlation Heatmap

This demo demonstrates how to compute and visualize channel correlation matrices for quality control and artifact detection.

This demo demonstrates how to compute and visualize channel correlation matrices for quality control and artifact detection.

### What is a Correlation Matrix?

A correlation matrix shows the Pearson correlation coefficient between all channel pairs:

- **Values range from -1 to +1**
- **High positive correlation** (near +1): Channels move together
- **Low correlation** (near 0): Channels are independent
- **Negative correlation** (near -1): Channels move in opposite directions

### Why Use Correlation Matrices?

**Quality Control**:

- **Identify bad channels**: Low correlation with all neighbors
- **Detect bridging**: Unexpectedly high correlation between channels
- **Assess reference choice**: Impact on spatial correlation structure

### Visualization Types

**Full Correlation Heatmap**:

Shows all channel-to-channel correlations as a matrix where:

- **Diagonal**: Always 1.0 (channel correlated with itself)
- **Off-diagonal**: Correlation between channel pairs
- **Color scale**: Red (high correlation) to blue (low correlation)

**Dual Selection Heatmap**:

Correlates two different channel sets:

- **Rows**: One set of channels (e.g., all EEG channels)
- **Columns**: Another set (e.g., EOG channels)

**Layout Overlay**:

Displays channel-to-channel correlations directly on the 2D sensor layout:

- Select a reference channel
- See its correlation with all other channels as colored dots
- **Use case**: Quick visual inspection of spatial correlation patterns

### Expected Patterns

**Healthy Data**:

- **High correlation with neighbors** (~0.7-0.9): Nearby channels capture similar signals
- **Decreasing correlation with distance**: Far channels less correlated
- **Symmetric matrix**: Correlation(A,B) = Correlation(B,A)

**Problematic Patterns**:

**Single bad channel**:

- Low correlation (<0.3) with all neighbors
- Stands out as a blue row/column in heatmap
- **Solution**: Interpolate or exclude

**Bridged channels**:

- Abnormally v. high correlation (>0.99) between adjacent channels
- Indicates electrical connection between electrodes
- **Solution**: Usually, more careful setup is required to avoid this

**EOG contamination**:

- High correlation between frontal channels and EOG
- Expected for Fp1/Fp2 (near eyes)
- **Solution**: Use ICA for removal if excessive

### Common Use Cases

**1. Post-Acquisition Quality Check**:

```julia
cm = correlation_matrix(dat)
plot_correlation_heatmap(cm)
```

Quick overview to identify problematic channels

**2. Targeted Channel Inspection**:

```julia
cm = correlation_matrix(dat, channel_selection = channels([:Fp1, :Fp2, :C3, :C4]))
```

Focus on specific regions or suspected bad channels

**3. EOG Artifact Assessment**:

```julia
cm = correlation_matrix_dual_selection(
    dat,
    channel_selection1 = channels(),  # All EEG
    channel_selection2 = channels([:vEOG, :hEOG])  # EOG
)
```

Identify which channels are affected by eye movements


### Workflow Summary

This demo shows:

1. **Compute full correlation matrix** for all channels
2. **Visualize as heatmap** for overview
3. **Compute subset** targeting specific channels
4. **Create EOG channels** for artifact assessment
5. **Dual selection correlation** to identify EOG-contaminated channels
6. **Layout overlay** for spatial visualization

The correlation matrix is a powerful diagnostic tool for both quality control and understanding spatial relationships in your data!


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
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Correlation Matrix
cm = EegFun.correlation_matrix(dat)
EegFun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")

cm = EegFun.correlation_matrix(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :C3, :C4]))
EegFun.plot_correlation_heatmap(cm, title = "Selected Correlations")

# Calculate EOG signals
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

# Calculate correlations between all channels and EOG channels
cm_eog = EegFun.correlation_matrix_dual_selection(
    dat,
    sample_selection = EegFun.samples(),  # All samples
    channel_selection1 = EegFun.channels(),  # All EEG channels
    channel_selection2 = EegFun.channels([:vEOG, :hEOG]),  # EOG channels
)

# Visualize which EEG channels are correlated with EOG
EegFun.plot_correlation_heatmap(cm_eog, title = "EEG-EOG Correlations", xlabel = "EOG Channels", ylabel = "EEG Channels")

# Display channel correlations in layout
cm = EegFun.correlation_matrix(dat)
EegFun.plot_layout_2d(layout_file, correlation_matrix = cm)
```

:::

## See Also

- [API Reference](../reference/index.md)
