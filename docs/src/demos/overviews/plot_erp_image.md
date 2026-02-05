This demo shows how to create ERP image plots for visualizing single-trial EEG activity.

### What is an ERP Image?

An ERP image displays single-trial data as a 2D heatmap:

- **Y-axis**: Individual trials/epochs
- **X-axis**: Time points
- **Color**: Amplitude at each time-trial coordinate

This reveals trial-to-trial variability that is hidden when viewing only averaged ERPs.

### Why Use ERP Images?

**See variability beyond the average**:

Traditional ERPs show the average across trials, but ERP images show:
- Trial-to-trial amplitude fluctuations
- Temporal jitter in component latency
- Phase consistency vs. phase resetting

**Quality control**:

Quickly spot:
- Artifact-contaminated trials (bright/dark streaks)
- Drift across the experiment
- Individual outlier trials

**Understand components**:

Distinguish between:
- **Time-locked** activity: Consistent latency across trials (vertical bands)
- **Phase-locked** activity: Consistent phase but variable amplitude
- **Induced** activity: Not phase-locked to stimulus

### Layout Options

The demo shows three layout modes:

| Layout | Description |
|--------|-------------|
| **:single** | All selected channels in one column |
| **:grid** | Channels arranged in a grid |
| **:topo** | Channels positioned by scalp location |

### Visualization Features

**Boxcar averaging**:

Smooth the image by averaging across neighboring trials:
```julia
plot_erp_image(epochs, boxcar_average = 20)
```

Reduces noise while preserving overall patterns.

**Colorbar control**:

```julia
plot_erp_image(epochs, colorbar_plot = false)  # Hide colorbar
plot_erp_image(epochs, colorrange = (-50, 50))  # Custom range
```

**Optional ERP overlay**:

```julia
plot_erp_image(epochs, plot_erp = false)  # Hide the averaged ERP
```

By default, the averaged ERP is plotted above the image for reference.

### Interpretation

**Vertical bands**:

Strong vertical alignment indicates time-locked activity with consistent latency across trials.

**Gradual color changes**:

Smooth transitions suggest sustained activity or slow baseline drift.

**Random speckles**:

High-frequency noise or lack of consistent time-locked activity.

**Outlier rows**:

Trials with extreme values may indicate artifacts (blinks, movements, poor contact).

### Workflow Summary

This demo demonstrates:

1. **Load and preprocess** data (rereferencing, filtering)
2. **Extract epochs** (-2 to 4 seconds)
3. **Plot ERP images** with different layouts
4. **Apply boxcar averaging** to smooth trial-to-trial noise
5. **Customize visualization** with colorbars and color ranges

ERP images are a powerful complement to traditional ERP plots, revealing the underlying trial structure that averages can obscure.
