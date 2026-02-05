# ERP Measurements

This demo demonstrates extracting quantitative measurements from ERP waveforms for statistical analysis and reporting.

This demo demonstrates extracting quantitative measurements from ERP waveforms for statistical analysis and reporting.

### What are ERP Measurements?

ERP measurements quantify specific features of averaged event-related potentials:

- **Amplitude**: Voltage magnitude at specific times or peaks
- **Latency**: Timing of component peaks or onsets
- **Area/Integral**: Total activity over a time window
- **Peak-to-peak**: Voltage range between positive and negative deflections

### Measurement Types

**Mean Amplitude**:

- Average voltage in a time window
- Most robust amplitude measure
- Less sensitive to noise than peak measures
- Standard for most ERP components

**Peak Amplitude**:

- Maximum (or minimum) voltage in window
- Captures strongest response
- Can be affected by noise
- Useful for P300, N400, etc.

**Peak Latency**:

- Time of maximum deflection
- Identifies when component peaks
- Sensitive to individual differences
- Important for timing analyses

**Fractional Area Latency**:

- Time point dividing area (e.g., 50% of total area)
- More robust than peak latency
- Less affected by noise and waveform shape
- Better reflects central tendency

**Area/Integral**:

- Total voltage × time in window
- Captures sustained activity
- Less sensitive to brief noise
- Good for slow components

**Peak-to-Peak**:

- Voltage difference between positive and negative peaks
- Useful for biphasic components
- Captures full deflection range

### Interactive GUI

**`plot_erp_measurement_gui`**:

- Visual interface for exploring measurements
- Adjust time windows interactively
- Select measurement types
- Preview results before batch processing

### Batch Processing

**`erp_measurements`**:

- Process multiple files at once
- Apply consistent measurement parameters
- Export to CSV for statistical analysis
- Includes metadata (file, condition, channel)

### Best Practices

**Choose appropriate measures**:

- **Mean amplitude**: Default for most components
- **Peak measures**: When timing precision matters
- **Fractional latency**: For robust timing analysis
- **Area**: For sustained or variable waveforms

**Define time windows carefully**:

- Based on grand averages or literature
- Should capture component of interest
- Avoid overlapping components when possible

**Baseline correction**:

- Apply before measurements
- Use pre-stimulus interval
- Ensures consistent zero reference

**Multiple measurements**:

- Combine amplitude and latency
- Use area for validation
- Report multiple metrics when appropriate

### Typical Workflow

1. **Visualize ERPs** to identify components
2. **Use GUI** to explore measurement parameters
3. **Define time windows** based on grand average
4. **Batch process** all files with `erp_measurements`
5. **Export to CSV** for statistical analysis

## Workflow Summary

This demo shows ERP measurement extraction:

### 1. Interactive Exploration

- Launch `plot_erp_measurement_gui`
- Visualize ERPs
- Explore different measurement types
- Adjust time windows interactively

### 2. Batch Processing

- Define measurement parameters
- Process multiple files with `erp_measurements`
- Select conditions and channels
- Specify analysis and baseline intervals

### 3. Export Results

- Measurements saved to CSV
- Includes metadata (file, condition, channel)
- Ready for statistical analysis
- Reproducible parameters documented


## Code Examples

::: details Show Code

```julia
"""
Tutorial: ERP Measurement Options

This script provides an introduction to the ERP measurement capabilities 
in EegFun for extracting quantitative features from ERP data.

1. Amplitude measurements (mean, peak)
2. Latency measurements (peak, fractional)
3. Area/integral measurements
4. Peak-to-peak measurements
"""

using EegFun
dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")


EegFun.plot_erp_filter_gui(dat)

EegFun.plot_erp_measurement_gui(dat)
EegFun.plot_erp_measurement_gui(dat[1])

# ----------------------------------------------------------------------------
# Amplitude Measurements
# ----------------------------------------------------------------------------

# dat = EegFun.read_data("./data/files/erps/example1_erps_good.jld2")
# EegFun.plot_erp(dat, condition_selection = EegFun.conditions([1]), channel_selection = EegFun.channels([:Pz]), baseline_interval = (0, 0))

input_dir = "./resources/data/erps"
file_pattern = "erps_good"

# Mean amplitude in a time window
mean_amp = EegFun.erp_measurements(
    file_pattern,
    "max_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    channel_selection = EegFun.channels(),
    # channel_selection = EegFun.channels([:Pz, :Cz, :Fz]),
    analysis_interval = (0.6, 0.8),
    baseline_interval = (0.0, 0.0),
)


dat = EegFun.read_data("./resources/data/erps/example1_erps_good.jld2")
EegFun.plot_erp(dat, condition_selection = EegFun.conditions([1]), channel_selection = EegFun.channels([:Pz]), baseline_interval = (0, 0))


EegFun.plot_erp_measurement_gui(dat[1])
EegFun.plot_erp_measurement_gui(dat)```

:::

## See Also

- [API Reference](../reference/index.md)
