# Plot Filter

## Overview

## Overview

This demo demonstrates visualization of filter frequency and phase responses.

### Filter Response Plots

Inspect filter characteristics before applying:
- **Magnitude response**: How much each frequency is attenuated
- **Phase response**: Time delays introduced at different frequencies
- **Impulse response**: Filter behavior in time domain

### Why Visualize Filters?

- Verify cutoff frequencies are correct
- Check transition band steepness
- Identify potential artifacts (ringing, phase distortion)
- Document filter parameters for methods sections

### Filter Types

- **High-pass**: Remove slow drifts
- **Low-pass**: Remove high-frequency noise
- **Band-pass**: Isolate specific frequencies
- **Notch**: Remove line noise (50/60 Hz)


## Code Examples

::: details Show Code

```julia
using EegFun

# Create lowpass IIR filter using create_filter
filter_info = EegFun.create_lowpass_filter(30.0, 256.0; filter_method = "iir")

# Plot filter response
EegFun.plot_filter_response(filter_info)

# Test with custom parameters
EegFun.plot_filter_response(
    filter_info,
    title = "Custom Lowpass Filter Plot",
    actual_color = :blue,
    actual_linewidth = 3,
    reference_lines = [-3, -12, -24],
    reference_color = :red,
    n_points = 1000,
)

# Create highpass IIR filter using create_filter
filter_info = EegFun.create_highpass_filter(1.0, 256.0; filter_method = "iir")
EegFun.plot_filter_response(filter_info, title = "High-pass Filter", actual_color = :green)

# Create FIR filter using create_filter
filter_info = EegFun.create_lowpass_filter(40.0, 256.0; filter_method = "fir")
EegFun.plot_filter_response(filter_info, title = "FIR Lowpass Filter", actual_color = :purple)

# Test additional filter with separate plotting
filter_info = EegFun.create_highpass_filter(0.5, 256.0; filter_method = "iir")
EegFun.plot_filter_response(filter_info, title = "High-pass Filter with Plot")





```

:::

## See Also

- [API Reference](../reference/index.md)
