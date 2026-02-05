# ERP Measurements

## Overview

## Overview

This demo shows how to extract quantitative measurements from ERP waveforms.

### ERP Measurements

Extract meaningful features from averaged waveforms:

**Amplitude Measures:**
- **Mean amplitude**: Average voltage in a time window
- **Peak amplitude**: Maximum/minimum voltage
- **Peak-to-peak**: Voltage difference between peaks

**Latency Measures:**
- **Peak latency**: Time of maximum deflection
- **Onset latency**: When component begins
- **Fractional area latency**: Time point dividing area

**Area Measures:**
- **Integral**: Total area under the curve
- **Rectified area**: Absolute area (ignores polarity)

### Use Cases

- Quantify component amplitudes (N1, P3, N400, etc.)
- Compare experimental conditions statistically
- Report ERP characteristics in publications


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
