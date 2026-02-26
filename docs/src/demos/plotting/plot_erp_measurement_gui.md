# Plot ERP Measurement GUI

This demo shows how to use the interactive ERP measurement GUI for visual exploration and parameter tuning.

### What Does It Do?

`plot_erp_measurement_gui` opens an interactive window that lets you explore ERP measurements in real time. You can adjust the analysis interval, baseline window, measurement type, and channel — all with live visual feedback. This is useful for:

- **Teaching** — demonstrating how measurement choices affect results
- **Parameter exploration** — finding the right analysis window before batch processing
- **Visual validation** — confirming that automated measurements are sensible

### Relationship to `erp_measurements`

Use the GUI to explore and decide on your measurement parameters, then apply those parameters programmatically with `EegFun.erp_measurements()` for your final analysis.

## Workflow Summary

This demo covers:

### Basic Usage

- Launch the GUI with all conditions overlaid

### Single Condition

- Explore one condition at a time

### Pre-configured Settings

- Open the GUI with specific initial channel, measurement type, and intervals


## Code Examples

::: details Show Code

```julia
# Demo: Interactive ERP Measurement GUI
# Explore ERP measurements interactively with live visual feedback.
# Useful for teaching, parameter exploration, and visual validation
# before batch processing with erp_measurements().

using EegFun

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")

# All conditions overlaid
EegFun.plot_erp_measurement_gui(dat)

# Single condition
EegFun.plot_erp_measurement_gui(dat[1])

# Select some specific initial settings
EegFun.plot_erp_measurement_gui(
    dat[1],
    channel = :Cz,
    analysis_type = "mean_amplitude",
    analysis_interval = (0.3, 0.5),
    baseline_interval = (-0.2, 0.0),
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
