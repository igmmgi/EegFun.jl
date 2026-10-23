# Plot ERP Measurement GUI

This demo shows how to use the interactive ERP measurement GUI for visual exploration and parameter tuning.

## What Does It Do?

`plot_erp_measurement_gui` opens an interactive window that lets you explore ERP measurements in real time. You can adjust the analysis interval, baseline window, measurement type, and channel — all with live visual feedback. This is useful for:

- **Teaching** — demonstrating how measurement choices affect results
- **Parameter exploration** — finding the right analysis window before batch processing
- **Visual validation** — confirming that automated measurements are sensible

## Relationship to `erp_measurements`

Use the GUI to explore and decide on your measurement parameters, then apply those parameters programmatically with `EegFun.erp_measurements()` for your final analysis.

## Workflow Summary

This demo covers:

## Basic Usage

- Launch the quick demo with default data and settings (dynamically 30 Hz low-pass filtered, `:Cz`, 300–500 ms):
  ```julia
  EegFun.plot_erp_measurement_gui()
  ```

## Custom Data & Settings

- Launch the GUI with specific initial channel, measurement type, and intervals:
  ```julia
  dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
  EegFun.plot_erp_measurement_gui(
      dat,
      channel = :Cz,
      analysis_type = "mean_amplitude",
      analysis_interval = (0.3, 0.5),
      baseline_interval = (-0.2, 0.0),
  )
  ```

## Code Examples

::: details Show Code

```julia
# Demo: Interactive ERP Measurement GUI
# Explore ERP measurements interactively with live visual feedback.
# Useful for teaching, parameter exploration, and visual validation
# before batch processing with erp_measurements().

using EegFun
using GLMakie

# Launch with default demo data (dynamically 30 Hz low-pass filtered, :Cz, 300-500 ms):
gui = EegFun.plot_erp_measurement_gui()

# Or pass custom data or specific settings:
# dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
# gui = EegFun.plot_erp_measurement_gui(
#     dat,
#     channel = :Cz,
#     analysis_type = "mean_amplitude",
#     analysis_interval = (0.3, 0.5),
#     baseline_interval = (-0.2, 0.0),
# )
```

:::

## See Also

- [API Reference](../../reference/index.md)
