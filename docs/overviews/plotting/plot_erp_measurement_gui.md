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

### 1. Basic Usage

- Launch the GUI with all conditions overlaid

### 2. Single Condition

- Explore one condition at a time

### 3. Pre-configured Settings

- Open the GUI with specific initial channel, measurement type, and intervals
