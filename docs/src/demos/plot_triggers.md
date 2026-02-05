# Plot Triggers

## Overview

## Overview

This demo demonstrates visualization of event markers and triggers in continuous data.

### Trigger Visualization

Display when experimental events occurred:
- **Trigger markers**: Vertical lines at event times
- **Trigger codes**: Different colors/styles for different event types
- **Timing validation**: Verify event sequences and ISIs

### Use Cases

- **Quality control**: Verify triggers were recorded correctly
- **Timing analysis**: Check inter-stimulus intervals
- **Event sequence**: Confirm experimental protocol
- **Troubleshooting**: Identify missing or spurious triggers

### Features

- Overlay on continuous data
- Color-coded by trigger type
- Zoom to inspect timing precision
- Summary statistics of trigger counts


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# count from raw file
count = EegFun.trigger_count(dat)

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

count = EegFun.trigger_count(dat)

EegFun.plot_trigger_overview(dat)
EegFun.plot_trigger_overview(dat; ignore_triggers = [3])

EegFun.plot_trigger_timing(dat)
EegFun.plot_trigger_timing(dat; ignore_triggers = [3])
```

:::

## See Also

- [API Reference](../reference/index.md)
