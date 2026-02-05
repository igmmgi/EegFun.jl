# Plot ERP Image

## Overview

## Overview

This demo shows how to create ERP image plots for visualizing trial-by-trial variations.

### ERP Images

ERP images display single-trial data as 2D color-coded images:
- **Rows**: Individual trials
- **Columns**: Time points
- **Color**: Amplitude
- **Sorting**: Trials can be sorted by RT, amplitude, or other variables

### Why Use ERP Images?

- **Visualize variability**: See trial-to-trial differences beyond the average
- **Identify artifacts**: Spot trials with transient noise
- **Phase consistency**: Observe whether components are time-locked vs. phase-locked
- **Sorting effects**: Reveal dynamics related to behavior or stimulus properties

### Applications

- Supplement traditional ERPs with single-trial information
- Identify outlier trials
- Explore relationships between neural activity and behavior


## Code Examples

::: details Show Code

```julia
using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

EegFun.plot_erp_image(epochs[1], layout = :single)
EegFun.plot_erp_image(epochs[1], layout = :single, channel_selection = EegFun.channels([:Fp1, :Fp2]))

EegFun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = false)
EegFun.plot_erp_image(epochs[1], layout = :grid, colorbar_plot = true)

EegFun.plot_erp_image(epochs[1], layout = :topo)
EegFun.plot_erp_image(epochs[1], channel_selection = EegFun.channels([:PO7]))
EegFun.plot_erp_image(epochs[1], channel_selection = EegFun.channels([:Fp1]), plot_erp = false)


EegFun.plot_erp_image(epochs[1], layout = :single)

(; fig, axes) = EegFun.plot_erp_image(
    epochs[1],
    # channel_selection = EegFun.channels([:Fp1, :Fp2]),
    # channel_selection = EegFun.channels([:Fp1, :Fp2]),
    layout = :topo,
    boxcar_average = 20,
    colorrange = (-50, 50),
)


# TODO: no electrode labels in title
EegFun.plot_erp_image(epochs[1], layout = :topo)
EegFun.plot_erp_image(epochs[1], layout = :grid)
```

:::

## See Also

- [API Reference](../reference/index.md)
