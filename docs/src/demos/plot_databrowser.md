# Plot Databrowser

## Overview

## Overview

This demo demonstrates the interactive databrowser for exploring EEG data.

### Interactive Data Browser

The databrowser provides real-time visualization and exploration:
- **Continuous data**: Scroll through raw recordings
- **Epoch data**: Navigate individual trials
- **ICA data**: Inspect component time courses

### Key Features

**Mouse Interactions:**
- Scroll/pan through time
- Select regions for analysis
- Select channels for repair
- Zoom in/out

**Keyboard Shortcuts:**
- Arrow keys: Navigate
- 'i': Show help
- 'r': Open channel repair menu
- 'c': Clear selections

### Use Cases

- Initial data inspection
- Identify artifacts visually
- Select regions for spectral analysis
- Interactive channel repair
- Validate preprocessing steps


## Code Examples

::: details Show Code

```julia
using EegFun
using JLD2

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Basic databrowser
EegFun.plot_databrowser(dat);

# return some analysis settings
(; fig, ax, analysis_settings) = EegFun.plot_databrowser(dat)
dat_new = EegFun.apply_analysis_settings(dat, analysis_settings)

EegFun.plot_databrowser(dat_new)

# Add additional columns that can be viewed in the databrowser
EegFun.is_extreme_value!(dat, 100);
EegFun.mark_epoch_windows!(dat, [1, 2], [-0.2, 1.0]) # simple epoch marking with trigger 1 and 3

EegFun.plot_databrowser(dat)

# try some custom styling
EegFun.plot_databrowser(dat; :channel_line_width => 2, :selection_color => (:red, 0.5))

# with ica 
ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100), percentage_of_data = 25)

EegFun.plot_databrowser(dat, ica_result)

# EPOCHS
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

# databrowser for epochs
EegFun.plot_databrowser(epochs[1])  # can now browse epochs
EegFun.plot_databrowser(epochs[1], ica_result)

# plot_databrowser epochs from file
jldsave("dat.jld2", data = dat)
jldsave("epochs.jld2", data = epochs)
jldsave("ica.jld2", data = ica_result)

# Plot by loading from file
EegFun.plot_databrowser("dat.jld2")
EegFun.plot_databrowser("epochs.jld2")
EegFun.plot_databrowser("dat.jld2", "ica.jld2")
EegFun.plot_databrowser("epochs.jld2", "ica.jld2")
```

:::

## See Also

- [API Reference](../reference/index.md)
