# Epoch Extraction

This demo explores the core functionality of segmenting continuous data into epochs (trials) based on event markers, and basic manipulation of these segments.

### What are Epochs?

Epoching is the process of extracting specific time windows around events (stimuli or responses) from a continuous recording. This allows for trial-based analysis and averaging to reveal Event-Related Potentials (ERPs).

**Key features of EegFun's epoching**:
- Precise time-relative segmentation
- Multi-condition definition via trigger sequences
- In-place baseline correction
- Seamless conversion to ERP averages

### Capabilities

- **Flexible Extraction**: Define windows relative to trigger onset (e.g., -200ms to +1000ms).
- **Condition Matching**: Match triggers or sequences of triggers to specific condition names.
- **Baseline Correction**: Subtract mean activity from a pre-stimulus interval to handle slow drifts.
- **Visualization**: Specialized plots for trial-by-trial data and ERP averages.

## Workflow Summary

1. **Data Preparation**: Load continuous data, apply layout, and perform basic preprocessing (filtering, re-referencing).
2. **Define Conditions**: Use `EpochCondition` to specify which triggers belong to which experimental condition.
3. **Extraction**: Use `extract_epochs()` to create segments. This returns a collection of `EpochedData` objects.
4. **Baseline Correction**: Correct for voltage offsets using `baseline!()`.
5. **Averaging and Comparison**: Generate ERPs by averaging across trials and plot them to compare conditions.


## Code Examples

::: details Show Code

```julia
using EegFun

# Read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# Read and prepare layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# Create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

# Minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Define simple epoch conditions
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
]

# Extract epochs (-200ms to 1000ms around triggers)
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0));

# Plot individual epochs for first condition for three channels
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1, :Cz, :Oz]))

# Plot different channel selections
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Fp1]))
EegFun.plot_epochs(epochs[2], channel_selection = EegFun.channels([:Cz]))

# Baseline correction to stimulus onset (t=0)
EegFun.baseline!(epochs, (-0.2, 0.0))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels([:Cz]))

# Compare both conditions on same plot
erps = EegFun.average_epochs(epochs)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Cz]))
```

:::

## See Also

- [API Reference](../reference/index.md)
