# Epoch Extraction

This demo explores the core functionality of segmenting continuous data into epochs (or trials) based on event markers, and basic manipulation of these segments.

### What are Epochs?

Epoching is the process of extracting specific time intervals around events (e.g., stimuli or responses) from a continuous recording. This allows for epoch-based analysis and averaging to reveal Event-Related Potentials (ERPs).

**Key features of EegFun's epoching**:

- Time-relative segmentation
- Multi-condition definition via trigger sequences

### Capabilities

- **Flexible Extraction**: Define intervals relative to trigger onset (e.g., -200ms to +1000ms).
- **Condition Matching**: Match triggers or sequences of triggers to specific condition names.

## Workflow Summary

1. **Data Preparation**: Load continuous data, apply layout, and perform basic preprocessing (filtering, re-referencing).
2. **Define Conditions**: Use `EpochCondition` to specify which triggers belong to which experimental condition.
3. **Extraction**: Use `extract_epochs()` to create segments. This returns a collection of `EpochData` objects.


## Code Examples

::: details Show Code

```julia
# Demo: Epoch Extraction
# Shows how to extract epochs from continuous data using trigger sequences.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# Read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# Read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# Create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Define simple epoch condition
epoch_cfg = EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]])
# Extract epochs (-200ms to 1000ms around triggers)
epoch = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0));

# Or multiple epoch conditions
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

- [API Reference](../../reference/index.md)
