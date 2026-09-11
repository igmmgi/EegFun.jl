# Plot Decoding

This demo shows how to visualise MVPA (multivariate pattern analysis) decoding results.

## Why Plot Decoding Results?

- **Temporal dynamics** — see when brain signals become discriminative
- **Error shading** — visualise variability across cross-validation folds or subjects
- **Significance markers** — overlay statistical test results on the accuracy curve

## Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_decoding(decoded)` | Plot single-subject accuracy | Quick inspection |
| `plot_decoding(decoded_list)` | Multi-subject subplot grid | Individual differences |
| `plot_decoding(decoded, stats)` | Accuracy + significance | Publication figure |

## What You'll Learn

1. Plotting decoding accuracy over time with error shading
2. Customising colours, line width and titles
3. Comparing subjects in a subplot grid
4. Overlaying significance markers from `test_against_chance`


## Code Examples

::: details Show Code

```julia
# Demo: Plotting Decoding Results
# Shows how to visualise MVPA decoding accuracy over time,
# compare individual subjects, and overlay statistical significance.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# LOAD DATA AND PREPARE FOR DECODING
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# basic preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# extract epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-1, 2))
EegFun.baseline!(epochs, (-0.2, 0.0))

# decode condition from EEG patterns
decoded = EegFun.decode_libsvm(epochs)

# plot accuracy over time with error shading
EegFun.plot_decoding(decoded)
```

:::

## See Also

- [API Reference](../../reference/index.md)
