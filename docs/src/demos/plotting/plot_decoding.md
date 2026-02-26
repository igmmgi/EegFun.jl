# Plot Decoding

This demo shows how to visualise MVPA (multivariate pattern analysis) decoding results.

### Why Plot Decoding Results?

- **Temporal dynamics** — see when brain signals become discriminative
- **Error shading** — visualise variability across cross-validation folds or subjects
- **Significance markers** — overlay statistical test results on the accuracy curve

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_decoding(decoded)` | Plot single-subject accuracy | Quick inspection |
| `plot_decoding(decoded_list)` | Multi-subject subplot grid | Individual differences |
| `plot_decoding(decoded, stats)` | Accuracy + significance | Publication figure |

### What You'll Learn

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

using EegFun

#######################################################################
# LOAD DATA AND PREPARE FOR DECODING
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# basic preprocessing
EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

# epoch
epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

#######################################################################
# RUN DECODING
#######################################################################

# decode condition from EEG patterns
decoded = EegFun.decode(epochs)

#######################################################################
# BASIC DECODING PLOT
#######################################################################

# plot accuracy over time with error shading
EegFun.plot_decoding(decoded)

#######################################################################
# CUSTOM STYLING
#######################################################################

# change colour, line width, title
EegFun.plot_decoding(decoded,
    color = :red,
    linewidth = 3,
    title = "Face vs. Object Decoding"
)

# hide error shading
EegFun.plot_decoding(decoded, show_error = false)

#######################################################################
# MULTI-SUBJECT SUBPLOT GRID
#######################################################################

# when you have a list of decoded results (one per subject),
# each is plotted in its own subplot
# EegFun.plot_decoding(all_decoded, title = "Individual Subjects")

#######################################################################
# DECODING WITH SIGNIFICANCE
#######################################################################

# after running statistics on the decoded data:
# stats = EegFun.test_against_chance(decoded_list, alpha = 0.05,
#              correction_method = :bonferroni)
# EegFun.plot_decoding(grand_avg, stats, show_significance = true)
```

:::

## See Also

- [API Reference](../../reference/index.md)
