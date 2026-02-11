# Plot RSA

This demo shows how to visualise Representational Similarity Analysis (RSA) results.

### What is RSA?

- **RDM** (Representational Dissimilarity Matrix) captures pairwise differences between neural representations
- **RDM Timecourse** shows how representational structure evolves over time
- **Model Correlations** compare observed RDMs against theoretical models

### Key Functions

| Function | Purpose |
| --- | --- |
| `plot_rdm_heatmap(rsa_result)` | Display an RDM as a heatmap |
| `plot_rdm_timecourse(rsa_result)` | Dissimilarity over time |
| `plot_model_correlations(rsa_result)` | Model comparison over time |

### Key Parameters

| Parameter | Function | Description |
| --- | --- | --- |
| `time_point` | `plot_rdm_heatmap` | Specific time (seconds or index); `nothing` = average |
| `condition_pairs` | `plot_rdm_timecourse` | Which pairs to plot; `:all` = all |
| `show_colorbar` | All | Show/hide colour bar |
| `colormap` | All | Colour map for the heatmap |

### What You'll Learn

1. Plotting RDM heatmaps at specific time points or averaged
2. Visualising dissimilarity timecourses for selected condition pairs
3. Comparing neural RDMs against theoretical model predictions


## Code Examples

::: details Show Code

```julia
# Demo: Plotting RSA Results
# Shows how to visualise Representational Similarity Analysis results:
# RDM heatmaps, dissimilarity timecourses, and model correlations.

using EegFun

#######################################################################
# 1. LOAD DATA AND COMPUTE RSA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.1)
EegFun.lowpass_filter!(dat, 30.0)

epochs = EegFun.epoch_data(dat, [:trigger1, :trigger2], (-0.2, 0.8))
EegFun.baseline!(epochs, (-0.2, 0.0))

# compute RSA
rsa_result = EegFun.rsa(epochs)

#######################################################################
# 2. RDM HEATMAP — AVERAGE ACROSS TIME
#######################################################################

# visualise overall representational structure
EegFun.plot_rdm_heatmap(rsa_result)

#######################################################################
# 3. RDM HEATMAP — SPECIFIC TIME POINT
#######################################################################

# RDM at 300 ms post-stimulus
EegFun.plot_rdm_heatmap(rsa_result, time_point = 0.3)

# RDM at time index 50
EegFun.plot_rdm_heatmap(rsa_result, time_point = 50)

#######################################################################
# 4. DISSIMILARITY TIMECOURSE
#######################################################################

# plot dissimilarity over time for all condition pairs
EegFun.plot_rdm_timecourse(rsa_result)

# only specific pairs
EegFun.plot_rdm_timecourse(rsa_result, condition_pairs = [(1, 2), (1, 3)])

#######################################################################
# 5. MODEL CORRELATIONS
#######################################################################

# compare RSA results against theoretical models
# (requires having run compare_models beforehand)
# EegFun.plot_model_correlations(rsa_result)
```

:::

## See Also

- [API Reference](../../reference/index.md)
