# RSA

## Overview

## Overview

This demo demonstrates Representational Similarity Analysis (RSA) for analyzing neural representational geometries.

### What is RSA?

RSA quantifies how neural patterns represent different stimuli/conditions:
1. **Compute dissimilarity**: Calculate distance between conditions (e.g., correlation distance)
2. **Create RDM**: Representational Dissimilarity Matrix showing all pairwise distances
3. **Compare geometries**: Relate neural RDMs to behavior, models, or other brain regions

### RSA Workflow

1. **Extract patterns**: Get spatial or spatiotemporal patterns for each condition
2. **Compute RDM**: Calculate dissimilarity between all condition pairs
3. **Statistical testing**: Assess significance of representational structure
4. **Model comparison**: Compare neural RDMs to theoretical predictions

### Applications

- Compare representational geometry across time
- Test computational models of cognition
- Relate neural to behavioral similarity
- Cross-modality comparisons

### Advantages

- Multivariate (uses all channels)
- Model-agnostic representational space
- Can compare across different measurement modalities


## Code Examples

::: details Show Code

```julia
using EegFun

@info EegFun.section("Representational Similarity Analysis (RSA)")
println("Basic RSA analysis...")

# Create test data for multiple participants
all_participant_epochs = Vector{Vector{EegFun.EpochData}}()
for p = 1:5
    # Use built-in test data creation function
    participant_epochs = EegFun.create_test_epoch_data_vector(
        conditions = 1:3,
        n_epochs = 100,
        n_channels = 3,
        fs = 250,
        n = 250,  # 1 second of data at 250 Hz
    )
    push!(all_participant_epochs, participant_epochs)
end

# Perform RSA for each participant using batch method
all_rsa_results = EegFun.rsa(all_participant_epochs; dissimilarity_measure = :correlation, average_trials = true)
all_rsa_results =
    EegFun.rsa(all_participant_epochs; dissimilarity_measure = :correlation, average_trials = true, normalize_method = :zscore)

# Compute grand average with noise ceiling
grand_avg_rsa = EegFun.grand_average(all_rsa_results)

# RSA Plots
EegFun.plot_rdm_heatmap(grand_avg_rsa, time_point = 0.0, title = "Grand Average RDM")
EegFun.plot_rdm_heatmap(grand_avg_rsa, title = "Grand Average RDM (time-averaged)")
EegFun.plot_rdm_timecourse(grand_avg_rsa, title = "Dissimilarity Over Time")


```

:::

## See Also

- [API Reference](../reference/index.md)
