# RSA

RSA workflow for analyzing representational geometries across conditions and time.

## Overview

Demonstrates RSA workflow for analyzing representational geometries across conditions and time.

## Source Code

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

## See Also

- [API Reference](../reference/index.md)
