
# ERP Analysis Tutorial {#ERP-Analysis-Tutorial}

Learn how to extract event-related potentials (ERPs) from continuous EEG data.

## Overview {#Overview}

ERP analysis involves:
1. Defining epoch conditions based on experimental events
  
2. Extracting time-locked epochs around events
  
3. Averaging epochs to compute ERPs
  
4. Analyzing and visualizing ERP components
  

## Defining Epoch Conditions {#Defining-Epoch-Conditions}

Create conditions based on trigger codes:

```julia
using EegFun

# Define a simple condition
condition = EegFun.EpochCondition(
    name = "Target",
    trigger_sequences = [[1]]  # Trigger code 1
)

# Multiple trigger codes
condition = EegFun.EpochCondition(
    name = "AllStimuli", 
    trigger_sequences = [[1, 2, 3]]
)
```

> 
> **TODO**: Add complex trigger sequence examples
> 


## Extracting Epochs {#Extracting-Epochs}

Extract time-locked segments around events:

```julia
# Extract epochs from -200ms to 800ms around triggers
epochs = EegFun.extract_epochs(
    dat,           # Preprocessed continuous data
    1,             # Sampling rate (TODO: clarify this parameter)
    condition,     # Epoch condition
    -0.2,          # Start time (seconds)
    0.8            # End time (seconds)
)
```


## Computing ERPs {#Computing-ERPs}

Average epochs to obtain ERPs:

```julia
# Average all epochs for this condition
erps = EegFun.average_epochs(epochs)
```


## Visualization {#Visualization}

Plot ERP waveforms:

```julia
using GLMakie

# Plot ERP waveforms
fig, ax = EegFun.plot_erp(erps)

# TODO: Add topographic map example
# TODO: Add channel selection example
```


## Extracting ERP Measurements {#Extracting-ERP-Measurements}

Measure peak amplitudes and latencies:

```julia
# TODO: Add erp_measurements example
# measurements = EegFun.erp_measurements(
#     erps,
#     analysis_interval = (0.3, 0.5),
#     ...
# )
```


## Complete Example {#Complete-Example}

```julia
# TODO: Add complete end-to-end example
```


## Next Steps {#Next-Steps}
- Learn about [epoch creation](../how-to/create-epochs.md) in detail
  
- Understand [ERP measurement methods](../explanations/erp-measurements.md)
  
- Explore [statistical testing](statistical-analysis.md)
  
