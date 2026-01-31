
# How to Create Epochs {#How-to-Create-Epochs}

Extract time-locked segments of EEG data around experimental events.

## Basic Epoch Extraction {#Basic-Epoch-Extraction}

Extract epochs around a single trigger code:

```julia
using EegFun

# Define condition
condition = EegFun.EpochCondition(
    name = "Target",
    trigger_sequences = [[1]]
)

# Extract epochs from -200ms to 800ms
epochs = EegFun.extract_epochs(
    dat,
    1,           # TODO: Clarify this parameter
    condition,
    -0.2,        # Pre-stimulus time (seconds)
    0.8          # Post-stimulus time (seconds)
)
```


## Multiple Trigger Codes {#Multiple-Trigger-Codes}

Combine multiple trigger codes into one condition:

```julia
# Respond to any of triggers 1, 2, or 3
condition = EegFun.EpochCondition(
    name = "AllStimuli",
    trigger_sequences = [[1, 2, 3]]
)
```


## Complex Trigger Sequences {#Complex-Trigger-Sequences}
> 
> **TODO**: Add examples of:
> - Sequential trigger patterns
>   
> - Time-constrained sequences
>   
> - Conditional triggering
>   
> 


## Baseline Correction {#Baseline-Correction}

```julia
# TODO: Add baseline correction example
# epochs_baseline = EegFun.baseline_correct(epochs, (-0.2, 0.0))
```


## Epoch Rejection {#Epoch-Rejection}

Mark and handle bad epochs:

```julia
# TODO: Add epoch rejection examples
# - Amplitude-based rejection
# - Channel-based rejection
# - Manual rejection
```


## Extracting Metadata {#Extracting-Metadata}

Access epoch metadata:

```julia
# TODO: Show how to access trigger info, timing, etc.
```


## Common Patterns {#Common-Patterns}

### Pattern 1: Standard ERP Epochs {#Pattern-1:-Standard-ERP-Epochs}

```julia
# 200ms baseline, 1000ms post-stimulus
epochs = EegFun.extract_epochs(dat, 1, condition, -0.2, 1.0)
```


### Pattern 2: Long Time-Frequency Epochs {#Pattern-2:-Long-Time-Frequency-Epochs}

```julia
# Longer epochs for TF analysis
epochs = EegFun.extract_epochs(dat, 1, condition, -0.5, 2.0)
```


## See Also {#See-Also}
- [ERP analysis tutorial](../tutorials/erp-analysis.md) - Complete ERP workflow
  
- [Epoch data structures](../explanations/data-structures.md) - Understand epoch organization
  
