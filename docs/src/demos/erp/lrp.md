# Lateralised Readiness Potential

This demo shows how to calculate the Lateralised Readiness Potential (LRP) from ERP data.

### What is the LRP?

The LRP is an ERP measure of lateralised motor preparation. It is computed from the difference in activity between hemispheres contralateral and ipsilateral to the responding hand, and is widely used in cognitive and motor neuroscience to study response preparation timing.

### How is it Calculated?

For each channel pair (e.g., C3/C4), the LRP is:

```text
LRP_C3 = 0.5 × ((C3_right − C4_right) + (C4_left − C3_left))
```

Where "left" and "right" refer to the hand used for the response.

### Key Functions

| Function | Purpose |
| --- | --- |
| `lrp(erp_left, erp_right)` | Single-participant LRP |
| `lrp(file_pattern, pairs)` | Batch LRP across participants |

### Channel Pairing

- By default, all odd/even channel pairs are detected automatically (e.g., C3/C4, CP1/CP2)
- Use `channel_selection` to restrict to specific pairs

## Workflow Summary

### 1. Single-Participant LRP

- Calculate LRP from left and right response conditions

### 2. Batch Processing

- Process all participants with specified condition pairs

### 3. Visualisation

- Plot LRP waveforms with `plot_erp`


## Code Examples

::: details Show Code

```julia
# Demo: Lateralised Readiness Potential (LRP)
# Shows how to calculate the LRP from left-hand and right-hand ERP conditions
# and perform batch processing across participants.

using EegFun

#######################################################################
# 1. LOAD ERP DATA
#######################################################################

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")


#######################################################################
# 2. SINGLE-PARTICIPANT LRP
#######################################################################

# Calculate LRP from two conditions:
#   - dat[1] = left-hand response ERPs
#   - dat[2] = right-hand response ERPs
# Automatically pairs odd/even channels (e.g., C3/C4, CP3/CP4)

lrp_data = EegFun.lrp(dat[1], dat[2])

# Select specific channels (left hemisphere only, auto-pairs with right)
lrp_data = EegFun.lrp(dat[1], dat[2], channel_selection = EegFun.channels([:C3, :CP3]))


#######################################################################
# 3. BATCH LRP (MULTIPLE PARTICIPANTS)
#######################################################################

# Process all participant ERP files in a directory.
# condition_pairs specifies (left_condition, right_condition)

# EegFun.lrp("erps_cleaned", [(1, 2)])

# Multiple condition pairs
# EegFun.lrp("erps_cleaned", [(1, 2), (3, 4)])

# Specific channels only
# EegFun.lrp("erps_cleaned", [(1, 2)], channel_selection = EegFun.channels([:C3]))


#######################################################################
# 4. VISUALIZE
#######################################################################

EegFun.plot_erp(lrp_data)
```

:::

## See Also

- [API Reference](../../reference/index.md)
