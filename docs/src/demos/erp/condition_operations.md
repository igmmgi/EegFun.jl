# Condition Operations

This demo shows how to combine epoch conditions, compute ERP difference waves, and create condition average waves.

### Key Functions

| Function | Purpose | Input |
| --- | --- | --- |
| `condition_combine` | Merge epoch conditions into new groups | Epoch data |
| `condition_difference` | Subtract one condition's ERP from another | ERP data |
| `condition_average` | Average multiple conditions into one ERP | ERP data |

### When to Use

- **`condition_combine`**: When you have many conditions and want to pool some together before averaging (e.g., merge "Go left" and "Go right" into a single "Go" condition)
- **`condition_difference`**: When you need difference waves for analysis (e.g., target minus standard, congruent minus incongruent)
- **`condition_average`**: When you want to average across conditions after computing ERPs (e.g., collapsing across irrelevant factors)

### Important Notes

- `condition_combine` works on epoch data (before averaging); `condition_difference` and `condition_average` work on ERP data (after averaging)
- All three also have batch versions that process all matching JLD2 files in a directory

## Workflow Summary

### Condition Combining

- Merge multiple epoch conditions into new groups (operates on epochs, before averaging)

### Condition Differences

- Create difference waves from ERP condition pairs (e.g., `[(1, 2), (3, 4)]`)

### Condition Averaging

- Average ERP conditions into combined waves (e.g., `[[1, 2], [3, 4]]` or `[[1, 2, 3, 4]]`)

### Typical Pipeline

- Preprocess → Combine conditions → Average → Difference/Average waves → Grand average


## Code Examples

::: details Show Code

```julia
# Demo: Condition Operations
# Shows how to combine epoch conditions, compute ERP difference waves,
# and create condition average waves.
# Covers condition_combine, condition_difference, and condition_average.

using EegFun

#######################################################################
# CREATE EXAMPLE DATA
#######################################################################

# Load raw data and do minimal preprocessing
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)

dat = EegFun.create_eegfun_data(dat, layout)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Extract epochs for four conditions
epoch_cfg = [
    EegFun.EpochCondition(name = "Trigger1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Trigger2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "Trigger3", trigger_sequences = [[3]]),
    EegFun.EpochCondition(name = "Trigger4", trigger_sequences = [[5]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.baseline!(epochs, (-0.2, 0.0))


#######################################################################
# CONDITION COMBINE (EPOCHS)
#######################################################################

# Combine conditions 1+2 into one group and 3+4 into another
epochs_combined = EegFun.condition_combine(epochs, [[1, 2], [3, 4]])


#######################################################################
# CONDITION DIFFERENCE WAVES (ERPS)
#######################################################################

# Create ERPs from the original four conditions
erps = EegFun.average_epochs(epochs)

# Difference wave: condition 1 minus condition 2
diff_waves = EegFun.condition_difference(erps, [(1, 2)])
EegFun.plot_erp(diff_waves, channel_selection = EegFun.channels([:Cz]))

# Multiple difference pairs
diff_waves = EegFun.condition_difference(erps, [(1, 2), (3, 4)])


#######################################################################
# CONDITION AVERAGE WAVES (ERPS)
#######################################################################

# Average conditions 1 and 2 together, and 3 and 4 together
avg_waves = EegFun.condition_average(erps, [[1, 2], [3, 4]])

# Average all four conditions into one
avg_all = EegFun.condition_average(erps, [[1, 2, 3, 4]])
```

:::

## See Also

- [API Reference](../../reference/index.md)
