# Realignment

This demo shows how to realign stimulus-locked epochs to a different time point, such as response time, for response-locked ERP analysis.

### When to Use Realignment

Stimulus-locked epochs are aligned to stimulus onset (t=0). Realignment shifts t=0 to a different event — typically the participant's response — so that you can study activity time-locked to that event instead.

Common use cases:

- **Response-locked ERPs** — study motor preparation relative to button press
- **Saccade-locked ERPs** — study activity relative to eye movement onset
- **Any event-locked analysis** — realign to any column in the epoch data

### How it Works

1. Each epoch's time vector is shifted so that the realignment value becomes t=0
2. All epochs are cropped to the common time interval that is valid across all trials
3. A uniform time vector is regenerated to ensure consistency

### Key Functions

| Function | Purpose |
| --- | --- |
| `realign!(epochs, :rt)` | Realign in place (mutating) |
| `realign(epochs, :rt)` | Return a realigned copy |
| `realign(file_pattern, :rt)` | Batch realign across participants |

## Workflow Summary

### 1. Single-Participant Realignment

- Realign epochs to response time column

### 2. Batch Realignment

- Process all participant files in a directory

### 3. Typical Pipeline

- Extract stimulus-locked epochs → realign to RT → average → LRP → jackknife


## Code Examples

::: details Show Code

```julia
# Demo: Response-Locked Realignment
# Shows how to realign stimulus-locked epochs to a different time point
# (e.g., response time) for response-locked ERP analysis.

using EegFun

#######################################################################
# 1. SINGLE-PARTICIPANT REALIGNMENT
#######################################################################

# Load stimulus-locked epoched data
# epochs = EegFun.read_data("participant_1_epochs.jld2")

# Realign to response times stored in the :rt column
# EegFun.realign!(epochs, :rt)

# Non-mutating version (returns a new copy)
# realigned = EegFun.realign(epochs, :rt)


#######################################################################
# 2. BATCH REALIGNMENT
#######################################################################

# Process all participant epoch files in a directory.
# Each epoch must have a column with the realignment time
# (e.g., response time relative to stimulus onset).

# EegFun.realign("epochs_cleaned", :rt)

# Specify input directory
# EegFun.realign("epochs_cleaned", :rt, input_dir = "/path/to/epochs/")

# Specific participants
# EegFun.realign("epochs_cleaned", :rt, participant_selection = EegFun.participants([1, 2, 3]))


#######################################################################
# 3. TYPICAL WORKFLOW
#######################################################################

# Response-locked LRP analysis:
#
# 1. Extract stimulus-locked epochs with RT column
# 2. Realign to response time:
#      realign("epochs_cleaned", :rt)
# 3. Average to ERPs:
#      average_epochs in realigned directory
# 4. Calculate LRP:
#      lrp("realigned_erps", [(1, 2)])
# 5. Jackknife average:
#      jackknife_average("lrp")
```

:::

## See Also

- [API Reference](../../reference/index.md)
