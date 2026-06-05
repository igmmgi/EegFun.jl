# Channel Operations

This demo shows how to manipulate channels by averaging groups together, computing differences between them, and deleting unwanted channels.

### Key Functions

| Function | Purpose | Example |
| --- | --- | --- |
| `channel_difference!` | Subtract one channel group from another | EOG calculation |
| `calculate_eog_channels!` | Calculate both vEOG and hEOG from config | Convenience wrapper |
| `channel_average!` | Average channel groups into new columns | ROI means |
| `channel_delete!` | Remove channels from data and layout | Drop non-EEG channels |

### Common Use Cases

- **EOG channels**: `channel_difference!` calculates vEOG (Fp1/Fp2 minus IO1/IO2) and hEOG (F9 minus F10)
- **ROI averages**: `channel_average!` creates region-of-interest means (e.g., frontal, parietal)
- **Cleanup**: `channel_delete!` removes non-EEG channels (e.g., photodiode, reference)

## Workflow Summary

### Channel Difference

- Calculate EOG channels from electrode pairs
- Create arbitrary difference channels

### EOG Configuration

- Use `EogConfig` to calculate both EOG channels in one call

### Channel Averaging

- Average channel groups with custom labels
- Reduce dataset to only averaged channels

### Channel Deletion

- Remove single or multiple channels (mutating and non-mutating versions)


## Code Examples

::: details Show Code

```julia
# Demo: Channel Operations
# Shows how to average, delete, and compute differences between channels.
# Covers channel_average!, channel_delete!, and channel_difference!.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

#######################################################################
# LOAD DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)

dat = EegFun.create_eegfun_data(dat, layout)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)


#######################################################################
# CHANNEL DIFFERENCE (EOG CALCULATION)
#######################################################################

# Vertical EOG: average of upper channels minus average of lower channels
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
)
# Horizontal EOG: left minus right
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
)



#######################################################################
# EOG WITH CONFIGURATION
#######################################################################

# Calculate both EOG channels from a configuration
eog_cfg = EegFun.EogConfig(
    vEOG_criterion = 50.0,
    hEOG_criterion = 30.0,
    vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
)
EegFun.calculate_eog_channels!(dat, eog_cfg)


#######################################################################
# CHANNEL AVERAGING
#######################################################################

# Average all channels into a single column
EegFun.channel_average!(dat)

dat.data[!, :avg] # average of all channels

# Average specific channel groups (average of Fp1 and Fp2, average of O1 and O2 with 
# default output labels of :Fp1_Fp2 and :O1_O2)
EegFun.channel_average!(dat, channel_selections = [EegFun.channels([:Fp1, :Fp2]), EegFun.channels([:O1, :O2])])

# Average with custom output labels
EegFun.channel_average!(
    dat,
    channel_selections = [EegFun.channels([:Fp1, :Fp2]), EegFun.channels([:O1, :O2])],
    output_labels = [:avg1, :avg2],
)

# Reduce to just the averaged channels (drops originals)
dat_reduced = EegFun.channel_average(dat, channel_selections = [EegFun.channels([:Fp1, :Fp2])], reduce = true)


#######################################################################
# CHANNEL DELETION
#######################################################################

# Delete a single channel
EegFun.channel_delete!(dat, :avg)
EegFun.channel_delete!(dat, [:Fp1_Fp2, :O1_O2])

# Delete multiple channels
EegFun.channel_delete!(dat, [:vEOG, :hEOG])

# Non-mutating version (returns a copy)
dat_subset = EegFun.channel_delete(dat, [:Fp1, :Fp2])
```

:::

## See Also

- [API Reference](../../reference/index.md)
