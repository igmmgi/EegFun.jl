# Selection Helpers

This demo shows how to use EegFun's selection helper functions for filtering, subsetting, and targeting specific parts of your data.

### Selection Helpers

EegFun uses predicate-generator functions that create selection criteria:

- **`channels()`**: Select channels by name, vector, or regex pattern
- **`times()`**: Select time windows by start/end in seconds
- **`epochs()`**: Select epochs by index or range
- **`samples()`**: Select samples using custom predicates on metadata columns

### How They Work

Selection helpers return functions (predicates) that are passed to `subset()`, and other functions via keyword arguments:

```julia
# channels(:Cz, :Pz) returns a function that selects those channels
subset(erp, channel_selection = channels(:Cz, :Pz))
```

### Composing Selections

Multiple selections can be combined in a single call:

```julia
subset(epochs[1],
    channel_selection = channels(:Cz),
    epoch_selection = epochs(1:3),
    interval_selection = times(0.0, 0.5),
)
```

This pattern is consistent across subsetting, plotting, and analysis functions.

## Workflow Summary

This demo covers:

- Channel selection by name, vector, and regex pattern
- Time window selection with `times(start, end)`
- Epoch selection with `epochs(range)`
- Sample-level predicates with `samples()`
- Combining multiple selections in subset and plotting calls


## Code Examples

::: details Show Code

```julia
# Demo: Selection Helpers
# Shows how to use EegFun's selection helper functions (channels, conditions,
# participants, epochs, times, samples) for filtering and subsetting data.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

#######################################################################
# SETUP
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Create epochs and ERPs
epoch_cfg =
    [EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]), EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
erps = EegFun.average_epochs(epochs)


#######################################################################
# CHANNELS() — Channel Selection
#######################################################################

# Select all channels (default)
EegFun.subset(erps[1], channel_selection = EegFun.channels())

# Select specific channels by name
EegFun.subset(erps[1], channel_selection = EegFun.channels(:Cz, :Pz, :Fz))

# Select channels by vector
EegFun.subset(erps[1], channel_selection = EegFun.channels([:Fp1, :Fp2]))



#######################################################################
# TIMES() — Time Window Selection
#######################################################################

# Select all time points (default)
EegFun.subset(erps[1], interval_selection = EegFun.times())

# Select a time window (in seconds)
EegFun.subset(erps[1], interval_selection = EegFun.times(0.0, 0.5))

# Combine channel and time selection
sub = EegFun.subset(erps[1], channel_selection = EegFun.channels(:Cz, :Pz), interval_selection = EegFun.times(0.1, 0.4))


#######################################################################
# EPOCHS() — Epoch Selection
#######################################################################

# Select all epochs (default)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs())

# Select specific epoch indices
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:5))

# Subset epochs from an EpochData object
EegFun.subset(epochs[1], epoch_selection = EegFun.epochs(1:3))


#######################################################################
# SAMPLES() — Sample-level Predicates
#######################################################################

# Custom predicate on metadata columns
EegFun.subset(dat, sample_selection = x -> x.sample .<= 5000)

# Time-based predicate
EegFun.subset(dat, sample_selection = x -> x.time .<= 5.0)


#######################################################################
# COMBINING SELECTIONS
#######################################################################

# Multiple selections at once
sub = EegFun.subset(
    epochs[1],
    channel_selection = EegFun.channels(:Cz),
    epoch_selection = EegFun.epochs(1:3),
    interval_selection = EegFun.times(0.0, 0.5),
)

# Use with plotting
EegFun.plot_erp(erps, channel_selection = EegFun.channels(:Cz, :Pz), interval_selection = EegFun.times(-0.1, 0.8))
```

:::

## See Also

- [API Reference](../../reference/index.md)
