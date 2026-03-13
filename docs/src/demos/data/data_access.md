# Data Access

This demo shows how to inspect and navigate EegFun data structures using common access utilities.

### Data Access Functions

EegFun provides several functions to access different parts of a data structure:

- **`all_data`**: Returns the complete DataFrame (metadata + channels + extra columns)
- **`meta_data`**: Returns only metadata columns (time, sample, triggers)
- **`channel_data`**: Returns only EEG channel columns
- **`extra_data`**: Returns derived/extra columns (EOG flags, artifact markers)
- **`channel_labels`**: Returns the channel names as a vector of Symbols

### Quick Preview

For quick data inspection without viewing the entire dataset:

- **`head(dat)`**: Shows the first N rows (default: 5)
- **`tail(dat)`**: Shows the last N rows (default: 5)
- **`viewer(dat)`**: Opens data in VS Code's table viewer (falls back to console display)

### Works Across Data Types

All access functions work consistently across ContinuousData, EpochData, and ErpData. For EpochData, you can add an `epoch_selection` parameter to access specific epochs.

## Workflow Summary

This demo covers:

- Accessing all data, metadata, channel data, and extra columns
- Using head and tail for quick data preview
- Using viewer for VS Code integration
- Accessing data with epoch selection
- Querying data dimensions (n_epochs, n_values)


## Code Examples

::: details Show Code

```julia
# Demo: Data Access and Inspection
# Shows how to inspect and navigate EegFun data structures using
# head, tail, viewer, and other access utilities.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

#######################################################################
# LOAD SOME DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)


#######################################################################
# DATA ACCESS FUNCTIONS
#######################################################################

# Get all data as a DataFrame
EegFun.all_data(dat)

# Get only metadata columns (time, sample, triggers, etc.)
EegFun.meta_data(dat)

# Get only EEG channel columns
EegFun.channel_data(dat)

# Get extra/derived columns (EOG flags, etc.) — empty until calculated
EegFun.extra_data(dat)

# Get specific channel labels
EegFun.channel_labels(dat)

# Get the filename
EegFun.file_name(dat)


#######################################################################
# HEAD AND TAIL — Quick Data Preview
#######################################################################

# View the first rows (default: 5)
EegFun.head(dat)

# View the first 10 rows
EegFun.head(dat, n = 10)

# View the last rows (default: 5)
EegFun.tail(dat)

# View the last 20 rows
EegFun.tail(dat, n = 20)


#######################################################################
# VIEWER — VS Code Integration
#######################################################################

# If in VS Code, opens data in the table viewer; otherwise prints to console
EegFun.viewer(EegFun.all_data(dat))


#######################################################################
# DATA ACCESS WITH EPOCHS
#######################################################################

# Create some epochs
epoch_cfg = [EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Access with epoch selection
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:3))
EegFun.channel_data(epochs, epoch_selection = EegFun.epochs(1))

# Head and tail work on all data types
EegFun.head(epochs[1])
EegFun.tail(epochs[1])

# Get number of epochs, channels, etc.
EegFun.n_epochs(epochs[1])
EegFun.n_channels(epochs[1])
```

:::

## See Also

- [API Reference](../../reference/index.md)
