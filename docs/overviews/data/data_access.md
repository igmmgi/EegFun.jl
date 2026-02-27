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
