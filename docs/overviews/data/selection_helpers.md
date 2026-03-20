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
