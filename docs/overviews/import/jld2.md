This demo shows how to save and load processed EEG data using the JLD2 format and EegFun's built-in loading utilities.

### Why JLD2?

JLD2 is the standard approach for persisting Julia objects. When used with EegFun.jl, it preserves all data types and metadata without conversion, so you can save and reload continuous data, epochs, ERPs, and ICA results exactly as they are.

### Key Functions

| Function | Purpose |
| --- | --- |
| `jldsave` | Save one or more objects to a JLD2 file |
| `load` | Load specific keys or all data from a JLD2 file |
| `read_data` | Smart EegFun loader — auto-detects data types from JLD2 |
| `read_all_data` | Batch load all matching files from a directory |
| `group_by_condition` | Organise loaded data into an OrderedDict by condition |

### Key Patterns

**Saving data** — use Julia's `jldsave` to store one or more objects:

```julia
jldsave("epochs.jld2"; data = epochs)
jldsave("results.jld2"; epochs = epochs, erps = erps, ica = ica_result)
```

**Loading data** — use `load` to read specific keys:

```julia
epochs = load("epochs.jld2", "data")
results = load("results.jld2")    # loads everything as a Dict
```

**EegFun utilities** — smart loading and batch operations:

```julia
loaded = EegFun.read_data("erps.jld2")                    # auto-detect types
all_erps = EegFun.read_all_data("erps", "derivatives/erps/")  # load all matching
grouped = EegFun.group_by_condition(all_erps)               # organise by condition
```

**Direct file paths** — many EegFun plot functions accept JLD2 file paths directly:

```julia
EegFun.plot_databrowser("continuous_data.jld2")
EegFun.plot_databrowser("continuous_data.jld2", "ica.jld2")
```

## Workflow Summary

### Basic Saving and Loading

- Save and reload continuous data, epochs, ERPs, and ICA results
- Verify round-trip consistency

### Multiple Items per File

- Store related objects (epochs, ERPs, ICA) in a single file
- Load individual keys by name

### EegFun Loading Utilities

- `read_data` for smart single-file loading with auto type detection
- `read_all_data` for batch loading a cohort with participant selection
- `group_by_condition` for organising loaded data by condition number

### Batch Loading

- Use file patterns to load a cohort of participants
- Feed into group-level analysis (e.g., grand average)

### File Organization

- Recommended project directory structure
- Naming conventions for easy batch processing
- Memory considerations for large datasets
