This demo shows how to save and load processed EEG data using `EegFun.save_data`, `EegFun.read_data`, and the JLD2 format.

### Why JLD2?

JLD2 is the standard approach for persisting Julia objects. When used with EegFun.jl, it preserves all data types and metadata without conversion, so you can save and reload continuous data, epochs, ERPs, ICA decompositions, and statistical results as-is.

### Key Functions

| Function | Purpose |
| --- | --- |
| `save_data` | Save an `EegFunData` object (or vector) to a `.jld2` file with format metadata |
| `read_data` | Smart EegFun loader — auto-detects data types from a `.jld2` file |
| `read_all_data` | Batch load all matching files from a directory |
| `group_by_condition` | Organise loaded data into an OrderedDict by condition |
| `JLD2.jldsave` / `JLD2.load` | Lower-level JLD2 access for multi-key custom bundles |

### Key Patterns

**Saving and loading EegFun data** — use `EegFun.save_data` and `EegFun.read_data` (no `using JLD2` needed):

```julia
EegFun.save_data("epochs.jld2", epochs)
loaded_epochs = EegFun.read_data("epochs.jld2")
```

**Batch loading** — load a cohort across participants and group by condition:

```julia
all_erps = EegFun.read_all_data("./derivatives/erps/erps_final")
grouped = EegFun.group_by_condition(all_erps)
```

**Custom multi-key bundles** — use `JLD2.jldsave` and `JLD2.load` directly when storing several named keys in one file:

```julia
using JLD2
jldsave("results.jld2"; epochs = epochs, erps = erps, ica = ica_result)
results = load("results.jld2")    # loads everything as a Dict
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
