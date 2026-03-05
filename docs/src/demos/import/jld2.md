# Data Persistence (JLD2)

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


## Code Examples

::: details Show Code

```julia
# Demo: Data Persistence with JLD2
# Shows how to save and load processed EEG data using JLD2 format,
# including EegFun's built-in loading utilities and file organization.

using EegFun
using JLD2


#######################################################################
# BASIC SAVING AND LOADING
#######################################################################

# Read and prepare some data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Some preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Save continuous data
jldsave("continuous_data.jld2"; data = dat)

# Load it back
loaded_dat = load("continuous_data.jld2", "data")

# Verify it's the same
EegFun.all_data(loaded_dat) == EegFun.all_data(dat)


#######################################################################
# SAVING EPOCHED DATA
#######################################################################

# Extract epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "Condition1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Condition2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))

# Save all conditions
jldsave("epochs_all.jld2"; data = epochs)

# Save individual conditions
jldsave("epochs_condition1.jld2"; data = epochs[1])
jldsave("epochs_condition2.jld2"; data = epochs[2])

# Load epochs back
loaded_epochs = load("epochs_all.jld2", "data")


#######################################################################
# SAVING ERPs
#######################################################################

# Create ERPs
erps = EegFun.average_epochs(epochs)

# Save ERPs
jldsave("erps.jld2"; data = erps)

# Load ERPs
loaded_erps = load("erps.jld2", "data")


#######################################################################
# SAVING ICA RESULTS
#######################################################################

# Run ICA
EegFun.is_extreme_value!(dat, 100)
ica_result = EegFun.run_ica(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))

# Save ICA decomposition
jldsave("ica_decomposition.jld2"; data = ica_result)

# Load ICA
loaded_ica = load("ica_decomposition.jld2", "data")

# Can use loaded ICA with databrowser
EegFun.plot_databrowser(dat, loaded_ica)


#######################################################################
# SAVING MULTIPLE ITEMS IN ONE FILE
#######################################################################

# Save multiple related objects together
jldsave("analysis_results.jld2"; epochs = epochs, erps = erps, ica = ica_result)

# Load specific items
my_epochs = load("analysis_results.jld2", "epochs")
my_erps = load("analysis_results.jld2", "erps")

# Or load everything
results = load("analysis_results.jld2")
results["epochs"]
results["erps"]
results["ica"]


#######################################################################
# EEGFUN DATA LOADING UTILITIES
#######################################################################

# read_data — smart loader that auto-detects EegFun data types from JLD2 files
# Returns the data directly (single variable) or a Dict (multiple variables)
loaded = EegFun.read_data("erps.jld2")

# read_all_data — batch load all files matching a pattern from a directory
# Useful for loading an entire cohort for group analysis
# all_erps = EegFun.read_all_data("derivatives/erps/erps")

# With participant selection
# all_erps = EegFun.read_all_data("derivatives/erps/erps", EegFun.participants([1, 2, 3]))

# group_by_condition — organise loaded data by condition number
# Returns an OrderedDict{Int, Vector{ErpData}}
# grouped = EegFun.group_by_condition(all_erps)


#######################################################################
# BATCH LOADING FROM DIRECTORY
#######################################################################

# Suppose you have multiple participant files:
# - participant_01_epochs.jld2
# - participant_02_epochs.jld2
# - participant_03_epochs.jld2

# List all matching files
using Glob
epoch_files = glob("participant_*_epochs.jld2", "my_output_directory/")

# Load all participants
all_participant_epochs = [load(file, "data") for file in epoch_files]

# Now you can do group analysis
# grand_avg = EegFun.grand_average(all_participant_epochs)


#######################################################################
# READING DATA IN PLOT FUNCTIONS
#######################################################################

# Many EegFun plot functions accept file paths directly
EegFun.plot_databrowser("continuous_data.jld2")
EegFun.plot_databrowser("epochs_all.jld2")
EegFun.plot_databrowser("continuous_data.jld2", "ica_decomposition.jld2")

# This is convenient for quick visualization without loading into memory


#######################################################################
# FILE ORGANIZATION
#######################################################################

# Best practice: organize outputs by processing stage
# 
# project/
# ├── data/
# │   └── raw/              # Original BDF/VHDR files
# ├── derivatives/
# │   ├── continuous/       # Filtered, rereferenced continuous data
# │   ├── ica/              # ICA decompositions
# │   ├── epochs/           # Epoched data
# │   ├── epochs_clean/     # Epochs after artifact rejection
# │   └── erps/             # Averaged ERPs
# └── results/
#     └── grand_averages/   # Group-level results

# Example: saving with organized structure
output_dir = mktempdir()
participant_id = "sub-01"
jldsave(joinpath(output_dir, "$(participant_id)_epochs.jld2"); data = epochs)


#######################################################################
# NAMING CONVENTIONS
#######################################################################

# Use consistent naming for easy batch processing:
# - Descriptive prefixes: "sub-", "participant_", etc.
# - Stage indicators: "_continuous", "_epochs", "_erps", "_clean"
# - Zero-padding for sorting: "sub-01" not "sub-1"

# Examples:
# jldsave("sub-01_continuous.jld2"; data = dat)
# jldsave("sub-01_epochs_raw.jld2"; data = epochs)
# jldsave("sub-01_epochs_clean.jld2"; data = epochs_cleaned)
# jldsave("sub-01_erps.jld2"; data = erps)


#######################################################################
# MEMORY CONSIDERATIONS
#######################################################################

# For large datasets, load only what you need
# Instead of: results = load("big_file.jld2")
# Use: epochs = load("big_file.jld2", "epochs")

# For very large cohorts, process one participant at a time
# rather than loading all into memory at once
```

:::

## See Also

- [API Reference](../../reference/index.md)
