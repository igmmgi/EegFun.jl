# Data Persistence (JLD2)

This demo shows how to save and load processed EEG data using the JLD2 format and EegFun's built-in loading utilities.

## Why JLD2?

JLD2 is the standard approach for persisting Julia objects. When used with EegFun.jl, it preserves all data types and metadata without conversion, so you can save and reload continuous data, epochs, ERPs, and ICA results exactly as they are.

## Key Functions

| Function | Purpose |
| --- | --- |
| `jldsave` | Save one or more objects to a JLD2 file |
| `load` | Load specific keys or all data from a JLD2 file |
| `read_data` | Smart EegFun loader — auto-detects data types from JLD2 |
| `read_all_data` | Batch load all matching files from a directory |
| `group_by_condition` | Organise loaded data into an OrderedDict by condition |

## Key Patterns

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
all_erps = EegFun.read_all_data("./derivatives/erps/erps_final")  # load all matching
grouped = EegFun.group_by_condition(all_erps)               # organise by condition
```

**Direct file paths** — many EegFun plot functions accept JLD2 file paths directly:

```julia
EegFun.plot_databrowser("continuous_data.jld2")
EegFun.plot_databrowser("continuous_data.jld2", "ica.jld2")
```

## Workflow Summary

## Basic Saving and Loading

- Save and reload continuous data, epochs, ERPs, and ICA results
- Verify round-trip consistency

## Multiple Items per File

- Store related objects (epochs, ERPs, ICA) in a single file
- Load individual keys by name

## EegFun Loading Utilities

- `read_data` for smart single-file loading with auto type detection
- `read_all_data` for batch loading a cohort with participant selection
- `group_by_condition` for organising loaded data by condition number

## Batch Loading

- Use file patterns to load a cohort of participants
- Feed into group-level analysis (e.g., grand average)

## File Organization

- Recommended project directory structure
- Naming conventions for easy batch processing
- Memory considerations for large datasets


## Code Examples

::: details Show Code

```julia
# Demo: Data Persistence with JLD2
# Shows how to save and load processed EEG data using JLD2 format,
# including EegFun's built-in loading utilities and file organization.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using JLD2

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)


#######################################################################
# BASIC SAVING AND LOADING
#######################################################################

# Read and prepare some data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Some preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Save continuous data
jldsave(joinpath(DEMO_OUTPUT, "continuous_data.jld2"); data = dat)

# Load it back
loaded_dat = load(joinpath(DEMO_OUTPUT, "continuous_data.jld2"), "data")

# Verify it's the same
EegFun.all_data(loaded_dat) == EegFun.all_data(dat)

EegFun.plot_databrowser(loaded_dat);


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
jldsave(joinpath(DEMO_OUTPUT, "epochs_all.jld2"); data = epochs)

# Save individual conditions
jldsave(joinpath(DEMO_OUTPUT, "epochs_condition1.jld2"); data = epochs[1])
jldsave(joinpath(DEMO_OUTPUT, "epochs_condition2.jld2"); data = epochs[2])

# Load epochs back
loaded_epochs = load(joinpath(DEMO_OUTPUT, "epochs_all.jld2"), "data")


#######################################################################
# SAVING ERPs
#######################################################################

# Create ERPs
erps = EegFun.average_epochs(epochs)

# Save ERPs
jldsave(joinpath(DEMO_OUTPUT, "erps.jld2"); data = erps)

# Load ERPs
loaded_erps = load(joinpath(DEMO_OUTPUT, "erps.jld2"), "data")


#######################################################################
# SAVING ICA RESULTS
#######################################################################

# Run ICA
EegFun.is_extreme_value!(dat, 100)
ica_result = EegFun.run_ica(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))

# Save ICA decomposition
jldsave(joinpath(DEMO_OUTPUT, "ica_decomposition.jld2"); data = ica_result)

# Load ICA
loaded_ica = load(joinpath(DEMO_OUTPUT, "ica_decomposition.jld2"), "data")

# Can use loaded ICA with databrowser
EegFun.plot_databrowser(dat, loaded_ica)


#######################################################################
# SAVING MULTIPLE ITEMS IN ONE FILE
#######################################################################

# Save multiple related objects together
jldsave(joinpath(DEMO_OUTPUT, "analysis_results.jld2"); epochs = epochs, erps = erps, ica = ica_result)

# Load specific items
my_epochs = load(joinpath(DEMO_OUTPUT, "analysis_results.jld2"), "epochs")
my_erps = load(joinpath(DEMO_OUTPUT, "analysis_results.jld2"), "erps")

# Or load everything
results = load(joinpath(DEMO_OUTPUT, "analysis_results.jld2"))
results["epochs"]
results["erps"]
results["ica"]


#######################################################################
# EEGFUN DATA LOADING UTILITIES
#######################################################################

# read_data — smart loader that auto-detects EegFun data types from JLD2 files
# Returns EegFunData, Vector{EegFunData}, or nothing
loaded = EegFun.read_data(joinpath(DEMO_OUTPUT, "erps.jld2"))

# read_all_data — batch load all files matching a pattern from a directory
# Useful for loading an entire cohort for group analysis
# all_erps = EegFun.read_all_data("derivatives/erps/erps")

# With participant selection
# all_erps = EegFun.read_all_data("derivatives/erps/erps", EegFun.participants([1, 2, 3]))

#######################################################################
# BATCH LOADING FROM DIRECTORY
#######################################################################

# Suppose you have multiple participant files in a directory:
# - sub-01_epochs.jld2
# - sub-02_epochs.jld2
# - sub-03_epochs.jld2

# Load all participants at once using EegFun's built-in batch loader
# all_participant_epochs = EegFun.read_all_data("my_output_directory/epochs")

# Or select specific participants
# selected_epochs = EegFun.read_all_data("my_output_directory/epochs", EegFun.participants([1, 2, 3]))


#######################################################################
# READING DATA IN PLOT FUNCTIONS
#######################################################################

# Many EegFun plot functions accept file paths directly
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "continuous_data.jld2"))
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "epochs_all.jld2"))
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "continuous_data.jld2"), joinpath(DEMO_OUTPUT, "ica_decomposition.jld2"))

```

:::

## See Also

- [API Reference](../../reference/index.md)
