This demo shows how to save and load processed EEG data using the JLD2 format.

### Why JLD2?

JLD2 is the standard approach for persisting Julia objects. When used with EegFun.jl, it preserves all data types and metadata without conversion, so you can save and reload continuous data, epochs, ERPs, and ICA results exactly as they are.

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

**Direct file paths** — many EegFun plot functions accept JLD2 file paths directly, which is convenient for quick visualization without loading data:

```julia
EegFun.plot_databrowser("continuous_data.jld2")
EegFun.plot_databrowser("continuous_data.jld2", "ica.jld2")
```

## Workflow Summary

This demo covers:

### 1. Basic Saving and Loading

- Save and reload continuous data
- Verify round-trip consistency

### 2. Epoched Data and ERPs

- Save all conditions in one file or individually
- Save averaged ERPs

### 3. ICA Results

- Persist ICA decompositions for reuse across sessions
- Load and apply to new data

### 4. Multiple Items per File

- Store related objects (epochs, ERPs, ICA) in a single file
- Load individual keys by name

### 5. Batch Loading

- Use file patterns to load a cohort of participants
- Feed into group-level analysis

### 6. File Organization

- Recommended project directory structure
- Naming conventions for easy batch processing
