# Selection Patterns

EegFun.jl uses **predicate generators** — functions that return selection rules — to filter channels, samples, epochs, components, and time intervals. This approach keeps analysis code readable, composable, and independent of the data at hand.

## Why Predicates?

You might wonder why we use `channels([:Cz, :Pz])` instead of just `[:Cz, :Pz]`.

1. **Lazy evaluation** — the selection rule is defined once but resolved when the data is actually accessed, so it works regardless of which channels a file happens to contain
2. **Order preservation** — `channels([:Pz, :Cz])` returns data in the order you specified, not file order
3. **Composability** — you can combine rules and pass them to any function that accepts a selection keyword

## File Selection

In TOML configuration files, file selectors use strings or regex patterns:

```toml
### Basic patterns ###

# Single file
raw_data_files = "subject_01.bdf"

# All BDF files in the directory
raw_data_files = "\\.bdf$"

# All BrainVision files
raw_data_files = "\\.vhdr$"

# Files starting with "sub-" followed by digits
raw_data_files = "^sub-\\d+\\.bdf$"

### Filtering by participant range ###

# Only participants 1–9 (single digit)
raw_data_files = "sub-0[1-9]\\.bdf$"

# Two-digit participant IDs (01–99)
raw_data_files = "sub-\\d{2}\\.bdf$"

### Advanced patterns ###

# Exclude specific participants (03, 07, 12)
raw_data_files = "sub-(?!03|07|12)\\d+\\.bdf$"

# Match multiple naming conventions with alternation
raw_data_files = "(sub|pp)[-_]\\d+\\.bdf$"

# Files containing a specific task name anywhere
raw_data_files = "task-flanker.*\\.bdf$"
```

In Julia scripts, use `EegFun.get_files` to resolve a pattern against a directory:

```julia
# Pattern matching — same regex syntax as the TOML config
files = EegFun.get_files("./raw_data", "\\.bdf$")

# Only participants with two-digit IDs
files = EegFun.get_files("./raw_data", "sub-\\d{2}\\.bdf\$")

# Explicit list of files
files = EegFun.get_files("./raw_data", ["sub-01.bdf", "sub-05.bdf"])
```

Files are returned in natural sort order (e.g. `sub-2` before `sub-10`).

## Channel Selection

Used in functions like `subset`, `rereference!`, `repair_channels!`, and `detect_bad_epochs_automatic`.

| Predicate | Selects |
| --- | --- |
| `channels()` | All channels |
| `channels(:Cz)` | A single channel by name |
| `channels([:Fp1, :Fp2])` | A list of channels |
| `channels(1:10)` | Channels by index |
| `channels_not(:hEOG)` | Everything except the named channel(s) |
| `channels_not([:vEOG, :hEOG])` | Everything except a list of channels |

```julia
# Rereference only scalp channels (exclude EOG)
EegFun.rereference!(dat, :average, channel_selection = channels_not([:vEOG, :hEOG]))

# Run artifact detection on a subset of channels
EegFun.detect_bad_epochs_automatic(epochs, channel_selection = channels([:Cz, :Pz, :Fz]))
```

### Custom Channel Selection

Pass a function that receives a `Vector{Symbol}` of all available channel labels:

```julia
# Select all frontal channels (labels containing "F")
frontal = channels(x -> contains.(string.(x), "F"))
EegFun.subset(dat, channel_selection = frontal)
```

## Sample Selection

Used for filtering rows of continuous data by metadata columns (e.g. artifact flags). Common in `run_ica`, `correlation_matrix`, and `channel_summary`.

| Predicate | Selects |
| --- | --- |
| `samples()` | All samples |
| `samples(:is_epoch_interval)` | Samples where the flag column is `true` |
| `samples_not(:is_extreme)` | Samples where the flag is `false` |
| `samples_and([:is_epoch, :is_clean])` | Samples where **all** flags are `true` |
| `samples_or([:is_extreme, :is_muscle])` | Samples where **any** flag is `true` |
| `samples_and_not([:is_extreme, :is_step])` | Samples where **not all** flags are `true` |
| `samples_or_not([:is_extreme, :is_step])` | Samples where **none** of the flags are `true` |

```julia
# Run ICA only on clean, non-extreme data
ica = EegFun.run_ica(dat, sample_selection = samples_not(:is_extreme))

# Compute channel summary for the epoch interval only
EegFun.channel_summary(dat, sample_selection = samples(:epoch_interval))
```

### Custom Sample Selection

You can also pass a custom function where `x` refers to the data:

```julia
# Select samples between 100 and 1000
EegFun.subset(dat, sample_selection = samples(x -> 100 .<= x.sample .<= 1000))
```

## Time Selection

Used to define an interval in seconds for functions like `subset`, `erp_measurements`, and `plot_erp`.

| Predicate | Selects |
| --- | --- |
| `times()` | All time points |
| `times(0.3, 0.5)` | 300–500 ms |
| `times(-0.2, 0.0)` | Baseline period |

```julia
# Measure mean amplitude in the N400 interval
EegFun.erp_measurements(erps, interval_selection = times(0.3, 0.5))

# Subset data to the baseline period
EegFun.subset(dat, interval_selection = times(-0.2, 0.0))
```

## Epoch Selection

Used in `subset` and `reject_epochs` to filter by epoch index.

| Predicate | Selects |
| --- | --- |
| `epochs()` | All epochs |
| `epochs(1:20)` | Epochs 1–20 |
| `epochs(5)` | A single epoch |
| `epochs_not(1:5)` | All epochs except 1–5 |

```julia
# Keep only the first 50 epochs
EegFun.subset(epochs_data, epoch_selection = epochs(1:50))
```

### Custom Epoch Selection

Pass a function that receives a range of epoch indices:

```julia
# Select only even-numbered epochs
EegFun.subset(epochs_data, epoch_selection = epochs(x -> x .% 2 .== 0))
```

## Component Selection

Used with ICA functions to select or exclude specific independent components.

| Predicate | Selects |
| --- | --- |
| `components()` | All components |
| `components(1:3)` | Components 1–3 |
| `components(5)` | A single component |
| `components_not(1:3)` | All components except 1–3 |

```julia
# Remove identified artifact components
EegFun.remove_ica_components!(dat, ica, component_selection = components([1, 3, 7]))
```

## Summary

| Type | Select | Exclude | Combine |
| --- | --- | --- | --- |
| **Channels** | `channels(...)` | `channels_not(...)` | — |
| **Samples** | `samples(...)` | `samples_not(...)` | `samples_and(...)`, `samples_or(...)` |
| **Time** | `times(...)` | — | — |
| **Epochs** | `epochs(...)` | `epochs_not(...)` | — |
| **Components** | `components(...)` | `components_not(...)` | — |
| **Files** | regex / list | — | — |

## Next Steps

- **Using selections in preprocessing** — [Manual Preprocessing](manual-preprocessing.md) shows channel and sample selectors applied throughout the pipeline
- **Using selections with artifacts** — [Artifact Handling](artifact-handling.md) uses `channel_selection` and `sample_selection` in detection functions
- **File patterns in batch pipelines** — [Batch Processing](batch-processing.md) uses the regex file selector in TOML configuration
