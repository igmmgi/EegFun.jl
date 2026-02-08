# Selection Patterns

`EegFun.jl` provides a unified system for selecting files, channels, samples, and conditions. This system uses **Functional Predicates**—functions that return rules—to make your analysis scripts both readable and powerful.

## 1. File Selection

File selection is driven by strings and regular expressions in your configuration or scripts.

### Basic Selection
```julia
# Single file
raw_data_files = "subject_01.bdf"

# List of files
raw_data_files = ["sub-01.bdf", "sub-05.bdf"]
```

### Regular Expressions (Regex)
Any selector that isn't a direct filename is treated as a Perl-style regex:
```julia
# All BDF files
raw_data_files = "\\.bdf$"

# Exclude specific IDs using negative lookarounds
# (Exclude 7, 14, 28; include any other ID from 1-39)
raw_data_files = "Flank_C_(?!7|14|28)([1-9]|[123][0-9])\\.bdf$"
```

## 2. Channel Selection

Functions like `subset`, `rereference`, and `repair_channels` use channel predicates.

- `channels()`: Select all available channels.
- `channels(:Cz)`: Select a single channel by name.
- `channels([:Fp1, :Fp2])`: Select a list of channels.
- `channels(1:10)`: Select channels by index.
- `channels_not(:hEOG)`: Select everything *except* the specified channels.

```julia
# Example: Rereferencing only scalp channels
rereference!(dat, :average, channel_selection = channels_not([:vEOG, :hEOG]))
```

## 3. Sample and Time Selection

You can select temporal windows using either row indices (Samples) or physical time (Seconds).

### Time Windows (`interval_selection`)
Use `times(start, stop)` to define a window in seconds. This is commonly used in `subset`, `plot_erp`, and `erp_measurements`.
```julia
# Select the 200ms baseline period
subset(dat, interval_selection = times(-0.2, 0.0))

# Measure mean amplitude in a 300-500ms window
erp_measurements(erps, interval_selection = times(0.3, 0.5))
```

### Sample Logic (`sample_selection`)
Use `samples()` predicates to filter by rows or metadata columns (like artifact flags). This is useful for excluding "bad" sections of data without removing the entire file.

- `samples()`: All samples.
- `samples(:is_artifact)`: Only samples where a specific flag is `true`.
- `samples_not(:is_extreme)`: Only samples where the flag is `false`.
- `samples_and([:is_epoch, :is_clean])`: Samples where *all* specified flags are `true`.
- `samples_or([:is_extreme, :is_muscle])`: Samples where *either* flag is `true`.

```julia
# Example: Running ICA only on clean data sections
ica = run_ica(dat, sample_selection = samples_not(:is_extreme))

# Example: Calculating correlation matrix only for the epoch window
cm = correlation_matrix_eog(dat, sample_selection = samples(:is_epoch_window))
```

## 4. Custom Selection with Lambdas

For maximum flexibility, you can pass anonymous functions (`x -> ...`) to the selection predicates.

### Custom Channel Selection
The function received by `channels()` is passed a `Vector{Symbol}` of all available channels.
```julia
# Select all channels starting with "F" (frontal)
frontal_chans = channels(x -> contains.(string.(x), "F"))
subset(dat, channel_selection = frontal_chans)
```

### Custom Sample Selection
The function received by `samples()` is passed the entire data `DataFrame`.
```julia
# Select samples where vEOG is within a specific range
custom_samples = samples(x -> 50 .> x.vEOG .> -50)
subset(dat, sample_selection = custom_samples)
```

### Custom Epoch Selection
The function received by `epochs()` is passed a range of indices `1:N`.
```julia
# Select only even-numbered epochs
even_epochs = epochs(x -> x .% 2 .== 0)
subset(epochs_data, epoch_selection = even_epochs)
```

## 5. Why use Predicates?

You might wonder why we use `channels([:Cz, :Pz])` instead of just `[:Cz, :Pz]`.

1.  **Lazy Evaluation**: The package doesn't need to know which channels exist *right now*; the rule is applied only when the data is actually accessed.
2.  **Order Preservation**: Predicates like `channels([...])` ensure that the resulting data follows the order you specified, not just the order in the file.
3.  **Composition**: You can easily combine rules (e.g., "Select all scalp channels AND they must be in the left hemisphere").

## 6. Summary Table

| Selection Type | Predicate Example | Usage Location |
| :--- | :--- | :--- |
| **Files** | `raw_data_files = "\\.bdf$"` | TOML Config |
| **Participants** | `participants(1:10)` | `read_all_data` |
| **Channels** | `channels_not([:EOG])` | `subset`, `rereference` |
| **Samples** | `samples_not(:is_bad)` | `run_ica`, `average_epochs` |
| **Time** | `times(0.1, 0.5)` | `subset`, `plot_erp` |
| **Epochs** | `epochs(1:20)` | `subset`, `reject_epochs` |
