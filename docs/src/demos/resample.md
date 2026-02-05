# Resample

This demo shows how to change the sampling rate of EEG data through resampling.

This demo shows how to change the sampling rate of EEG data through resampling.

### What is Resampling?

Resampling changes the number of samples per second in your data:

- **Downsampling**: Reduce sampling rate (e.g., 2048 Hz → 512 Hz)
- **Upsampling**: Increase sampling rate (e.g., 250 Hz → 500 Hz)
- **Resample by factor**: Divide or multiply by an integer factor

### Why Resample?

**Downsampling benefits**:

- **Reduce file size**: Fewer samples = less disk space
- **Speed up processing**: Faster filtering, epoching, time-frequency analysis
- **Match dataset requirements**: Some analyses expect specific rates
- **Remove unnecessary detail**: Most ERP information is below 50 Hz

**Upsampling use cases**:

- **Match sampling rates**: Combine datasets with different rates
- **Synchronization**: Align with external recordings
- **Interpolation**: Create smoother visualizations

> **Note**: Upsampling cannot add information that wasn't in the original signal - it only interpolates between existing samples.

### The Nyquist Theorem

The sampling rate must be **at least 2× the highest frequency** in your signal:

- **250 Hz sampling** → Can represent frequencies up to 125 Hz
- **500 Hz sampling** → Can represent frequencies up to 250 Hz
- **1000 Hz sampling** → Can represent frequencies up to 500 Hz

**For ERP research**:
- Most ERP components are below 30-50 Hz
- **250-500 Hz** sampling is typically adequate
- Higher rates needed for high-frequency oscillations (gamma: 30-100 Hz)

### Anti-Aliasing

When downsampling, **always lowpass filter first** to prevent aliasing:

**Bad practice** (aliasing risk):
```julia
dat_new = resample(dat, 4)  # Downsample by 4× WITHOUT filtering
```

**Good practice** (safe downsampling):
```julia
# If original rate is 2048 Hz and downsampling to 512 Hz:
# Filter at ~200 Hz (below new Nyquist of 256 Hz)
lowpass_filter!(dat, 200)
dat_new = resample(dat, 4)  # Now safe to downsample
```

EegFun's `resample` function applies anti-aliasing automatically, but it's still good practice to filter first if you have specific frequency requirements.

### Downsampling Factors

The demo shows downsampling by factors of 2 and 4:

```julia
dat_new = resample(dat, 2)  # Divide rate by 2 (e.g., 2048 → 1024 Hz)
dat_new = resample(dat, 4)  # Divide rate by 4 (e.g., 2048 → 512 Hz)
```

**Common downsampling targets**:

| Original Rate | Factor | Target Rate | Use Case |
|---------------|--------|-------------|----------|
| 2048 Hz | 4 | 512 Hz | General EEG/ERP |
| 2048 Hz | 8 | 256 Hz | Standard ERP rate |
| 1024 Hz | 4 | 256 Hz | Standard ERP rate |
| 1024 Hz | 2 | 512 Hz | General EEG |
| 512 Hz | 2 | 256 Hz | Standard ERP rate |

### Trigger Preservation

The demo verifies that triggers are preserved during resampling:

```julia
trigger_count(dat)      # Original trigger count
trigger_count(dat_new)  # Should match after resampling
```

Trigger timing is automatically adjusted to match the new sampling rate.

### When to Resample

**Resample early** in your pipeline:
1. Load raw data
2. **Resample** (if needed)
3. Filter
4. Epoch
5. Analyze

This minimizes processing time for subsequent steps.

**Don't resample after epoching** unless necessary - it's more efficient to resample continuous data first.

### Workflow Summary

This demo demonstrates:

1. **Load continuous data** at original sampling rate
2. **Check current rate** with `sample_rate()`
3. **Downsample by factor** (2× and 4×)
4. **Verify new rate** matches expected value
5. **Verify triggers preserved** with `trigger_count()`

### Best Practices

**Choose appropriate target rate**:
- **250 Hz**: Minimum for standard ERP work
- **500 Hz**: Good balance for most EEG applications
- **1000+ Hz**: Needed for high-frequency analyses

**Filter before downsampling**:
- Prevents aliasing artifacts
- Use lowpass filter at ~80% of new Nyquist frequency

**Document the change**:
- Always note original and resampled rates in your analysis notes
- Important for interpretation and reproducibility


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.sample_rate(dat)   # current sample rate
EegFun.trigger_count(dat) # current triggers in file

dat_new = EegFun.resample(dat, 2) # downsample by a factor of 2
EegFun.sample_rate(dat_new)       # should = original ÷ 2
EegFun.trigger_count(dat_new)     # triggers should be preserved

dat_new = EegFun.resample(dat, 4) # downsample by a factor of 4
EegFun.sample_rate(dat_new)       # should = original ÷ 4
EegFun.trigger_count(dat_new)     # triggers should be preserved
```

:::

## See Also

- [API Reference](../reference/index.md)
