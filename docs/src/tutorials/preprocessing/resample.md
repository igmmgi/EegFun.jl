# Resample

This demo shows how to change the sampling rate of EEG data through polyphase resampling.

## What is Resampling?

Resampling changes the number of samples per second (Hz) in your recording. `EegFun.resample` and `EegFun.resample!` resample continuous, epoched, or ERP data to a positive integer target sampling rate (e.g., 2048 Hz → 512 Hz or 256 Hz) using `DSP.resample` polyphase anti-aliasing filtering with mirror padding.

## Why Resample?

**Downsampling benefits**:

- **Reduce memory & disk footprint**: Fewer samples per channel
- **Speed up downstream steps**: Faster ICA, filtering, time-frequency decomposition, and MVPA decoding
- **Harmonize multi-site cohorts**: Bring recordings with different acquisition rates onto a common sampling grid

## The Nyquist Theorem

The target sampling rate must be **strictly greater than 2× the highest frequency** of interest:

- **256 Hz sampling** → Nyquist = 128 Hz (well-suited for standard ERP analyses < 40 Hz)
- **512 Hz sampling** → Nyquist = 256 Hz (well-suited for time-frequency and gamma oscillations up to ~100 Hz)
- **1000 Hz sampling** → Nyquist = 500 Hz

## Anti-Aliasing & Mirror Padding

`EegFun.resample` automatically constructs a polyphase anti-aliasing filter (`DSP.resample_filter`) and applies mirror padding at segment boundaries to suppress edge transients:

```julia
dat_512 = EegFun.resample(dat, 512) # Resample to 512 Hz
dat_256 = EegFun.resample(dat, 256) # Resample to 256 Hz
```

## Trigger & Metadata Preservation

Event triggers and extra channels (such as EOG or boolean artifact flags) are mapped to the new sample grid automatically:

```julia
EegFun.trigger_count(dat)      # Original trigger counts
EegFun.trigger_count(dat_256)  # Preserved after resampling
```

## Workflow Summary

1. **Load continuous data** at the acquisition sampling rate
2. **Inspect current rate** with `EegFun.sample_rate(dat)`
3. **Resample to target Hz** (`EegFun.resample(dat, 512)` or `EegFun.resample!(dat, 256)`)
4. **Verify triggers preserved** with `EegFun.trigger_count(dat_256)`


## Code Examples

::: details Show Code

```julia
# Demo: Resampling
# Shows how to resample EEG data to a target sampling rate (Hz).

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

EegFun.sample_rate(dat)      # current sample rate (Hz)
EegFun.trigger_count(dat)    # current triggers in file
EegFun.plot_databrowser(dat) # view current data

dat_512 = EegFun.resample(dat, 512) # resample to 512 Hz (polyphase anti-aliasing)
EegFun.sample_rate(dat_512)         # 512 Hz
EegFun.trigger_count(dat_512)       # triggers are preserved
EegFun.plot_databrowser(dat_512)    # view resampled data

dat_256 = EegFun.resample(dat, 256) # resample to 256 Hz
EegFun.sample_rate(dat_256)         # 256 Hz
EegFun.trigger_count(dat_256)       # triggers are preserved
EegFun.plot_databrowser(dat_256)    # view resampled data
```

:::

## See Also

- [API Reference](../../reference/index.md)
