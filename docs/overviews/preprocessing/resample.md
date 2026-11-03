This demo shows how to change the sampling rate of EEG data through polyphase resampling.

### What is Resampling?

Resampling changes the number of samples per second (Hz) in your recording. `EegFun.resample` and `EegFun.resample!` resample continuous, epoched, or ERP data to a positive integer target sampling rate (e.g., 2048 Hz → 512 Hz or 256 Hz) using `DSP.resample` polyphase anti-aliasing filtering with mirror padding.

### Why Resample?

**Downsampling benefits**:

- **Reduce memory & disk footprint**: Fewer samples per channel
- **Speed up downstream steps**: Faster ICA, filtering, time-frequency decomposition, and MVPA decoding
- **Harmonize multi-site cohorts**: Bring recordings with different acquisition rates onto a common sampling grid

### The Nyquist Theorem

The target sampling rate must be **strictly greater than 2× the highest frequency** of interest:

- **256 Hz sampling** → Nyquist = 128 Hz (well-suited for standard ERP analyses < 40 Hz)
- **512 Hz sampling** → Nyquist = 256 Hz (well-suited for time-frequency and gamma oscillations up to ~100 Hz)
- **1000 Hz sampling** → Nyquist = 500 Hz

### Anti-Aliasing & Mirror Padding

`EegFun.resample` automatically constructs a polyphase anti-aliasing filter (`DSP.resample_filter`) and applies mirror padding at segment boundaries to suppress edge transients:

```julia
dat_512 = EegFun.resample(dat, 512) # Resample to 512 Hz
dat_256 = EegFun.resample(dat, 256) # Resample to 256 Hz
```

### Trigger & Metadata Preservation

Event triggers and extra channels (such as EOG or boolean artifact flags) are mapped to the new sample grid automatically:

```julia
EegFun.trigger_count(dat)      # Original trigger counts
EegFun.trigger_count(dat_256)  # Preserved after resampling
```

### Workflow Summary

1. **Load continuous data** at the acquisition sampling rate
2. **Inspect current rate** with `EegFun.sample_rate(dat)`
3. **Resample to target Hz** (`EegFun.resample(dat, 512)` or `EegFun.resample!(dat, 256)`)
4. **Verify triggers preserved** with `EegFun.trigger_count(dat_256)`
