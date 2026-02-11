# Filter

This demo shows how to apply temporal filters to EEG data to remove unwanted frequency components.

### Why Filter?

- **High-pass filter** removes slow drifts and DC offsets that obscure the signal of interest
- **Low-pass filter** removes high-frequency noise (muscle activity, line noise harmonics)
- Together they define the frequency band you want to analyse

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `highpass_filter!` | Remove frequencies below cutoff | 0.1 Hz for ERPs, 1 Hz for ICA |
| `lowpass_filter!` | Remove frequencies above cutoff | 30–40 Hz for ERPs |
| `create_highpass_filter` | Create a filter object for inspection | Visualize filter response |
| `create_lowpass_filter` | Create a filter object for inspection | Compare filter settings |
| `plot_filter` | Plot filter frequency response | Verify filter characteristics |

### IIR vs FIR

- **IIR** (default): Butterworth filter — fast, minimal parameters, good for most uses
- **FIR**: Windowed-sinc filter — linear phase, better when precise timing matters

## Workflow Summary

### 1. High-Pass Filtering

- Standard 0.1 Hz for ERP analysis
- Stronger 1 Hz for ICA preprocessing

### 2. Low-Pass Filtering

- 30 Hz for typical ERP analysis
- Anti-aliasing before resampling

### 3. Channel-Specific Filtering

- Apply filters to selected channels only

### 4. Filter Inspection

- Create filter objects and visualise frequency response
- Compare IIR vs FIR characteristics


## Code Examples

::: details Show Code

```julia
# Demo: Filtering
# Shows how to apply high-pass and low-pass filters to EEG data,
# compare IIR and FIR methods, and visualize filter characteristics.

using EegFun

#######################################################################
# 1. LOAD DATA
#######################################################################

dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)


#######################################################################
# 2. HIGH-PASS FILTER
#######################################################################

# Remove slow drifts with a 0.1 Hz high-pass (standard for ERPs)
EegFun.highpass_filter!(dat, 0.1)

# Stronger high-pass for ICA preprocessing (1 Hz)
dat_ica = copy(dat)
EegFun.highpass_filter!(dat_ica, 1.0)


#######################################################################
# 3. LOW-PASS FILTER
#######################################################################

# Remove high-frequency noise (30 Hz low-pass)
EegFun.lowpass_filter!(dat, 30.0)

# Anti-aliasing filter before downsampling (40 Hz)
EegFun.lowpass_filter!(dat, 40.0, order = 3)


#######################################################################
# 4. FILTER SPECIFIC CHANNELS
#######################################################################

# Apply filter to only certain channels
EegFun.lowpass_filter!(dat, 30.0, channel_selection = EegFun.channels([:Cz, :Pz, :Fz]))


#######################################################################
# 5. IIR vs FIR FILTERS
#######################################################################

# Default: IIR (Butterworth) — fast, good for most uses
EegFun.highpass_filter!(dat, 0.1, filter_method = "iir", order = 1)

# FIR (windowed sinc) — linear phase, better for time-domain analysis
EegFun.highpass_filter!(dat, 0.1, filter_method = "fir")


#######################################################################
# 6. VISUALIZE FILTER RESPONSE
#######################################################################

# Create a filter object for inspection
filt = EegFun.create_highpass_filter(0.1, 2048.0)
EegFun.plot_filter(filt)

filt_lp = EegFun.create_lowpass_filter(30.0, 2048.0)
EegFun.plot_filter(filt_lp)

# Compare IIR vs FIR
filt_iir = EegFun.create_highpass_filter(0.1, 2048.0, filter_method = "iir")
filt_fir = EegFun.create_highpass_filter(0.1, 2048.0, filter_method = "fir")
EegFun.plot_filter(filt_iir)
EegFun.plot_filter(filt_fir)
```

:::

## See Also

- [API Reference](../../reference/index.md)
