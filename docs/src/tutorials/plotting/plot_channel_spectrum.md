# Plot Channel Spectrum

This demo demonstrates frequency domain analysis of EEG channels using power spectral density (PSD) visualization.

### What is Power Spectral Density?

Power Spectral Density (PSD) shows the distribution of signal power across different frequencies:

- **Decomposes time-domain signal** into frequency components
- **Reveals dominant rhythms** (alpha, theta, beta, etc.)
- **Detects artifacts** like line noise (50/60 Hz)
- **Quantifies frequency band power** for analysis

### EEG Frequency Bands

Standard frequency bands in EEG analysis:

| Band | Frequency Range |
|------|----------------|
| **Delta** | 0.5-4 Hz |
| **Theta** | 4-8 Hz |
| **Alpha** | 8-13 Hz |
| **Beta** | 13-30 Hz |
| **Gamma** | >30 Hz |

### Visualization Options

**Single Channel Spectrum**:

- View PSD for one channel
- Identify dominant frequency peaks
- Detect line noise or artifacts

**Multiple Channel Comparison**:

- Compare spectral profiles across electrodes
- Identify spatial patterns (e.g., posterior alpha)
- Assess differential effects across regions

**All Channels**:

- Overview of entire dataset
- Quickly spot problematic channels
- Identify global artifacts

### Common Use Cases

**1. Quality Control**:

- **Line noise detection**: Look for sharp peaks at 50/60 Hz
- **Artifact identification**: High power in unexpected frequency ranges
- **Channel malfunction**: Abnormal spectral profile compared to neighbors

### Interpretation Guidelines

**Posterior alpha peak**:

- Should be dominant in occipital channels
- Typically 8-13 Hz
- Suppressed with eyes open

**Line noise**:

- Sharp peak at exactly 50 or 60 Hz
- Often includes harmonics (100, 150, 200 Hz)
- Indicates need for notch filtering

### Technical Details

**FFT Parameters**:

The function uses Welch's method for PSD estimation:

- **Window size**: Automatic based on data length
- **Overlap**: 50% between windows
- **Window type**: Hamming window
- **Frequency resolution**: Determined by window size

**Best Practices**:

**Data requirements**:

- **Minimum duration**: A few seconds (longer = better frequency resolution)
- **Sample rate**: Should capture frequencies of interest (≥2× max frequency)
- **Artifact-free data**: Remove extreme values before spectral analysis

**Filtering considerations**:

- **High-pass**: Remove DC offset and slow drifts (≥0.1 Hz)
- **Low-pass**: Anti-aliasing already applied during acquisition
- **Notch filtering**: Consider before spectrum if you want to see line noise

**Selecting channels**:

- **Single channel**: For targeted analysis (e.g., Oz for alpha)
- **Regional selection**: Compare frontal vs posterior, left vs right
- **All channels**: Initial quality control and overview

### Workflow Summary

This demo shows a simple workflow:

1. **Load data** (raw continuous data)
2. **Create EEG structure** with channel layout
3. **Visualize spectra**:
   - All channels at once (overview)
   - Single channel (e.g., Fp1 for EOG artifacts)
   - Channel groups (e.g., frontal channels for beta analysis)


## Code Examples

::: details Show Code

```julia
# Demo: Channel Power Spectrum Plotting
# Shows how to visualise power spectra from both raw continuous data 
# and pre-computed SpectrumData, featuring interactive controls for 
# axis scaling, unit conversions, and frequency band overlays.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

#######################################################################
# LOAD AND PREPARE DATA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

EegFun.highpass_filter!(dat, 0.5)

#######################################################################
# WORKFLOW 1: COMPUTING ON THE FLY FROM RAW DATA
#######################################################################
# plot_channel_spectrum can compute the spectrum dynamically from raw data

# Plot all channels
EegFun.plot_channel_spectrum(dat)

# Plot a single channel with a specific title
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1]), title = "Fp1 Power Spectrum")

# Plot multiple frontal channels
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :F4]), title = "Frontal Channels Power Spectrum")


#######################################################################
# WORKFLOW 2: PLOTTING PRE-COMPUTED SPECTRUM DATA
#######################################################################
# If you have already run freq_spectrum, you can pass the result directly

# Compute spectrum first
spectrum = EegFun.freq_spectrum(dat)

# Plot all channels from pre-computed data
EegFun.plot_channel_spectrum(spectrum)

# Single channel from pre-computed data
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels(:Cz))

#######################################################################
# INTERACTIVE/INITIAL PARAMETER OVERRIDES
#######################################################################
# Note: You can toggle these visually in the GUI, but you can also 
# set them as starting defaults via keyword arguments.

# Start with log-scale axes
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels(:Cz), x_scale = :log10, y_scale = :log10)

# Start with Decibel (dB) units instead of linear power
EegFun.plot_channel_spectrum(spectrum, channel_selection = EegFun.channels([:Cz, :Pz]), unit = :dB, max_freq = 100.0)

# Custom styling for lines
EegFun.plot_channel_spectrum(
    spectrum,
    channel_selection = EegFun.channels([:Cz, :Oz]),
    linewidth = 3,
    line_alpha = 0.6,
    title = "Custom Styled Power Spectrum",
)
```

:::

## See Also

- [API Reference](../../reference/index.md)
