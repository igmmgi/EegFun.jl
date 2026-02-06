# Baseline

This demo demonstrates baseline correction methods for ERP data, including different baseline windows and correction approaches for continuous, epoch, and ERP data.

This demo demonstrates baseline correction methods for ERP data, including different baseline windows and correction approaches for continuous, epoch, and ERP data.

### What is Baseline Correction?

Baseline correction removes the pre-stimulus mean from each trial, ensuring that activity is measured relative to a neutral reference period. This is essential for ERP analysis because:

- **Removes slow drifts**: Eliminates DC offsets that vary across trials
- **Standardizes pre-stimulus activity**: Ensures comparable starting points across conditions
- **Enables amplitude interpretation**: Makes measurements interpretable relative to baseline activity

### Common Baseline Windows

**Pre-stimulus period** (most common):

- Standard: -200 to 0 ms before stimulus onset
- Captures typical pre-stimulus activity
- Assumes stable activity before stimulus

### Baseline Methods

**Mean correction** (default):

- Subtract the mean voltage during baseline period
- Standard approach in ERP research
- Appropriate when baseline is stable

## Best Practices

**Baseline window selection**:

- Choose based on experimental design and paradigm timing
- Standard: -200 to 0 ms for most ERP paradigms
- Document chosen window in methods section
- Keep consistent across conditions and participants

**Timing considerations**:

- Ensure baseline period is free from stimulus overlap
- Avoid including activity from previous trials
- Verify baseline period is stable (visually inspect)

**When to baseline**:

- Apply after filtering and before averaging
- Can apply to epoched or average ERPs

## IMPORTANT

Baseline correction assumes the pre-stimulus period represents neutral brain activity. This assumption may be violated in designs with:

- Very short inter-trial intervals
- Anticipatory activity (e.g., motor preparation)
- Sustained activity from previous trials

## Workflow Summary

This demo demonstrates baseline correction for different data types:

### 1. Baseline Continuous Data 

- Load raw BioSemi data
- Apply baseline to entire continuous recording
- Visualize DC offset removal in databrowser
- Baseline to specific timepoint (t=0)

### 2. Baseline Epoch Data

- Extract epochs around experimental events
- Apply baseline at different timepoints (t=0, t=0.5)
- Visualize effect on individual trials
- Compare baseline choices

### 3. Baseline ERP Data

- Average epochs into ERPs
- Baseline to single timepoint (t=0, t=0.5)
- Baseline to time window (-200 to 0 ms)
- Compare effects across baseline choices


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file);

# Plot 
EegFun.plot_databrowser(dat) # DC shift visible

# Baseline stuff and replot 
EegFun.baseline!(dat)
EegFun.plot_databrowser(dat) # now zero mean over all samples (reduced DC shift)

# NB. such DC shift is removed when applying a high-pass filter (e.g., 0.1 Hz)

# baseline to timepoint 
EegFun.baseline!(dat, (0, 0)); # timepoint = 0
EegFun.plot_databrowser(dat)

# EpochData
# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(epochs, (0.0, 0.0)) # baseline to t=0
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(epochs, (0.5, 0.5)) # baseline to t=0.5 (just for demo purposes!)
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(epochs, (-0.2, 0.0)) # baseline from -200 to 0 ms t=0 (common pre-event baseline)
EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))


# Create some ERP data
erps = EegFun.average_epochs(epochs)

EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(erps, (0.0, 0.0)) # baseline to t=0
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(erps, (0.5, 0.5)) # baseline to t=0.5 (just for demo purposes!)
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

EegFun.baseline!(erps, (-0.2, 0)) # baseline to t=-0.2 to 0.0
EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

# We can see the influence of baseline window interactively using the plot_erp_measurement_gui
EegFun.plot_erp_measurement_gui(erps)

```

:::

## See Also

- [API Reference](../reference/index.md)
