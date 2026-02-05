# Artifacts

This demo demonstrates artifact detection functionality, including identification, repair, and rejection strategies.

This demo demonstrates artifact detection functionality, including identification, repair, and rejection strategies.

### What are EEG Artifacts?

EEG artifacts are unwanted signals in the recording that originate from sources other than brain activity. They can severely compromise data quality. Common artifact sources include:

**Physiological Artifacts:**

- **Eye movements (EOG)** - Horizontal and vertical eye movements create large amplitude deflections, particularly in frontal electrodes
- **Eye blinks** - Produce characteristic sharp, high-amplitude signals in frontal channels
- **Muscle activity (EMG)** - Jaw clenching, facial movements, and neck tension create high-frequency noise
- **Cardiac signals (ECG)** - Heart electrical activity can contaminate EEG

**External Artifacts:**

- **Electrode issues** - Poor contact, bridging, or movement can create extreme values
- **Movement artifacts** - Head or body movements produce low-frequency drifts and sharp transients
- **Environmental interference** - Line noise (50/60 Hz), electromagnetic interference from equipment

It is common in EEG experiments to lose trials due to artifacts, and the percentage of data remaining after artifact rejection is often reported in method sections. In some cases, entire datasets may be considered too noisy, resulting in too few trials remaining for a participant to be included in the final analysis. Again, such details are often (should be) reported in the methods section of published papers.

### Detection Strategies

Available approaches include:

1. **Threshold-based detection** - Identify samples/epochs exceeding absolute amplitude criteria (e.g., ±100 μV)
2. **Statistical detection** - Use z-scores to flag outliers based on variance, range, kurtosis, and other metrics
3. **EOG-specific detection** - Compute differential channels and detect characteristic eye movement patterns
4. **Visual inspection** - Interactive tools for manual review and refinement of automatic detection
5. **Epoch windows** - Mark time periods around triggers for focused artifact assessment

### Repair vs. Rejection

After detecting artifacts, there are two main options:

**Artifact Repair (Interpolation):**

- **When to use**: Bad channels in otherwise good epochs, isolated channel artifacts
- **Methods**:
  - `:neighbor_interpolation`
  - `:spherical_spline`
- **Advantage**: Retains maximum data, maintains epoch count
- **Limitation**: Cannot fix all artifacts (e.g., extreme eye movements affecting multiple channels)

**Artifact Rejection:**

- **When to use**: Entire epochs, or a large number of individual channels contaminated, movement artifacts, or when interpolation is inadequate
- **Advantage**: Removes contribution of noisy epoch to analysis
- **Limitation**: Reduces trial count, may bias remaining data if rejection is systematic

**Best/Recommended Practice**: Combine both approaches - repair isolated channel issues, reject severely contaminated epochs.

## IMPORTANT
While artifact detection and correction techniques (including ICA, interpolation, and rejection) are valuable tools, no computational method can rescue extremely poor data quality. **Garbage in, garbage out** applies to EEG analysis. The most effective approach is careful data collection: proper electrode preparation, clear participant instructions, comfortable seating, and monitoring during recording. Time invested in quality data acquisition will far exceed any benefit from post-hoc correction of noisy data.

## Workflow Summary

This demo demonstrates some features available in EegFun for artifact detection and correction:

### 1. Data Preparation

- Load raw BioSemi data and electrode layout
- Apply average reference and high-pass filtering (0.5 Hz)
- Visualize data in interactive databrowser

### 2. EOG Detection (Continuous Data)

- Compute vertical EOG channel (Fp1/Fp2 - IO1/IO2)
- Compute horizontal EOG channel (F9 - F10)
- Automatically detect eye movement onsets with threshold criteria

### 3. Extreme Value Detection (Continuous Data)

- Detect samples exceeding absolute voltage thresholds
- Apply detection globally or to specific channel/sample subsets
- Use predicate functions for flexible channel selection

### 4. Epoch Window Marking

- Mark time windows around experimental triggers
- Apply artifact detection specifically within epoch windows
- Visualize epoch boundaries and detected artifacts in databrowser

### 5. Bad Epoch Detection (Epoched Data)

- Extract epochs around experimental events
- **Automatic detection**: Combine absolute and z-score criteria across multiple metrics
  - Default metrics: variance, max, min, abs, range, kurtosis
  - Customizable thresholds and metric selection
- **Visual inspection**: Interactive grid view for manual verification
- **Hybrid approach**: Start with automatic detection, refine manually

### 6. Artifact Information

- Query unique rejection reasons, affected channels, and bad epochs
- Retrieve comprehensive artifact summaries per condition

### 7. Artifact Repair

- Interpolate bad channels using neighbor interpolation or spherical spline
- Visualize before/after comparison with repair plots
- Compare repair quality across different methods

### 8. Artifact Rejection

- Remove contaminated epochs entirely
- Visualize rejection impact on epoch count and data distribution


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
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

# Plot data in data browser
EegFun.plot_databrowser(dat)
# NB. The databrowser is interactive, with mouse and keyboard combinations. 
# Press "i" for information about the interactions.

# Available functionality
# Region selection -> plot topography, plot spectrum
# Channel selectoin -> repair
# Rereference
# Filtering (NB. if high-pass already applied as above, not available again in browser)
# Additiional channels
# Channel selection
# Highlight extreme activity channel wise
# View trigger onsets


# Calculate a vertical and a horizontal EOG channel
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

# Automatically detect EOG onsets
EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 50, :hEOG, :is_hEOG)

# How many vEOG/hEOG onsets were detected?
EegFun.n_values(dat, :is_vEOG)
EegFun.n_values(dat, :is_hEOG)

# We should now be able to see newly calculated EOG signals in 
# the databrowser under "Extra Channels"
EegFun.plot_databrowser(dat)

# Look for extreme values (adds Bool column to dataframe :is_extreme_value_100)
EegFun.is_extreme_value!(dat, 100); # 100 mV criterion

# How many extreme values
EegFun.n_extreme_value(dat, 100)                    # across all electrodes
EegFun.n_extreme_value(dat, 100, mode = :separate)  # separately for each electrode

# Additional examples with channel selection predicate
EegFun.is_extreme_value!(dat, 100; channel_selection = EegFun.channels_not([:Fp1, :Fp2]));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> endswith.(string.(x), "z"));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")));
EegFun.is_extreme_value!(dat, 100; channel_selection = x -> .!(endswith.(string.(x), "z")), sample_selection = x -> x.sample .> 1000);


# Artifact Detection in Epoched Data
EegFun.trigger_count(dat) # our markers or triggers
EegFun.mark_epoch_windows!(dat, [1, 2, 3, 5, 6], [-0.2, 1.0]) # mark 200 ms and 1000 ms around trigger values

# Can use marked epoch window as a basis for automatic artifact detection
EegFun.is_extreme_value!(dat, 50; sample_selection = EegFun.samples(:epoch_window), channel_out = :is_extreme_value_epoch)

# We can now see the epoch windows and artifacts within them in the databrowser under "Epoch Windows"
EegFun.plot_databrowser(dat)


# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

# Detect bad epochs
# Default is combination of absolute and z-score based criteria with z-score based on all epochs using the 
# following metrics: [:variance, :max, :min, :abs, :range, :kurtosis]
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs)
bad_epochs[1]

# But we can change the metrics (here only an absolute criterion is used)
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

# Or here, using only z-score based criteria for variance
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2, z_measures = [:variance])

# We can inspect the automatic detection visually
EegFun.plot_artifact_detection(epochs[1], bad_epochs[1]) # first epoch
EegFun.plot_artifact_detection(epochs[2], bad_epochs[2]) # second epoch

# We can also perform visual artifact detection
bad_epochs_manual = EegFun.detect_bad_epochs_interactive(epochs[1], dims = (4, 4))


# Our a combination of both
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)
bad_epochs_manual_visual = EegFun.detect_bad_epochs_interactive(epochs[1], artifact_info = bad_epochs[1], dims = (4, 4))


# We can get information about the bad epochs that were automatically or visually detected
# Artifact Info
# First condition
EegFun.unique_rejections(bad_epochs[1])
EegFun.unique_channels(bad_epochs[1])
EegFun.unique_epochs(bad_epochs[1])
EegFun.get_rejected(bad_epochs[1])

# repair
bad_epochs = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 100)

epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :neighbor_interpolation)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :spherical_spline)
# epochs_repaired = EegFun.repair_artifacts(epochs, bad_epochs, method = :reject)

EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1])
EegFun.plot_artifact_repair(epochs[1], epochs_repaired[1], bad_epochs[1], ylim = (-100, 100))

# reject 
epochs_rejected = EegFun.reject_epochs(epochs, bad_epochs)
EegFun.plot_artifact_rejection(epochs[1], epochs_rejected[1], bad_epochs[1])
```

:::

## See Also

- [API Reference](../reference/index.md)
