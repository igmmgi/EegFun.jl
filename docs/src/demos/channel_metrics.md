# Channel Metrics

This demo demonstrates how to calculate and visualize channel quality metrics for identifying bad channels and detecting artifacts.

This demo demonstrates how to calculate and visualize channel quality metrics for identifying bad channels and detecting artifacts.

### What are Channel Metrics?

Channel quality metrics quantify various characteristics of the EEG signal at each electrode. These metrics help identify problematic channels that may need repair or rejection before analysis.

**Variance**:

- Measures overall signal power/amplitude variability at each channel
- High variance may indicate artifacts (muscle, movement)
- Very low variance suggests poor contact or disconnected electrode

**Kurtosis**:

- Measures the distribution shape (tailedness) of the signal
- High kurtosis indicates spiky signals with extreme values (artifacts)

**Correlation**:

- Average correlation between a channel and its spatial neighbors
- Low correlation suggests the channel is behaving differently from nearby electrodes
- Detects bad contacts, bridging, or electrode-specific noise

**Joint Probability**:

- Statistical outlier detection across channels simultaneously
- Identifies channels that deviate from the typical multi-channel pattern
- Combines information from variance and correlation metrics

### Use Cases

**Quality control**:

- Identify bad channels before epoching or averaging
- Screen data quality during or after acquisition
- Track bad channels across participants - if a specific electrode is consistently bad across multiple participants, this may indicate a broken electrode or faulty electrode set that needs replacement

## Workflow Summary

This demo shows channel quality assessment workflows:

### 1. Basic Channel Metrics

- Load and preprocess raw data (average reference, high-pass filter)
- Calculate channel joint probability metrics
- Identify channels with extreme values

### 2. EOG Correlation Analysis

- Compute vertical and horizontal EOG channels
- Detect EOG onsets automatically
- Calculate correlation between EEG channels and EOG
- Partition bad channels into EOG-related vs. non-EOG artifacts

### 3. Metric Interpretation

- Add z-score columns for standardized thresholds
- Identify channels that are artifact-related vs. bad contacts
- Guide decisions on repair (interpolation) vs. rejection


## Code Examples

::: details Show Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eegfun_data(dat, layout_file)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# extreme value bool column
EegFun.is_extreme_value!(dat, 100)

# channel joint probability
channel_joint_probability = EegFun.channel_joint_probability(dat)
channel_joint_probability = EegFun.channel_joint_probability(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))

# Calculate EOG signals
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

EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)


# Calculate correlations between all channels and EOG channels
cm = EegFun.correlation_matrix_dual_selection(
    dat,
    sample_selection = EegFun.samples(),  # All samples
    channel_selection1 = EegFun.channels(),  # All EEG channels
    channel_selection2 = EegFun.channels([:vEOG, :hEOG]),  # EOG channels
)
EegFun.add_zscore_columns!(cm)

bad_channels = [:Fp1, :AF3]
non_eog_related, eog_related = EegFun.partition_channels_by_eog_correlation(bad_channels, cm)
```

:::

## See Also

- [API Reference](../reference/index.md)
