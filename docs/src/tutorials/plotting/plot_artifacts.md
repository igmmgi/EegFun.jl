# Plot Artifacts

This demo demonstrates methods for detecting and visualizing artifacts in EEG data.

## What are EEG Artifacts?

Artifacts are unwanted signals from non-neural sources:

**Physiological:**

- **Eye movements/blinks**: Large amplitude deflections in frontal channels
- **Muscle activity (EMG)**: High-frequency noise from jaw, face, or neck
- **Cardiac signals (ECG)**: Heart electrical activity contamination

**External:**

- **Electrode issues**: Poor contact, bridging, movement
- **Line noise**: 50/60 Hz environmental interference
- **Movement**: Low-frequency drifts and transients

## Detection Methods

**Threshold-based:**

- Absolute amplitude criteria (e.g., ±100 μV)
- Simple and interpretable
- Quick identification of extreme values

**Statistical:**

- Z-score outlier detection
- Variance, range, and kurtosis metrics
- Flags unusual distributions

**Visual inspection:**

- Grid-based epoch review
- Manual validation of automatic detection
- Refinement of rejection criteria

## Visualization Approaches

This demo shows multiple complementary views:

**Channel-wise:**

- Which electrodes are most affected
- Spatial patterns of artifact contamination
- EOG correlation for eye artifact identification

**Epoch-wise:**

- Trial-by-trial rejection statistics
- Temporal distribution of artifacts
- Condition-specific artifact rates

**Rejection summaries:**

- Percentage of data flagged per condition
- Channel quality metrics
- Decision support for repair vs. rejection

## Artifact Management Strategy

**Repair (Interpolation):**

- Use for isolated bad channels in good epochs
- Maintains trial count
- Cannot fix widespread contamination

**Rejection:**

- Use for severely contaminated epochs
- Removes entire trials
- Reduces but cleans the dataset

**Hybrid approach (recommended):**

1. Repair isolated channel issues
2. Reject epochs with widespread artifacts

## Workflow Summary

This demo shows:

1. **Loading data**: Continuous or epoched EEG
2. **Detection**: Applying threshold and statistical criteria
3. **Visualization**: Plotting detected artifacts across channels and epochs
4. **Assessment**: Evaluating artifact patterns to guide cleaning strategy
5. **Reporting**: Quantifying data attrition for methods documentation

## Critical Principle

> **Garbage in, garbage out**: No computational method can rescue poor quality data. The most effective strategy is careful data collection with proper electrode preparation and participant instruction. Artifact management is a secondary cleaning strategy, not a substitute for quality acquisition.


## Code Examples

::: details Show Code

```julia
# Demo: Artifact Detection Visualization
# Shows plots for visualizing artifact detection and rejection.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));
dat = EegFun.create_eegfun_data(dat)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# extract epochs
epoch_cfg = EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs, artifacts)

# Or with custom criteria
artifacts = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 100, z_criterion = 0)
EegFun.plot_artifact_detection(epochs, artifacts)

artifacts = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 0, z_criterion = 3)
EegFun.plot_artifact_detection(epochs, artifacts)

# extract multiple epochs
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))

# Artifacts
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Plot artifacts
EegFun.plot_artifact_detection(epochs, artifacts)
EegFun.plot_artifact_detection(epochs[1], artifacts[1])


# Interactive artifact detection
EegFun.detect_bad_epochs_interactive(epochs, dims = (4, 4))

# Or combine with automatic artifact detection
artifacts = EegFun.detect_bad_epochs_automatic(epochs)
EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts, dims = (4, 4))

# Or combine with automatic artifact detection
EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts, dims = (4, 4), colormap = :GnBu)




```

:::

## See Also

- [API Reference](../../reference/index.md)
