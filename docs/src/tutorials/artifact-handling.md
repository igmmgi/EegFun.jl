# Automated Artifact Handling

Artifacts are an unavoidable aspect of EEG research. `EegFun.jl` provides a structured workflow for identifying artifacts, repairing what can be salvaged, and rejecting what cannot.

## Continuous Data: Sample-Level Detection

Before epoching, you can flag artifact-contaminated samples in continuous data. This is useful for marking segments to exclude from downstream analysis (e.g., ICA training) or for data quality reporting.

### Extreme Values

`is_extreme_value!` flags samples where the absolute voltage exceeds a threshold. This catches amplifier saturation, large movement artifacts, and disconnected electrodes.

```julia
# Flag samples exceeding ±100 μV (adds column :is_extreme_value_100)
EegFun.is_extreme_value!(dat, 100)

# Check how many samples were flagged
EegFun.n_values(dat, :is_extreme_value_100)

# Only check specific channels
EegFun.is_extreme_value!(dat, 100, channel_selection = EegFun.channels([:Fp1, :Fp2]))
```

### Step Values

`is_step_value!` flags samples where the voltage jump between consecutive samples exceeds a threshold. This detects sudden discontinuities (cable pulls, electrode pops) that extreme value detection might miss if the signal returns quickly.

```julia
# Flag jumps > 50 μV between consecutive samples (adds column :is_step_value_50.0)
EegFun.is_step_value!(dat, 50.0)

# Check how many step artifacts were detected
EegFun.n_values(dat, :is_step_value_50.0)
```

> [!TIP]
> Both functions also have non-mutating versions (`is_extreme_value`, `is_step_value`) that return a boolean mask without modifying the data.

Once flagged, these artifact columns can be visualised in `EegFun.plot_databrowser` — select the relevant column from the **Extra Channels** dropdown menu to highlight flagged samples on the EEG traces.

## Continuous Data: Bad Channel Detection

Before epoching, it is worth checking whether any channels are consistently noisy or broken across the entire recording. In particular, removing or repairing excessively bad channels **before ICA** is strongly recommended.

### Channel Summary Statistics

`channel_summary` computes per-channel descriptive statistics — minimum, maximum, standard deviation, variance, range, and critically the **z-scored variance** (`zvar`). A channel whose variance is far from the group mean is likely problematic.

```julia
summary_df = EegFun.channel_summary(dat)
```

You can restrict the summary to samples within marked epoch intervals (or any other defined interval), which often gives a better picture of data quality. Periods outside the task — such as block breaks where participants may stretch or the experimenter may adjust electrodes — can inflate channel summary values and thus, influence channel quality judgements:

```julia
summary_epoch = EegFun.channel_summary(dat, sample_selection = EegFun.samples(:epoch_interval))
```

### Channel Joint Probability

`channel_joint_probability` offers a complementary metric. It computes how likely each channel's data is under a probability model — channels with improbable distributions are flagged.

```julia
cjp = EegFun.channel_joint_probability(dat)
```

> [!TIP]
> Blink artefacts can inflate variance and joint probability for frontal channels (e.g. Fp1, Fp2), causing them to be falsely identified as bad. If you have already detected EOG onsets with `detect_eog_signals!`, you can mark an interval around each onset and exclude those samples when computing channel metrics:
>
> ```julia
> # Mark a −100 to +300 ms interval around each vEOG onset
> EegFun.mark_epoch_intervals!(dat, :is_vEOG, [-0.1, 0.3])
>
> # Compute channel metrics excluding blink intervals and extreme values
> clean = EegFun.samples_or_not([:is_vEOG_interval, :is_extreme_value_100])
> summary_df = EegFun.channel_summary(dat, sample_selection = clean)
> cjp = EegFun.channel_joint_probability(dat, sample_selection = clean)
> ```

### Identifying Bad Channels

`identify_bad_channels` combines both measures. Channels exceeding the z-variance criterion **or** flagged by joint probability are labelled as bad:

```julia
bad = EegFun.identify_bad_channels(summary_df, cjp)

# More conservative criterion
bad = EegFun.identify_bad_channels(summary_df, cjp, zvar_criterion = 2.5)
```

### Partitioning by EOG Correlation

Some channels may look "bad" simply because they pick up a lot of eye-movement activity. These are better handled later by ICA rather than labeled as "bad" or interpolated. `partition_channels_by_eog_correlation` separates the two groups:

```julia
eog_cm = EegFun.correlation_matrix_eog(dat, eog_cfg)
EegFun.add_zscore_columns!(eog_cm)

non_eog_bad, eog_bad = EegFun.partition_channels_by_eog_correlation(
    bad, eog_cm;
    eog_channels = [:hEOG, :vEOG],
    threshold = 0.3,
)
```

### Repairing Bad Channels

Channels that are genuinely broken (not EOG-related) can be repaired by interpolation. Two methods are available:

**Neighbour interpolation** — weighted average of nearby channels. Requires pre-computed neighbour information:

```julia
# Compute neighbours first (only needed once)
EegFun.get_neighbours_xyz!(dat.layout, 0.5)

repair_info = EegFun.create_continuous_repair_info(:neighbor_interpolation)
EegFun.channel_repairable!(repair_info, non_eog_bad, dat.layout)
EegFun.repair_channels!(dat, repair_info; method = :neighbor_interpolation)
```

**Spherical spline interpolation** — uses all remaining good channels weighted by their spatial relationship on the scalp (Perrin et al., 1989). Does not require neighbours, but does require 3D coordinates in the layout:

```julia
EegFun.repair_channels!(dat, non_eog_bad; method = :spherical_spline)
```

## Epoch-Level: Artifact Detection

After epoching, each trial needs to be checked for artifacts. EegFun provides automatic detection, interactive visual review, and a combined workflow.

### Automatic Detection

`detect_bad_epochs_automatic` flags epochs based on two complementary criteria:

- **Z-score criterion** — channels whose variance, kurtosis, or range is far from the group mean across epochs
- **Absolute criterion** — channels where the voltage exceeds a hard threshold (e.g. ±100 μV)

An epoch is flagged if **any** channel exceeds **either** criterion.

```julia
# Default: z_criterion = 3.0, abs_criterion = 100 μV
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# More aggressive z-score threshold
artifacts = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 2.0)

# Absolute threshold only (disable z-score)
artifacts = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 80.0)
```

The return value is an `EpochRejectionInfo` object that records which channels and epochs were flagged, and why.

### Interactive Review

`detect_bad_epochs_interactive` opens a Makie grid view where you can visually inspect each epoch and toggle it as good or bad:

```julia
state = EegFun.detect_bad_epochs_interactive(epochs)
```

Controls:

- **Toggles** — check to reject, uncheck to keep
- **◀ Prev / Next ▶** — navigate pages
- **Show bad channels only** — filter to channels flagged by automatic detection (requires `artifact_info`)

### Combined Workflow

The most effective approach runs automatic detection first, then passes the results into the interactive viewer so flagged epochs are pre-marked and you can filter to only the problematic channels:

```julia
# Step 1: automatic detection
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# Step 2: visual review with pre-marked epochs
state = EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts)
```

## Epoch-Level: Channel Repair

If only one or two channels are noisy in an otherwise clean trial, you can repair them rather than discarding the whole epoch.

### Analyse Repairability

Before repairing with neighbour interpolation, use `channel_repairable!` to check if bad channels have enough clean neighbours:

```julia
EegFun.channel_repairable!(artifacts, epochs.layout)
```

### Apply Repair

```julia
# Method 1: Weighted neighbour interpolation (fast, requires neighbours)
EegFun.repair_artifacts!(epochs, artifacts, method = :neighbor_interpolation)

# Method 2: Spherical spline interpolation (higher quality, requires 3D coordinates)
EegFun.repair_artifacts!(epochs, artifacts, method = :spherical_spline)
```

## Epoch-Level: Rejection

After repairing what can be salvaged, remove the remaining bad epochs:

```julia
EegFun.reject_epochs!(epochs, artifacts)
```

## Participant-Level Subsetting

If a participant has lost too much data across a condition, it may be better to exclude them entirely:

```julia
# Exclude participants with < 70% data retention
EegFun.subset_bad_data("preprocessed_files", 70.0)
```

## Summary Workflow

A typical epoch-level artifact handling workflow:

```julia
# 1. Detect
artifacts = EegFun.detect_bad_epochs_automatic(epochs)

# 2. Review (optional)
state = EegFun.detect_bad_epochs_interactive(epochs, artifact_info = artifacts)

# 3. Repair (salvage single bad channels)
EegFun.channel_repairable!(artifacts, epochs.layout)
EegFun.repair_artifacts!(epochs, artifacts, method = :neighbor_interpolation)

# 4. Reject (remove non-repairable trials)
EegFun.reject_epochs!(epochs, artifacts)
```

## Next Steps

- **Before epoching** — [Manual Preprocessing](manual-preprocessing.md) covers continuous-data artifact detection and ICA in full context
- **Defining epochs** — [Epoch Selection](epoch-selection.md) for trigger-based epoch extraction which feeds into artifact handling
- **Automating rejection** — [Batch Processing](batch-processing.md) runs detect → repair → reject automatically across a cohort
- **Interactive exploration** — [Plot GUI](plot-gui.md) for visually reviewing flagged epochs and channel data
