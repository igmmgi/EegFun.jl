# Preprocessing Workflow

This demo walks through a complete manual preprocessing pipeline from raw BioSemi data to cleaned ERPs. It mirrors the steps and order used in `EegFun.preprocess_v1`, making it a useful reference for understanding what the automated pipeline does or for building a custom workflow.

### When to Use This

- **Learning** — understand each preprocessing step and why it is applied in a specific order
- **Custom workflows** — adapt the pipeline to your specific experimental design
- **Debugging** — step through the pipeline interactively to diagnose issues with your data

### Processing Phases

The workflow follows three phases:

### Phase 1: Setup and Initial Preprocessing

- Load raw data and configure electrode layout with coordinates and neighbours
- Mark epoch intervals for targeted analysis
- Rereference to average
- Apply initial filters (0.1 Hz high-pass for ERPs)
- Calculate and detect EOG (vEOG, hEOG)
- Extract "original" epochs for later comparison

### Phase 2: Artifact Detection and Cleaning (Continuous Level)

- Compute channel summary statistics and joint probability
- Detect extreme values at two thresholds (250 μV for ICA exclusion, 75 μV for artifact rejection)
- Identify and partition bad channels (EOG-related vs. non-EOG-related)
- Run ICA on a separate 1 Hz high-pass filtered copy, identify and remove artifact components
- Repair bad channels via neighbour interpolation
- Recalculate EOG channels after ICA and repair

### Phase 3: Epoch Extraction and Epoch-Level Processing

- Extract epochs from cleaned continuous data
- Baseline correction
- Two rounds of automatic epoch detection with channel repair between rounds
- Reject remaining bad epochs
- Average to ERPs and compare against originals

### Key Takeaway

The demo produces both "original" (lightly preprocessed) and "good" (fully cleaned) ERPs so you can visually assess the effect of the entire pipeline.


## Code Examples

::: details Show Code

```julia
# Demo: End-to-End Preprocessing Workflow
# Shows a complete manual preprocessing pipeline from raw data to ERPs.
# This mirrors the steps and order used in pipeline_v1.jl

using EegFun
using JLD2

const DEMO_OUTPUT = "./demos/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
#
# PHASE 1: BASIC SETUP AND INITIAL PREPROCESSING
#
# This phase prepares the data with basic preprocessing before
# artifact detection and ICA
#
#######################################################################

#######################################################################
# LOAD RAW DATA
#######################################################################

@info EegFun.section("Raw Data")

# Read raw data file
raw_data = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")

# Load electrode layout and configure all coordinate systems
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout)   # For 2D topography plots
EegFun.polar_to_cartesian_xyz!(layout)  # For 3D interpolation (spherical spline)
EegFun.get_neighbours_xy!(layout, 0.35) # Calculate neighbors (normalized units for this layout)

# Optional: Visualize layout and neighbors
# EegFun.plot_layout_2d(layout, neighbours = true)

# Create EegFun data structure
dat = EegFun.create_eegfun_data(raw_data, layout)

# Optional: save original continuous data
# jldsave(joinpath(DEMO_OUTPUT, "continuous_original.jld2"); data = dat)


#######################################################################
# MARK EPOCH WINDOWS
#######################################################################

@info EegFun.section("Marking epoch intervals")

# Define epoch conditions
epoch_cfg =
    [EegFun.EpochCondition(name = "Standard", trigger_sequences = [[1]]), EegFun.EpochCondition(name = "Target", trigger_sequences = [[2]])]

# Mark which samples belong to epoch intervals
# This creates a :epoch_interval column for targeted analysis
EegFun.mark_epoch_intervals!(dat, epoch_cfg, [-0.2, 1.5]) # -200 to 1.5 seconds 

# Optional: visualize triggers
# EegFun.plot_trigger_overview(dat)
# EegFun.plot_trigger_timing(dat)


#######################################################################
# REREFERENCE
#######################################################################

@info EegFun.section("Rereference")

# Rereference to average 
EegFun.rereference!(dat, :avg)


#######################################################################
# INITIAL FILTERS
#######################################################################

@info EegFun.section("Initial Filters")

# Apply high-pass filter (0.1 Hz for ERPs)
EegFun.highpass_filter!(dat, 0.1) # NB a 1Hz filter is applied to the data for ICA later

#######################################################################
# EOG CALCULATION AND DETECTION
#######################################################################

@info EegFun.section("EOG")

# Calculate vertical EOG (vEOG)
@info EegFun.subsection("Calculating EOG (vEOG/hEOG) channels")
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
)

# Calculate horizontal EOG (hEOG)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
)

# Detect EOG onsets (blinks and saccades)
@info EegFun.subsection("Detecting EOG (vEOG/hEOG) onsets")
EegFun.detect_eog_onsets!(dat, 50.0, :vEOG, :is_veog)
EegFun.detect_eog_onsets!(dat, 30.0, :hEOG, :is_heog)

# Calculate channel x EOG correlation matrices
@info EegFun.subsection("Channel x vEOG/hEOG Correlation Matrix")

# Create EOG configuration to match pipeline structure
eog_cfg = EegFun.EogConfig(
    vEOG_criterion = 50.0,
    hEOG_criterion = 30.0,
    vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
    hEOG_channels = [["F9"], ["F10"], ["hEOG"]],
)

# Calculate correlations across whole dataset
hEOG_vEOG_cm = EegFun.correlation_matrix_eog(dat, eog_cfg)
# Optional: view correlations
# println(hEOG_vEOG_cm)

# Calculate correlations within epoch intervals only
@info EegFun.subsection("Channel x vEOG/hEOG Correlation Matrix (epoch interval)")
hEOG_vEOG_cm_epoch = EegFun.correlation_matrix_eog(dat, eog_cfg, sample_selection = EegFun.samples(:epoch_interval))
# This is used later to identify channels with high EOG correlation

#######################################################################
# INITIAL EPOCH EXTRACTION (FOR COMPARISON)
#######################################################################

@info EegFun.section("Initial Epoch Extraction")

# Extract epochs from "lightly preprocessed" continuous data
# (filtered and rereferenced, but NO ICA or artifact removal yet)
epochs_original = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
erps_original = EegFun.average_epochs(epochs_original)

# Optional: save for comparison
# jldsave(joinpath(DEMO_OUTPUT, "epochs_original.jld2"); data = epochs_original)
# jldsave(joinpath(DEMO_OUTPUT, "erps_original.jld2"); data = erps_original)


#######################################################################
#
# PHASE 2: ARTIFACT DETECTION AND CLEANING (CONTINUOUS LEVEL)
#
# Now we identify and remove artifacts from the continuous data
# before final epoch extraction
#
#######################################################################

#######################################################################
# CHANNEL SUMMARY AND BAD CHANNEL DETECTION
#######################################################################

@info EegFun.section("Artifact Detection: Continuous Data")
@info EegFun.subsection("Channel Summary")

# Calculate channel statistics on whole dataset
summary_whole_dataset = EegFun.channel_summary(dat)

# Calculate channel statistics on epoch interval only
summary_epoch_interval = EegFun.channel_summary(dat, sample_selection = EegFun.samples(:epoch_interval))


#######################################################################
# EXTREME VALUE DETECTION
#######################################################################

@info EegFun.subsection("Artifact Detection (extreme values)")

# Detect extreme values (250 μV) - will exclude from ICA
@info "Detecting extreme values: 250.0 μV"
EegFun.is_extreme_value!(dat, 250.0, channel_out = :is_extreme_value_250)

# Detect artifact values (75 μV) - for epoch rejection
@info EegFun.subsection("Artifact Detection (criterion values)")
@info "Detecting artifact values: 75.0 μV"
EegFun.is_extreme_value!(dat, 75.0, channel_out = :is_artifact_value_75)


#######################################################################
# CHANNEL JOINT PROBABILITY
#######################################################################

@info EegFun.subsection("Bad Channel Detection using Channel Joint Probability + Z-Score Variance")

# Calculate joint probability on whole dataset
cjp_whole_dataset = EegFun.channel_joint_probability(dat)

# Calculate joint probability on epoch interval
cjp_epoch_interval = EegFun.channel_joint_probability(dat, sample_selection = EegFun.samples(:epoch_interval))

# Identify bad channels
@info EegFun.subsubsection("Bad Channels")
bad_channels_whole_dataset = EegFun.identify_bad_channels(summary_whole_dataset, cjp_whole_dataset)
bad_channels_epoch_interval = EegFun.identify_bad_channels(summary_epoch_interval, cjp_epoch_interval)

# Take intersection (more conservative)
bad_channels = intersect(bad_channels_whole_dataset, bad_channels_epoch_interval)

# Partition into EOG-related vs non-EOG-related
# EOG-related bad channels are better handled by ICA
bad_channels_non_eog_related, bad_channels_eog_related = EegFun.partition_channels_by_eog_correlation(
    bad_channels,
    hEOG_vEOG_cm_epoch;
    eog_channels = [:hEOG, :vEOG],
    threshold = 0.3,
    use_z = false,
)

@info "Bad channels (non-EOG related): $(length(bad_channels_non_eog_related)) - $bad_channels_non_eog_related"

# Create repair info for continuous-level repair
continuous_repair_info = nothing
if !isempty(bad_channels_non_eog_related)
    continuous_repair_info = EegFun.create_continuous_repair_info(:neighbor_interpolation)
    EegFun.channel_repairable!(continuous_repair_info, bad_channels_non_eog_related, layout)
end


#######################################################################
# ICA
#######################################################################

@info EegFun.section("ICA")

# Create a copy for ICA preprocessing
dat_ica = copy(dat)

# Remove bad channels before ICA (if any were identified as repairable)
if !isnothing(continuous_repair_info) && !isempty(continuous_repair_info.repaired)
    @info EegFun.subsection("Removing repairable bad channels for ICA")
    dat_ica = EegFun.subset(dat_ica, channel_selection = EegFun.channels_not(continuous_repair_info.repaired), include_extra = true)
    @info "Removed $(length(continuous_repair_info.repaired)) channels: $(continuous_repair_info.repaired)"
end

# Apply ICA-specific filters (1 Hz high-pass)
@info EegFun.subsection("ICA filters")
EegFun.highpass_filter!(dat_ica, 1.0)  # Stronger high-pass for ICA
EegFun.lowpass_filter!(dat_ica, 40.0)  # Optional: wider bandwidth

# Run ICA on clean sections only
@info EegFun.subsection("Running ICA")
ica_result = EegFun.run_ica(
    dat_ica,
    sample_selection = EegFun.samples_not(:is_extreme_value_250),
    percentage_of_data = 100,  # Use all clean data (or 20-50 for speed)
)

# Identify artifact components
@info EegFun.subsection("Component Identification")
component_artifacts, component_metrics = EegFun.identify_components(
    dat,  # Use original data for identification (important!)
    ica_result,
    sample_selection = EegFun.samples_not(:is_extreme_value_250),
)

# plot examples to visually inspect automatic component selections
EegFun.plot_topography(ica_result, component_selection = EegFun.components(component_artifacts.eog[:vEOG]))
EegFun.plot_topography(ica_result, component_selection = EegFun.components(component_artifacts.eog[:hEOG]))
EegFun.plot_topography(ica_result, component_selection = EegFun.components([1, 2, 3, 4]))
EegFun.plot_topography(ica_result, component_selection = EegFun.components(component_artifacts.line_noise))
EegFun.plot_ica_component_spectrum(
    dat,
    ica_result,
    component_selection = EegFun.components(component_artifacts.line_noise),
    sample_selection = EegFun.samples_not(:is_extreme_value_250),
)
EegFun.plot_topography(ica_result, component_selection = EegFun.components(component_artifacts.channel_noise))
EegFun.plot_ica_component_spectrum(dat_ica, ica_result, component_selection = EegFun.components(component_artifacts.channel_noise))
EegFun.plot_ica_component_activation(dat_ica, ica_result, component_selection = EegFun.components([1, 2, 3, 4]))
EegFun.plot_ica_component_activation(dat, ica_result)
EegFun.plot_ica_component_activation(dat, ica_result, artifacts = component_artifacts)

# Optional: view component metrics
component_metrics[:eog_metrics]
component_metrics[:line_noise_metrics]

# Remove artifact components
@info EegFun.subsection("Removing ICA components")
all_removed_components = EegFun.get_all_ica_components(component_artifacts)
@info "Removed $(length(all_removed_components)) ICA components: $component_artifacts"
EegFun.remove_ica_components!(dat, ica_result, component_selection = EegFun.components(all_removed_components))

# Optional: save ICA decomposition
# jldsave(joinpath(DEMO_OUTPUT, "ica.jld2"); data = ica_result)


#######################################################################
# REPAIR BAD CHANNELS (CONTINUOUS LEVEL)
#######################################################################

if !isnothing(continuous_repair_info)
    @info EegFun.section("Channel Repair")
    EegFun.repair_channels!(dat, continuous_repair_info; method = :neighbor_interpolation)
    # EegFun.repair_channels!(dat, continuous_repair_info; method = :spline_interpolation)
    @info continuous_repair_info
end


#######################################################################
# RECALCULATE EOG AFTER ICA AND REPAIR
#######################################################################

@info EegFun.section("EOG Recalculation")
@info EegFun.subsection("Recalculating EOG (vEOG/hEOG) channels after ICA and repair")
# EOG channels must be recalculated because underlying data changed
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
)


#######################################################################
# FINAL ARTIFACT DETECTION FOR EPOCHING
#######################################################################

@info EegFun.section("Detecting artifact values in continuous data")

# Re-detect artifacts after ICA and repair
EegFun.is_extreme_value!(dat, 75.0, channel_out = :is_artifact_value_75_final)

# Optional: save cleaned continuous data
# jldsave(joinpath(DEMO_OUTPUT, "continuous_cleaned.jld2"); data = dat)


#######################################################################
#
# PHASE 3: EPOCH EXTRACTION AND EPOCH-LEVEL PROCESSING
#
# Extract final epochs from cleaned data and perform epoch-level
# artifact rejection with repair
#
#######################################################################

#######################################################################
# EXTRACT EPOCHS (CLEANED)
#######################################################################

@info EegFun.section("Extracting cleaned epoched data")
# Extract epochs from fully cleaned continuous data
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))


#######################################################################
# BASELINE WHOLE EPOCHS
#######################################################################

@info EegFun.section("Baseline whole epochs")
# Apply baseline correction (using entire epoch by default)
EegFun.baseline!(epochs)


#######################################################################
# DETECT BAD EPOCHS (ROUND 1)
#######################################################################

@info EegFun.section("Automatic epoch detection")
# First round of bad epoch detection
rejection_info_step1 = EegFun.detect_bad_epochs_automatic(
    epochs;
    z_criterion = 0.0,    # No z-score (pipeline_v1 uses 0.0)
    abs_criterion = 75.0, # Absolute voltage threshold
    name = "rejection_step1",
)
# Identify which channels in bad epochs are repairable
EegFun.channel_repairable!(rejection_info_step1, epochs[1].layout)
@info rejection_info_step1


#######################################################################
# CHANNEL REPAIR PER EPOCH
#######################################################################

@info EegFun.section("Channel Repair per Epoch")
# Repair bad channels within individual epochs
# This rescues trials with isolated bad channels
EegFun.repair_artifacts!(epochs, rejection_info_step1)

# Optional: save epochs with repairs
# jldsave(joinpath(DEMO_OUTPUT, "epochs_cleaned.jld2"); data = epochs)


#######################################################################
# RE-DETECT ARTIFACTS AFTER REPAIR (ROUND 2)
#######################################################################

@info EegFun.subsection("Re-detecting artifacts after repair")
# Second round of detection to check repair effectiveness
rejection_info_step2 = EegFun.detect_bad_epochs_automatic(epochs; z_criterion = 0.0, abs_criterion = 75.0, name = "rejection_step2")
EegFun.channel_repairable!(rejection_info_step2, epochs[1].layout)
@info rejection_info_step2

# Optional: compare before and after repair
# rejection_comparison = EegFun.compare_rejections(rejection_info_step1, rejection_info_step2)


#######################################################################
# SAVE ARTIFACT INFO
#######################################################################

# Collect all artifact information
artifact_info = EegFun.ArtifactInfo(
    !isnothing(continuous_repair_info) ? [continuous_repair_info] : EegFun.ContinuousRepairInfo[],
    vcat(rejection_info_step1, rejection_info_step2),
    component_artifacts,
)
# jldsave(joinpath(DEMO_OUTPUT, "artifact_info.jld2"); data = artifact_info)


#######################################################################
# LOG EPOCH COUNTS
#######################################################################

# Compare original vs. final epoch counts
epoch_count_comparison =
    EegFun.log_epochs_table(epochs_original, epochs, title = "Epoch counts per condition (after repair and rejection):")


#######################################################################
# REJECT BAD EPOCHS
#######################################################################

@info EegFun.subsection("Rejecting bad epochs")
# Remove remaining bad epochs
epochs_good = EegFun.reject_epochs(epochs, rejection_info_step2)


#######################################################################
# SAVE GOOD EPOCH DATA
#######################################################################

# jldsave(joinpath(DEMO_OUTPUT, "epochs_good.jld2"); data = epochs_good)


#######################################################################
# AVERAGE TO ERPs
#######################################################################

erps_good = EegFun.average_epochs(epochs_good)
# jldsave(joinpath(DEMO_OUTPUT, "erps_good.jld2"); data = erps_good)


#######################################################################
# VISUALIZE COMPARISON
#######################################################################

@info EegFun.section("End of Processing")

# Compare original vs. cleaned ERPs
EegFun.plot_erp(erps_original, channel_selection = EegFun.channels([:Cz, :Pz]))
EegFun.plot_erp(erps_good, channel_selection = EegFun.channels([:Cz, :Pz]))


#######################################################################
# SUMMARY
#######################################################################

# This workflow demonstrates the processing steps from pipeline_v1:
#
# PHASE 1: BASIC SETUP AND INITIAL PREPROCESSING
#   - Load raw data and configure layout
#   - Mark epoch intervals
#   - Rereference (before filtering!)
#   - Apply initial filters (0.1 Hz high, 30 Hz low)
#   - Calculate EOG and detect onsets
#   - Extract "original" epochs for comparison
#
# PHASE 2: ARTIFACT DETECTION AND CLEANING (CONTINUOUS LEVEL)
#   - Channel summary statistics
#   - Extreme value detection (two thresholds: 250 μV and 75 μV)
#   - Channel joint probability and bad channel identification
#   - Run ICA (1 Hz high-pass, exclude bad channels and extreme values)
#   - Repair bad channels (continuous level)
#   - Recalculate EOG after ICA and repair
#   - Final artifact detection
#
# PHASE 3: EPOCH EXTRACTION AND EPOCH-LEVEL PROCESSING
#   - Extract epochs from cleaned continuous data
#   - Baseline correction
#   - Detect bad epochs (round 1)
#   - Repair channels per epoch
#   - Re-detect bad epochs (round 2)
#   - Save artifact info
#   - Log epoch count comparison
#   - Reject remaining bad epochs
#   - Save good epochs
#   - Average to ERPs
#   - Visualize comparison
#
# For automated batch processing of multiple participants,
# see the preprocess_v1() function and the batch-processing tutorial.
```

:::

## See Also

- [API Reference](../../reference/index.md)
