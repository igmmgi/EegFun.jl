# Automated Batch Processing

EegFun.jl provides automated pipelines that process an entire cohort of EEG recordings from raw data to cleaned epochs. You provide a TOML configuration file; the pipeline handles file discovery, preprocessing, artifact management, and cohort-level reporting.

## Running a Pipeline

```julia
using EegFun

# Run the v1 pipeline
EegFun.preprocess("config.toml")

```

Relative paths in the TOML file are resolved relative to the config file's directory, so you can keep the config alongside your analysis script.

## TOML Configuration

A pipeline config must specify **input files**, **output settings**, and **preprocessing parameters**. You only need to include values that differ from the defaults — EegFun merges your config with `src/config/default.toml`.

### Minimal Example

```toml
[files.input]
directory = "./raw_data"
raw_data_files = "\\.bdf"          # regex pattern matching raw files
layout_file = "biosemi64.csv"
epoch_condition_file = "epochs.toml"

[files.output]
directory = "./preprocessed_files"

[preprocess]
epoch_start = -0.2
epoch_end = 0.8
reference_channel = "avg"
```

### Key Configuration Sections

| Section | Controls |
| --- | --- |
| `[files.input]` | Raw data directory, file pattern, layout file, epoch condition file |
| `[files.output]` | Output directory, which intermediate files to save (epochs, ERPs, ICA) |
| `[preprocess]` | Epoch interval, reference channel |
| `[preprocess.eeg]` | Artifact thresholds (absolute μV and z-score) |
| `[preprocess.eog]` | EOG channel definitions and detection criteria |
| `[preprocess.filter.highpass]` | High-pass filter (cutoff, order, method) |
| `[preprocess.filter.lowpass]` | Low-pass filter (cutoff, order, method) |
| `[preprocess.filter.ica_highpass]` | Separate high-pass for the ICA copy (typically 1 Hz) |
| `[preprocess.filter.ica_lowpass]` | Optional low-pass for the ICA copy |
| `[preprocess.ica]` | Whether to apply ICA, percentage of data to use |
| `[preprocess.layout]` | Neighbour distance criterion for channel repair |

> [!TIP]
> See `src/config/default.toml` for all available options and their defaults. Copy it as a starting point and modify only the values you need to change.

### Epoch Condition File

The epoch condition file is a separate TOML that defines which trigger sequences form each experimental condition:

```toml
[epochs]

[[epochs.conditions]]
name = "standard"
trigger_sequences = [[1]]
reference_index = 1

[[epochs.conditions]]
name = "deviant"
trigger_sequences = [[2]]
reference_index = 1

# Multi-trigger sequence with timing constraint
[[epochs.conditions]]
name = "response_locked"
trigger_sequences = [[1, 100], [2, 100]]
reference_index = 2
timing_pairs = [[1, 2]]
min_interval = 0.1
max_interval = 1.5
```

## Processing Lifecycle

### Phase 1 — Setup

- Loads and validates the TOML configuration
- Reads the electrode layout and computes spatial neighbours
- Discovers raw data files matching the file pattern
- Parses epoch condition definitions

### Phase 2 — Continuous Processing (per file)

- **Rereferencing** — standardises electrode voltages (e.g. average reference)
- **Filtering** — high-pass to remove drift, optional low-pass
- **Artifact detection** — marks extreme values and identifies bad channels via channel joint probability and z-score variance
- **ICA** — optionally runs ICA on a separately filtered copy of the data, identifies artifact components (eye blinks, muscle, line noise), and removes them from the original
- **Channel repair** — interpolates bad scalp channels using neighbour weights

### Phase 3 — Epoching and Cleanup

- **Epoch extraction** — cuts the continuous data around triggers
- **Baseline correction** — subtracts the pre-stimulus mean
- **Per-epoch artifact detection** — flags epochs exceeding amplitude thresholds
- **Per-epoch repair** — interpolates channels that are only bad in specific epochs
- **Epoch rejection** — removes epochs that cannot be repaired

### Phase 4 — Cohort Consolidation

After all files are processed:

- **Epoch count summary** — how many trials survived per participant and condition
- **Channel reliability** — which electrodes were frequently repaired
- **ICA statistics** — average components removed per participant

## Error Isolation

If a single file fails (corrupt data, missing triggers, etc.), the pipeline logs the error and continues with the remaining files. You get:

- A **study-level log** summarising the entire run
- **Per-file logs** with every table, warning, and decision for each participant
- **Traceability** — output files record the pipeline version and config used

## Output Files

The pipeline saves intermediate and final outputs as JLD2 files (controlled by `[files.output]` flags):

| File suffix | Contents |
| --- | --- |
| `_continuous_original` | Raw continuous data after import |
| `_continuous_cleaned` | Continuous data after filtering, ICA, and repair |
| `_epochs_original` | Epochs before artifact handling |
| `_epochs_cleaned` | Epochs after channel repair |
| `_epochs_good` | Epochs after rejection (final clean data) |
| `_erps_original` / `_erps_cleaned` / `_erps_good` | Averaged ERPs at each stage |
| `_ica` | ICA decomposition results |
| `_artifact_info` | Artifact tracking (repaired channels, rejected epochs, ICA components) |

## Custom Pipelines

If the built-in pipelines don't match your workflow, generate a template with the standard boilerplate (config loading, logging, error handling) already wired up:

```julia
EegFun.generate_pipeline_template("my_pipeline.jl", "my_preprocess")
```

This creates a Julia file with a complete pipeline skeleton. Edit the processing steps to suit your experiment, then run:

```julia
include("my_pipeline.jl")
my_preprocess("config.toml")
```

## Post-Processing

After batch preprocessing, use `subset_bad_data` to exclude participants with too much data loss:

```julia
EegFun.subset_bad_data("preprocessed_files", 70.0)
```

This moves participants with less than 70% data retention to an "excluded" subdirectory.

## Philosophy

The pipeline follows a **minimal-intervention** approach to preprocessing, guided by the principle that EEG data is often better left alone than aggressively cleaned (Delorme, 2023). Rather than applying many sequential transformations — each of which risks distorting the signal — the pipeline focuses on:

- **Filtering conservatively** — a gentle high-pass to remove drift, with optional low-pass
- **Using ICA primarily for clear EOG artifacts** — eye blinks and saccades are the main targets for component removal, rather than attempting to classify and remove every possible source of noise
- **Repairing rather than rejecting** where possible — interpolating isolated bad channels preserves trial counts
- **Rejecting only clearly contaminated epochs** — using statistical thresholds rather than overly aggressive criteria

Several processing steps draw on the statistical thresholding approach described in FASTER (Nolan, Whelan, & Reilly, 2010), and the practical workflow has been influenced by Dudschig, Mackenzie, Strozyk, Kaup, & Leuthold (2016).

**References:**

- Delorme, A. (2023). EEG is better left alone. *Scientific Reports*, 13(1), 2372.
- Nolan, H., Whelan, R., & Reilly, R. B. (2010). FASTER: Fully automated statistical thresholding for EEG artifact rejection. *Journal of Neuroscience Methods*, 192(1), 152–162.
- Dudschig, C., Mackenzie, I. G., Strozyk, J., Kaup, B., & Leuthold, H. (2016). The sounds of sentences: Differentiating the influence of physical sound, sound imagery, and linguistically implied sounds on physical sound processing. *Cognitive, Affective, & Behavioral Neuroscience*, 16(5), 940–961.

## Next Steps

- **Understanding each step** — [Manual Preprocessing](manual-preprocessing.md) explains the rationale behind every pipeline stage
- **Artifact detection options** — [Artifact Handling](artifact-handling.md) for interactive review and fine-grained control
- **Epoch condition syntax** — [Epoch Selection](epoch-selection.md) for trigger patterns and TOML condition files
- **Selection predicates** — [Selection Patterns](selection-patterns.md) for channel, sample, and time filters used in config and scripts
