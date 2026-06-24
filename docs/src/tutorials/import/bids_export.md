# BIDS Export

Export your EegFun data (raw files and/or preprocessed outputs) to a
[BIDS-compliant](https://bids.neuroimaging.io/index.html) directory structure.

The generated dataset can be validated using the
[BIDS Validator](https://bids-standard.github.io/bids-validator/).

## Minimal Example (Raw Data Only)

```julia
using EegFun

EegFun.export_bids(
    task = "posner",
    raw_dir = ".",
    output_dir = "./bids_dataset",
)
```

This discovers all `.bdf` files in `raw_dir`, auto-extracts subject IDs from the
filenames (`example1.bdf` → `sub-01`), and creates the BIDS structure with all required
sidecar files.

## Full Example

```julia
EegFun.export_bids(
    # Required
    task = "posner",
    raw_dir = ".",
    output_dir = "./bids_dataset",
    layout_file = "biosemi72.csv",

    # Dataset metadata
    dataset_name = "Visual Attention (Posner Cueing)",
    dataset_authors = ["Author A", "Author B"],
    dataset_license = "CC0",

    # Equipment
    manufacturer = "BioSemi",
    manufacturer_model = "ActiveTwo",
    cap_manufacturer = "BioSemi",
    cap_model = "headcap with 72 active electrodes",

    # Task
    task_description = "Posner cueing task with valid and invalid spatial cues",
    instructions = "Fixate on the central cross. Press the button as quickly as possible when the target appears.",

    # Recording context
    institution_name = "Your University",
    institution_address = "123 University Street, City, Country",
    institution_department = "Department of Psychology",
    eeg_ground = "CMS/DRL",
    eeg_placement_scheme = "10-10 subset with 72 electrodes",
    power_line_frequency = 50,

    overwrite = true,
)
```

## Including Preprocessed Data (Derivatives)

If you have already run the EegFun preprocessing pipeline, pass the
`derivatives_dir` argument to include the preprocessed JLD2 files:

```julia
EegFun.export_bids(
    task = "posner",
    raw_dir = ".",
    derivatives_dir = "./preprocessed_files",
    output_dir = "./bids_dataset",
    layout_file = "biosemi72.csv",
    dataset_name = "Visual Attention (Posner Cueing)",
    overwrite = true,
)
```

Pipeline output files are automatically renamed to BIDS naming conventions:

| Pipeline Suffix | BIDS Label |
|-----------------|------------|
| `_continuous_original` | `desc-continuousOriginal` |
| `_continuous_cleaned` | `desc-continuousCleaned` |
| `_epochs_original` | `desc-epochsOriginal` |
| `_epochs_cleaned` | `desc-epochsCleaned` |
| `_epochs_good` | `desc-epochsGood` |
| `_erps_original` | `desc-erpsOriginal` |
| `_erps_cleaned` | `desc-erpsCleaned` |
| `_erps_good` | `desc-erpsGood` |
| `_ica` | `desc-ica` |
| `_artifact_info` | `desc-artifactInfo` |

## Custom Subject ID Mapping

By default, subject IDs are auto-extracted from filenames using the first number
found (`example12.bdf` → `sub-12`). For custom mappings, use the `participant_map`
argument:

```julia
EegFun.export_bids(
    task = "posner",
    raw_dir = ".",
    participant_map = Dict(
        "Pract6" => "sub-01",
        "ExpA2"  => "sub-02",
        "test_file_99" => "sub-03",
    ),
)
```

## Generated Directory Structure

```text
bids_dataset/
├── dataset_description.json
├── participants.tsv
├── participants.json
├── sub-01/
│   └── eeg/
│       ├── sub-01_task-posner_eeg.bdf
│       ├── sub-01_task-posner_eeg.json
│       ├── sub-01_task-posner_channels.tsv
│       ├── sub-01_task-posner_events.tsv
│       ├── sub-01_task-posner_events.json
│       ├── sub-01_task-posner_electrodes.tsv
│       └── sub-01_task-posner_coordsystem.json
├── sub-02/ ...
└── derivatives/
    └── EegFun/
        ├── dataset_description.json
        └── sub-01/
            └── eeg/
                ├── sub-01_task-posner_desc-epochsGood_eeg.jld2
                ├── sub-01_task-posner_desc-erpsGood_eeg.jld2
                └── ...
```

## Auto-Computed Fields

The following fields are extracted automatically from your data — you don't
need to provide them:

| Field | Source |
|-------|--------|
| `SamplingFrequency` | From raw data |
| `EEGReference` | From `AnalysisInfo` |
| `EEGChannelCount` / `EOGChannelCount` / etc. | Auto-classified from channel names |
| `RecordingDuration` | Computed from data length |
| `RecordingType` | Always `"continuous"` for raw data |
| Channel names, types, units | From layout and channel labels |
| Event onsets and trigger values | From trigger column |
| Electrode positions (x, y, z) | From layout file |

## Validation

After exporting, validate your dataset using the
[BIDS Validator](https://bids-standard.github.io/bids-validator/).
Providing more of the optional metadata arguments will reduce the number
of validator warnings.

## API Reference

```@docs
EegFun.export_bids
```
