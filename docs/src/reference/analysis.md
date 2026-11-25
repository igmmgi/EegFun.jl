```@meta
CollapsedDocStrings = true
```

# Analysis, I/O & Pipeline Functions

All public data I/O, layout, access/selection, preprocessing, ICA, ERP, RIDE, time-frequency, statistical, MVPA decoding, RSA, configuration, and pipeline functions in `EegFun.jl`, organized by subsystem.

## Index

```@index
Pages = ["analysis.md"]
```

## Data I/O & BIDS

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/io/read_data.jl",
    "src/io/eeglab.jl",
    "src/io/fieldtrip.jl",
    "src/io/xdf.jl",
    "src/io/bids_import.jl",
    "src/io/bids_export.jl",
]
```

## Electrode Layouts & Spatial Neighbors

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/layouts/layout.jl",
]
```

## Data Access, Selection & Utilities

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/types/core.jl",
    "src/utils/data/access.jl",
    "src/utils/data/from_raw.jl",
    "src/utils/data/rename.jl",
    "src/utils/data/subset.jl",
    "src/utils/data/subset_files.jl",
    "src/utils/batch.jl",
    "src/utils/create_data.jl",
    "src/utils/files.jl",
    "src/utils/gpu.jl",
    "src/utils/logging.jl",
    "src/utils/misc.jl",
    "src/utils/print.jl",
    "src/utils/viewer.jl",
]
```

## Preprocessing — Signal Operations

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/processing/signal/baseline.jl",
    "src/analysis/processing/signal/detrend.jl",
    "src/analysis/processing/signal/filter.jl",
    "src/analysis/processing/signal/cleanline.jl",
    "src/analysis/processing/signal/mirror.jl",
    "src/analysis/processing/signal/resample.jl",
]
```

## Preprocessing — Channel Operations

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/processing/channels/channel_average.jl",
    "src/analysis/processing/channels/channel_delete.jl",
    "src/analysis/processing/channels/channel_difference.jl",
    "src/analysis/processing/channels/channel_metrics.jl",
    "src/analysis/processing/channels/channel_repair.jl",
    "src/analysis/processing/channels/channel_summary.jl",
    "src/analysis/processing/channels/rereference.jl",
    "src/analysis/processing/channels/csd.jl",
]
```

## Preprocessing — Artifact Detection & Autoreject

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/processing/artifact_detection/eog_onsets.jl",
    "src/analysis/processing/artifact_detection/extreme.jl",
    "src/analysis/processing/artifact_detection/step.jl",
    "src/analysis/processing/artifact_detection/flatline.jl",
    "src/analysis/processing/artifact_detection/peak_to_peak.jl",
    "src/analysis/processing/artifact_detection/bridged.jl",
    "src/analysis/processing/artifact_detection/epoch_rejection_types.jl",
    "src/analysis/processing/artifact_detection/epoch_detection.jl",
    "src/analysis/processing/artifact_detection/repair.jl",
    "src/analysis/processing/artifact_detection/subset_bad_data.jl",
    "src/analysis/processing/artifact_detection/autoreject.jl",
]
```

## Preprocessing — Triggers, Epochs & Conditions

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/processing/epochs/triggers.jl",
    "src/analysis/processing/epochs/epochs.jl",
    "src/analysis/processing/epochs/apply.jl",
    "src/analysis/processing/conditions/condition_combine.jl",
    "src/analysis/processing/conditions/condition_difference.jl",
    "src/analysis/processing/conditions/condition_average.jl",
    "src/analysis/processing/conditions/grand_average.jl",
]
```

## Independent Component Analysis (ICA)

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/processing/ica/api.jl",
    "src/analysis/processing/ica/infomax.jl",
    "src/analysis/processing/ica/picard.jl",
    "src/analysis/processing/ica/amica.jl",
    "src/analysis/processing/ica/components.jl",
    "src/analysis/processing/ica/artifact_id.jl",
]
```

## ERP Analysis & RIDE

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/erps/erp_measurements.jl",
    "src/analysis/erps/gfp.jl",
    "src/analysis/erps/jackknife_average.jl",
    "src/analysis/erps/lrp.jl",
    "src/analysis/erps/realign.jl",
    "src/analysis/ride/ride.jl",
    "src/analysis/ride/decompose.jl",
    "src/analysis/ride/woody.jl",
    "src/analysis/ride/utils.jl",
]
```

## Time-Frequency & Spectral Analysis

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/time_frequency/freq_spectrum.jl",
    "src/analysis/time_frequency/tf_morlet.jl",
    "src/analysis/time_frequency/tf_multitaper.jl",
    "src/analysis/time_frequency/tf_stft.jl",
    "src/analysis/time_frequency/tf_itc.jl",
    "src/analysis/time_frequency/baseline.jl",
]
```

## Statistical Testing & Mixed Models

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/statistics/core.jl",
    "src/analysis/statistics/statistics.jl",
    "src/analysis/statistics/statistics_tf_helpers.jl",
    "src/analysis/statistics/permutations.jl",
    "src/analysis/statistics/lmm_stats.jl",
]
```

## MVPA Decoding & Representational Similarity Analysis (RSA)

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/analysis/decoding/preparation.jl",
    "src/analysis/decoding/liblinear.jl",
    "src/analysis/decoding/decoding.jl",
    "src/analysis/decoding/decoding_statistics.jl",
    "src/analysis/rsa/rsa.jl",
    "src/analysis/rsa/rsa_models.jl",
    "src/analysis/rsa/rsa_crossvalidation.jl",
    "src/analysis/rsa/rsa_noise_ceiling.jl",
]
```

## Configuration & Preprocessing Pipelines

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/config/config.jl",
    "src/pipelines/pipeline.jl",
    "src/pipelines/report.jl",
    "src/pipelines/utils/repair.jl",
    "src/pipelines/utils/repair_info.jl",
    "src/pipelines/utils/summaries.jl",
    "src/pipelines/utils/template.jl",
    "src/pipelines/utils/utils.jl",
]
```
