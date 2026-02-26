This demo is a comprehensive walkthrough of artifact handling in EegFun: continuous-data artifact detection, EOG channels, epoched artifact detection, interactive review, artifact repair, and epoch rejection.

### Key Functions

| Function | Purpose |
| --- | --- |
| `channel_difference!` | Create derived channels (e.g., vEOG, hEOG) |
| `detect_eog_onsets!` | Detect EOG events in continuous data |
| `is_extreme_value!` | Flag samples exceeding a voltage threshold |
| `is_step_value!` | Flag sudden voltage jumps between consecutive samples |
| `detect_bad_epochs_automatic` | Flag bad epochs using z-score and/or absolute criteria |
| `detect_bad_epochs_interactive` | Visual grid review with optional pre-filled detections |
| `repair_artifacts` | Interpolate bad channels within flagged epochs |
| `reject_epochs` | Remove flagged epochs entirely |
| `plot_artifact_detection` | Visualise detection results per condition |
| `plot_artifact_repair` | Compare original vs repaired epochs |
| `plot_artifact_rejection` | Compare original vs rejected epochs |

## Workflow Summary

### EOG Channel Creation

- Derive vertical and horizontal EOG from frontal and periocular channels using `channel_difference!`
- Detect blink/saccade onsets with `detect_eog_onsets!` and mark affected intervals with `mark_epoch_intervals!`

### Continuous-Data Artifact Detection

- `is_extreme_value!` — marks samples exceeding an absolute threshold (e.g., ±100 μV)
- `is_step_value!` — marks sudden jumps that extreme value detection might miss (e.g., cable disconnections)
- Both support `channel_selection` and `sample_selection` for targeted detection
- Use `n_extreme_value` / `n_step_value` to count detections (combined or per-channel with `mode = :separate`)

### Epoch Artifact Detection

- **`z_criterion`** (default: 3.0) — epochs with z-scores above this are rejected. Set to 0 to disable.
- **`abs_criterion`** (default: 100 μV) — epochs exceeding this absolute voltage are rejected. Set to 0 to disable.
- **`z_measures`** — which metrics to use for z-score: `[:variance, :max, :min, :abs, :range, :kurtosis]`
- Use `channel_selection` for region-specific detection (e.g., frontal channels only)

### Inspect and Review

- Use `unique_epochs`, `unique_channels`, `unique_rejections`, and `get_rejected` to assess rejection rates
- Visualise with `plot_artifact_detection`
- Optionally use `detect_bad_epochs_interactive` for manual verification

### Repair and Rejection

- `repair_artifacts` interpolates bad channels within flagged epochs (`:neighbor_interpolation`, `:spherical_spline`)
- `reject_epochs` removes flagged epochs entirely
- Visualise results with `plot_artifact_repair` and `plot_artifact_rejection`
