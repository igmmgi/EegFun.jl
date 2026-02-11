This demo shows how to automatically detect bad epochs using statistical criteria and inspect the results.

### Why Automatic Detection?

- **Z-score rejection** flags epochs that are statistical outliers across the dataset
- **Absolute threshold** catches epochs with extreme voltage values regardless of distribution
- Combining both gives robust detection with minimal false positives

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `detect_bad_epochs_automatic` | Flag bad epochs using z-score and/or absolute criteria | Primary detection step |
| `detect_bad_epochs_interactive` | Visual grid review with optional pre-filled detections | Manual verification |
| `unique_rejections` | Get all unique (channel, epoch) rejection pairs | Inspection |
| `unique_channels` | Get channels that contributed to rejections | Identify noisy channels |
| `unique_epochs` | Get indices of flagged epochs | Review specific trials |
| `get_rejected` | Get rejected epoch indices | Pass to repair/reject |
| `plot_artifact_detection` | Visualise detection results per condition | Quality check |

### Detection Criteria

- **`z_criterion`** (default: 3.0) — epochs with z-scores above this are rejected. Set to 0 to disable.
- **`abs_criterion`** (default: 100 μV) — epochs exceeding this absolute voltage are rejected. Set to 0 to disable.
- **`z_measures`** — which metrics to use for z-score: `[:variance, :max, :min, :abs, :range, :kurtosis]`

## Workflow Summary

### 1. Extract and Baseline Epochs

- Create epochs from continuous data before running detection

### 2. Run Automatic Detection

- Start with defaults, then tune `z_criterion` and `abs_criterion`
- Use `z_measures = [:variance]` for a simpler, variance-only approach

### 3. Inspect Results

- Use `unique_epochs`, `unique_channels`, and the printed summary to assess rejection rates
- Visualise with `plot_artifact_detection`

### 4. Optional Interactive Review

- Use `detect_bad_epochs_interactive` to visually confirm or adjust automatic detections
