This demo shows how to perform common operations on time-frequency data: baseline correction, channel averaging and differencing, condition averaging and differencing, and grand averaging.

### Key Functions

| Function | Purpose |
| --- | --- |
| `tf_baseline!` / `tf_baseline` | Apply baseline correction to TF data |
| `channel_average!` / `channel_average` | Average channel groups (TF-specific) |
| `channel_difference!` / `channel_difference` | Subtract channel groups (TF-specific) |
| `condition_average` | Average conditions together (TF) |
| `condition_difference` | Subtract conditions (TF) |
| `grand_average` | Average across participants (TF) |

### Important Notes

- All TF operations update **both** `data_power` and `data_phase` DataFrames consistently
- The API mirrors the ERP operations but is dispatched on `TimeFreqData`
- `tf_baseline!` supports multiple methods: `:db`, `:absolute`, `:relative`, `:relchange`, `:normchange`, `:zscore`, `:percent`

## Workflow Summary

### Baseline Correction

- Apply after TF decomposition and before statistical analysis
- Decibel (`:db`) is the most common method for EEG
- Cannot be applied twice — re-compute TF data to change baseline

### Channel Operations

- `channel_average` creates ROI averages (e.g., midline, frontal cluster)
- `channel_difference` computes laterality or other contrasts
- Use `reduce = true` to drop original channels

### Condition Operations

- `condition_difference` subtracts one condition's power from another
- `condition_average` pools conditions together

### Grand Average

- Load multi-participant data with `read_all_data`
- `grand_average` groups by condition and averages power across participants
