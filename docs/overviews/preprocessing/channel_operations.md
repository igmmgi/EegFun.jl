This demo shows how to manipulate channels by averaging groups together, computing differences between them, and deleting unwanted channels.

### Key Functions

| Function | Purpose | Example |
| --- | --- | --- |
| `channel_difference!` | Subtract one channel group from another | EOG calculation |
| `calculate_eog_channels!` | Calculate both vEOG and hEOG from config | Convenience wrapper |
| `channel_average!` | Average channel groups into new columns | ROI means |
| `channel_delete!` | Remove channels from data and layout | Drop non-EEG channels |

### Common Use Cases

- **EOG channels**: `channel_difference!` calculates vEOG (Fp1/Fp2 minus IO1/IO2) and hEOG (F9 minus F10)
- **ROI averages**: `channel_average!` creates region-of-interest means (e.g., frontal, parietal)
- **Cleanup**: `channel_delete!` removes non-EEG channels (e.g., photodiode, reference)

## Workflow Summary

### 1. Channel Difference

- Calculate EOG channels from electrode pairs
- Create arbitrary difference channels

### 2. EOG Configuration

- Use `EogConfig` to calculate both EOG channels in one call

### 3. Channel Averaging

- Average channel groups with custom labels
- Reduce dataset to only averaged channels

### 4. Channel Deletion

- Remove single or multiple channels (mutating and non-mutating versions)
