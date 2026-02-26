This demo shows how to realign stimulus-locked epochs to a different time point, such as response time, for response-locked ERP analysis.

### When to Use Realignment

Stimulus-locked epochs are aligned to stimulus onset (t=0). Realignment shifts t=0 to a different event — typically the participant's response — so that you can study activity time-locked to that event instead.

Common use cases:

- **Response-locked ERPs** — study motor preparation relative to button press
- **Saccade-locked ERPs** — study activity relative to eye movement onset
- **Any event-locked analysis** — realign to any column in the epoch data

### How it Works

1. Each epoch's time vector is shifted so that the realignment value becomes t=0
2. All epochs are cropped to the common time interval that is valid across all trials
3. A uniform time vector is regenerated to ensure consistency

### Key Functions

| Function | Purpose |
| --- | --- |
| `realign!(epochs, :rt)` | Realign in place (mutating) |
| `realign(epochs, :rt)` | Return a realigned copy |
| `realign(file_pattern, :rt)` | Batch realign across participants |

## Workflow Summary

### Single-Participant Realignment

- Realign epochs to response time column

### Batch Realignment

- Process all participant files in a directory

### Typical Pipeline

- Extract stimulus-locked epochs → realign to RT → average → LRP → jackknife
