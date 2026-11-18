This demo shows how to realign stimulus-locked epochs to a different time point, such as response time, for response-locked ERP analysis.

### When to Use Realignment

Stimulus-locked epochs are aligned to stimulus onset (t=0). Realignment shifts t=0 to a different event — typically the participant's response — so that you can study activity time-locked to that event instead.

Common use cases:

- **Response-locked ERPs** — study motor preparation relative to button press
- **Saccade-locked ERPs** — study activity relative to eye movement trigger
- **Any trigger-locked analysis** — realign to any trigger code(s) present within each epoch

### How it Works

1. Each epoch's time vector is shifted so that the first matching realignment trigger becomes `t = 0`
2. Epochs that do not contain any of the requested realignment triggers are dropped
3. Remaining epochs are cropped to the shared common time interval across trials

### Key Functions

| Function | Purpose |
| --- | --- |
| `realign!(epochs, [201, 202])` | Realign in place to matching trigger(s) |
| `realign(epochs, [201, 202])` | Return a realigned copy |
| `realign(file_pattern, [201, 202])` | Batch realign across participants |

## Workflow Summary

### Single-Participant Realignment

- Realign epochs to response trigger codes (e.g., `[201, 202]`)

### Batch Realignment

- Process all participant files in a directory

### Typical Pipeline

- Extract stimulus-locked epochs → realign to RT → average → LRP → jackknife
