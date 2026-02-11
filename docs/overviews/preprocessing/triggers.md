This demo shows how to inspect trigger (event marker) data and search for trigger sequences in EEG recordings.

### Why Triggers?

Triggers mark experimental events (stimulus onsets, responses, conditions) in the continuous EEG recording. Before epoching, you need to verify that triggers are present, correctly coded, and occur at the expected frequencies.

### Key Functions

| Function | Purpose |
| --- | --- |
| `trigger_count` | Count occurrences of each trigger value |
| `search_sequence` | Find indices where trigger patterns occur |
| `plot_trigger_overview` | Visualise trigger distribution |
| `plot_trigger_timing` | Visualise trigger timing |

### Sequence Search Features

- **Single values**: find all onsets of trigger `1`
- **Multi-trigger sequences**: find `[1, 2]` (cue followed by target)
- **Wildcards**: `:any` matches any trigger between specific values
- **Ranges**: `1:5` matches triggers 1 through 5
- **Multiple sequences (OR logic)**: `[[1, 2], [3, 4]]` finds either sequence

## Workflow Summary

### 1. Trigger Counting

- View raw and cleaned trigger counts in a formatted table

### 2. Sequence Searching

- Find single triggers, multi-trigger sequences, and patterns with wildcards

### 3. Visualisation

- Plot trigger timing and distribution
