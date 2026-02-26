This demo shows how to combine epoch conditions, compute ERP difference waves, and create condition average waves.

### Key Functions

| Function | Purpose | Input |
| --- | --- | --- |
| `condition_combine` | Merge epoch conditions into new groups | Epoch data |
| `condition_difference` | Subtract one condition's ERP from another | ERP data |
| `condition_average` | Average multiple conditions into one ERP | ERP data |

### When to Use

- **`condition_combine`**: When you have many conditions and want to pool some together before averaging (e.g., merge "Go left" and "Go right" into a single "Go" condition)
- **`condition_difference`**: When you need difference waves for analysis (e.g., target minus standard, congruent minus incongruent)
- **`condition_average`**: When you want to average across conditions after computing ERPs (e.g., collapsing across irrelevant factors)

### Important Notes

- `condition_combine` works on epoch data (before averaging); `condition_difference` and `condition_average` work on ERP data (after averaging)
- All three also have batch versions that process all matching JLD2 files in a directory

## Workflow Summary

### Condition Combining

- Merge multiple epoch conditions into new groups (operates on epochs, before averaging)

### Condition Differences

- Create difference waves from ERP condition pairs (e.g., `[(1, 2), (3, 4)]`)

### Condition Averaging

- Average ERP conditions into combined waves (e.g., `[[1, 2], [3, 4]]` or `[[1, 2, 3, 4]]`)

### Typical Pipeline

- Preprocess → Combine conditions → Average → Difference/Average waves → Grand average
