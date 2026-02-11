This demo shows how to combine and compute differences between experimental conditions at the batch level (processing all participants at once).

### Key Functions

| Function | Purpose | Input |
| --- | --- | --- |
| `condition_combine` | Merge epoch conditions into new groups | Epoch data |
| `condition_difference` | Subtract one condition's ERP from another | ERP data |

### When to Use

- **`condition_combine`**: When you have many conditions and want to pool some together before averaging (e.g., merge "Go left" and "Go right" into a single "Go" condition)
- **`condition_difference`**: When you need difference waves for analysis (e.g., target minus standard, congruent minus incongruent)

### Important Notes

- Both functions operate as batch processors — they find and process all matching JLD2 files in a directory
- `condition_combine` works on epoch data; `condition_difference` works on ERP data
- Output is saved to a new directory automatically

## Workflow Summary

### 1. Condition Combining

- Merge multiple epoch conditions into new groups across all participants

### 2. Condition Differences

- Create difference waves from ERP condition pairs across all participants

### 3. Typical Pipeline

- Preprocess → Combine conditions → Average → Difference waves → Grand average
