This demo shows how to compute Global Field Power (GFP) and Global Dissimilarity from ERP data.

### What are GFP and Global Dissimilarity?

- **Global Field Power (GFP)** is the standard deviation across all channels at each time point. It provides a reference-independent measure of overall response strength — high GFP means strong, synchronised activity across the scalp.

- **Global Dissimilarity** measures the rate of topographic change over time. Peaks in dissimilarity indicate moments when the scalp distribution is changing rapidly, which may mark transitions between ERP components or brain states.

### Key Functions

| Function | Returns |
| --- | --- |
| `gfp` | DataFrame with `:time` and `:gfp` columns |
| `global_dissimilarity` | DataFrame with `:time` and `:dissimilarity` columns |
| `gfp_and_dissimilarity` | DataFrame with both `:gfp` and `:dissimilarity` |

### Options

- `channel_selection` — restrict to specific channels
- `normalize` — scale to 0–100% for cross-dataset comparison
- `condition_selection` — compute for specific conditions

## Workflow Summary

### GFP Calculation

- All channels or a subset; raw or normalised

### Global Dissimilarity

- Identify topographic transition points

### Combined Calculation

- Compute both metrics efficiently in one call
