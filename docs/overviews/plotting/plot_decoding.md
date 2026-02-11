This demo shows how to visualise MVPA (multivariate pattern analysis) decoding results.

### Why Plot Decoding Results?

- **Temporal dynamics** — see when brain signals become discriminative
- **Error shading** — visualise variability across cross-validation folds or subjects
- **Significance markers** — overlay statistical test results on the accuracy curve

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_decoding(decoded)` | Plot single-subject accuracy | Quick inspection |
| `plot_decoding(decoded_list)` | Multi-subject subplot grid | Individual differences |
| `plot_decoding(decoded, stats)` | Accuracy + significance | Publication figure |

### What You'll Learn

1. Plotting decoding accuracy over time with error shading
2. Customising colours, line width and titles
3. Comparing subjects in a subplot grid
4. Overlaying significance markers from `test_against_chance`
