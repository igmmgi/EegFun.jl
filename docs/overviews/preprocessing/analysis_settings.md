This demo shows how to save and replay preprocessing settings using `AnalysisSettings`.

### What are Analysis Settings?

`AnalysisSettings` stores a preprocessing recipe that can be applied to data in one step. This is useful for:

- Reproducing the same preprocessing across multiple datasets
- Storing settings from interactive exploration (e.g., databrowser GUI)
- Standardising processing across participants in a study

### Tracked Settings

An `AnalysisSettings` object can contain:

- **`hp_filter`**: High-pass filter cutoff frequency
- **`lp_filter`**: Low-pass filter cutoff frequency
- **`reference`**: Re-referencing scheme (e.g., `:avg`, `:Cz`)
- **`repaired_channels`**: Channels to interpolate
- **`repair_method`**: Interpolation method
- **`removed_ica_components`**: ICA components to remove
- **`selected_regions`**: ROI definitions

### Applying Settings

- **`apply_analysis_settings!(dat, settings)`**: Apply in-place (mutating)
- **`apply_analysis_settings(dat, settings)`**: Return a processed copy (non-mutating)
- **`apply_analysis_settings!(dat, ica, settings)`**: Apply with ICA component removal

### Analysis Info Tracking

Every EegFun data object automatically tracks its preprocessing history via `analysis_info`. This records what filters, references, and other operations have been applied.

## Workflow Summary

This demo covers:

- Inspecting analysis info on data objects
- Creating AnalysisSettings with filter and reference parameters
- Applying settings to data (both mutating and non-mutating)
