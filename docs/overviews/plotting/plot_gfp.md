This demo shows how to visualise Global Field Power (GFP) and Global Dissimilarity for ERP data.

### What is GFP?

- **GFP** is the spatial standard deviation across all electrodes at each time point
- It summarises the overall strength of the scalp field, independent of polarity
- **Global Dissimilarity** measures how much the scalp topography changes between time points

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_gfp(erp)` | Plot GFP for one condition | Quick inspection |
| `plot_gfp([erp1, erp2])` | Compare GFP across conditions | Condition effects |
| `plot_gfp(gfp_dataframe)` | Plot pre-computed GFP | Custom workflows |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `normalize` | `true` | Percentage (true) or raw μV (false) |
| `show_erp_traces` | `false` | Add overlaid ERP waveforms |
| `show_dissimilarity` | `false` | Add Global Dissimilarity panel |
| `channel_selection` | `channels()` | Subset of channels for GFP |

### What You'll Learn

1. Plotting GFP for single and multiple conditions
2. Showing ERP traces and dissimilarity in a multi-panel layout
3. Using pre-computed GFP results
4. Selecting specific channels for GFP calculation
