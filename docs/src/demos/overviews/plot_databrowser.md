This demo demonstrates the interactive databrowser for exploring and analyzing EEG data in real-time.

### What is the Databrowser?

The databrowser is an interactive visualization tool for navigating through EEG data:

- **Continuous data**: Scroll through raw recordings dynamically
- **Epoch data**: Navigate individual trials interactively
- **ICA components**: Inspect component time courses and topographies
- **Analysis settings**: Select regions, channels, and configurations visually

### Key Features

**Interactive Navigation**:

- **Scroll through time**: Mouse wheel or drag to pan
- **Zoom in/out**: Adjust time scale for detailed inspection
- **Channel selection**: Click channels to select for analysis
- **Region selection**: Drag to select time regions
- **Epoch browsing**: Navigate through individual trials

**Visual Overlays**:

- **Trigger markers**: See event timing
- **Epoch windows**: Visualize extraction boundaries
- **Extreme values**: Highlight artifact-contaminated samples
- **Channel repair**: Mark and repair problematic channels

**ICA Integration**:

```julia
plot_databrowser(dat, ica_result)
```

View ICA components alongside raw data with linked topographies.

### Analysis Settings

The databrowser returns analysis settings for further processing:

```julia
fig, ax, analysis_settings = plot_databrowser(dat)
dat_new = apply_analysis_settings(dat, analysis_settings)
```

Extract user-selected:

- Time regions
- Channel subsets  
- Sample masks
- Trigger configurations

### Keyboard Shortcuts

| Key | Action |
|-----|--------|
| **i** | Show help/info |
| **r** | Open channel repair menu |
| **c** | Clear current selections |
| **Arrow keys** | Navigate time windows |

### Custom Styling

Customize visual appearance:

```julia
plot_databrowser(dat,
    channel_line_width = 2,
    selection_color = (:red, 0.5)
)
```

### Working with Files

Load data directly from saved files:

```julia
plot_databrowser("dat.jld2")
plot_databrowser("epochs.jld2", "ica.jld2")
```

Useful for quick inspection without loading into memory.

### Common Workflows

**1. Initial Quality Control**:

- Load raw data
- Scroll through recording
- Identify noisy channels
- Mark extreme values
- Plan preprocessing strategy

**2. Artifact Identification**:

- Add `is_extreme_value!()` column
- Visualize artifact-contaminated regions
- Select clean segments for ICA

**3. Epoch Validation**:

- Mark epoch windows with `mark_epoch_windows!()`
- Verify trigger timing
- Check epoch boundaries
- Inspect single trials

**4. ICA Component Review**:

- Run ICA on clean data
- Browse components interactively
- Check topographies
- Identify artifact components

### Use Cases

**Interactive channel repair**:

- Select bad channels visually
- Choose repair method (spherical spline, neighbors)
- Apply and verify results immediately

**Spectral analysis region selection**:

- Identify interesting time windows
- Select for further time-frequency analysis
- Export settings for batch processing

**Preprocessing validation**:

- Verify filtering effects
- Check rereferencing results
- Validate artifact rejection

### Workflow Summary

This demo shows:

1. **Basic databrowser**: Simple visualization
2. **Analysis settings**: Capture and apply user selections
3. **Additional overlays**: Extreme values, epoch markers
4. **Custom styling**: Visual customization
5. **ICA integration**: Component inspection
6. **Epoch browsing**: Individual trial navigation
7. **File loading**: Direct file inspection

The databrowser is essential for interactive EEG exploration and quality control!
