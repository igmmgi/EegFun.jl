This demo demonstrates visualizing individual epoch waveforms with flexible layout options.

### What is Epoch Plotting?

Epoch plots display individual trial timecourses, showing:

- **Single-trial waveforms**: Raw epoch data (not averaged)
- **Time-locked activity**: Aligned to stimulus or response events
- **Trial variability**: Individual differences across repetitions
- **ERP overlay**: Optional averaged waveform on top

### Layout Options

**`:single`** (default):

- All channels overlaid on one plot
- Best for comparing few channels
- Clear view of channel differences

**`:grid`**:

- Each channel in separate subplot
- Auto-calculated grid dimensions
- Best for many channels
- Synchronized axes

**`:topo`**:

- Channels arranged by scalp position
- Preserves spatial relationships
- Intuitive interpretation
- Optional scale axis

### Use Cases

**Trial-level inspection**:

- View variability across trials
- Identify outlier epochs
- Assess signal quality

**Verify epoching**:

- Check time window selection
- Confirm event alignment
- Validate baseline correction

**Artifact detection**:

- Spot trial-specific artifacts
- Find contaminated epochs
- Guide rejection decisions

**Condition comparison**:

- Overlay multiple conditions
- Compare trial-level activity
- Assess consistency

### Customization Options

**Channel selection**:

- Select specific channels or regions
- Focus on electrodes of interest

**Visual styling**:

- Individual trial transparency
- Average line width
- Color schemes
- Axis limits

**Layout parameters**:

- Grid dimensions
- Spacing and gaps
- Scale axis (topo layout)

### Interactive Features

When `interactive=true` (default):

- **Control panel** (press 'c'): Toggle conditions, adjust baseline
- **Keyboard shortcuts**: Scaling, help
- **Linked axes**: Synchronized zoom/pan across subplots

## Workflow Summary

This demo shows epoch visualization workflows:

### 1. Load and Prepare Data

- Read raw data
- Apply preprocessing
- Extract epochs

### 2. Basic Plotting

- Plot single condition
- Select specific channels
- Compare conditions with overlay

### 3. Layout Options

- Single plot: All channels together
- Grid layout: Individual subplots
- Topo layout: Spatial arrangement

### 4. Channel Selection

- Focus on specific channels
- Combine with different layouts
- Maintain consistent scaling
