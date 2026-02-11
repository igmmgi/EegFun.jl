This demo demonstrates visualizing electrode layout configurations in 2D and 3D space.

### What is Electrode Layout?

Electrode layout defines the spatial positions of EEG sensors on the scalp using:

- **2D coordinates** (x, y): Projected positions for topographic maps
- **3D coordinates** (x, y, z): Actual positions on spherical head model
- **Spherical coordinates** (incidence, azimuth): Angular positions

### Layout Visualization Functions

**plot_layout_2d**:

- Standard topographic view from above
- Shows electrode positions and labels
- Customizable head outline, markers, and labels
- Supports region-of-interest highlighting

**plot_layout_3d**:

- Full 3D spherical visualization
- Interactive rotation and zoom
- Shows true spatial arrangement

**Neighbor visualization**:

- Display spatial relationships between electrodes
- Helps verify interpolation neighborhoods
- Interactive hover to show connections

### Use Cases

**Verify electrode montage**:

- Confirm correct channel positions
- Check for mislabeled electrodes
- Validate custom layouts

**Publication figures**:

- Document electrode positions
- Highlight regions of interest
- Create publication-ready diagrams

**Analysis planning**:

- Understand spatial sampling density
- Plan interpolation strategies
- Define ROIs for analysis

### Customization Options

**Head outline**:

- Color, line width, radius

**Electrode markers**:

- Marker style, size, color
- Show/hide markers

**Labels**:

- Font size, color
- Position offsets
- Show/hide labels

**Regions of interest (ROIs)**:

- Highlight electrode groups
- Customizable borders and fills
- Multiple ROIs with different styles

### Neighbor Detection

Calculate spatial neighbors based on distance thresholds:

- **2D neighbors**: Using projected coordinates
- **3D neighbors**: Using spherical coordinates
- Export neighbor definitions for channel repair

## Workflow Summary

This demo shows electrode layout visualization:

### 1. Basic 2D Visualization

- Load electrode layout
- Convert to 2D coordinates
- Plot with default settings

### 2. Customize Appearance

- Adjust head outline (color, width, radius)
- Modify markers (style, size, color)
- Customize labels (size, color, offsets)

### 3. Add Regions of Interest

- Highlight electrode groups
- Customize ROI borders and fills
- Create multiple ROIs with different styles

### 4. Neighbor Visualization

- Calculate neighbor relationships
- Plot with neighbor connections
- Export neighbor definitions

### 5. Save Figures

- Export publication-ready figures
- Use CairoMakie for vector graphics
