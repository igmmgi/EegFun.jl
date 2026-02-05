This demo demonstrates the complete ICA workflow for decomposing EEG into independent components, identifying artifacts, and removing them.

### What is ICA?

Independent Component Analysis (ICA) is a blind source separation technique that decomposes multi-channel EEG into maximally independent sources:

- **Decomposes mixed signals**: Separates brain activity from artifacts based on statistical independence
- **Identifies artifact sources**: Eye movements, blinks, heartbeat, muscle activity, line noise

### How ICA Works (Conceptually)

EEG is a mixture of sources:

Recorded EEG = Brain + Eye movements + Muscle + Heartbeat + Noise

ICA finds the unmixing matrix that separates these:

Components = Unmixing × Recorded EEG

Each component represents one independent source (e.g., one component might be pure eye blinks, another pure alpha rhythm).

### ICA Workflow

**1. Preprocessing**:

- Apply high-pass filter (≥1 Hz recommended)
- Detect and reject extreme artifacts or extreme sections of data
- Optionally remove bad channels
- Create EOG channels for identification

**2. Run ICA**:

- `run_ica()` on continuous or epoched data
- Choose algorithm: `:infomax` (default) or `:infomax_extended` (for sub-gaussian sources)
- Can use subset of data to speed up (20-50% often sufficient)
- Exclude extreme artifacts from decomposition

**3. Identify Artifact Components**:

- `identify_components()` - all artifact types automatically
- Or identify individually:
  - `identify_eog_components()` - eye movements and blinks
  - `identify_ecg_components()` - heartbeat
  - `identify_line_noise_components()` - 50/60 Hz interference
  - `identify_spatial_kurtosis_components()` - channel noise

**4. Visualize Components**:

- `plot_topography()` - spatial patterns (component topographies)
- `plot_ica_component_activation()` - time courses
- `plot_ica_component_spectrum()` - frequency content
- `plot_databrowser()` - interactive exploration

**5. Remove Artifacts**:

- `remove_ica_components()` - subtract artifact components
- `restore_ica_components()` - add back if needed (for validation)

### Component Identification

**EOG (Eye) Components**:

- **Spatial pattern**: Strong frontal (blinks) or side (saccades) loading
- **Time course**: Correlates with EOG channels
- **Frequency**: Low frequency (<5 Hz)

**ECG (Heartbeat) Components**:

- **Spatial pattern**: Localized or bilateral temporal
- **Time course**: Regular rhythm (~60-100 bpm)
- **Frequency**: Spectral peak at heart rate

**Line Noise Components**:

- **Spatial pattern**: Diffuse or specific to noisy channels
- **Frequency**: Strong peak at 50 or 60 Hz and harmonics

**Muscle Components**:

- **Spatial pattern**: Temporal or frontal regions
- **Frequency**: High frequency (>20 Hz)
- **Time course**: Brief bursts

**Brain Components**:

- **Alpha rhythm**: Posterior, ~10 Hz, eyes-closed modulation
- **Mu rhythm**: Central, ~10 Hz, motor-related
- **Task-related**: Modulated by experimental conditions

### Algorithms

**Infomax** (default):

- Fast and robust
- Optimized for super-gaussian sources (most artifacts)
- Good general-purpose choice

**Infomax Extended**:

- Handles both super- and sub-gaussian sources
- Better for mixed artifact types
- Slightly slower

### Best Practices

**Data requirements**:

- **Channels**: ≥30 recommended (more data = better separation
)
- **Duration**: Several minutes minimum
- **Preprocessing**: High-pass filter ≥1 Hz to remove slow drifts, remove extreme artifacts

**Component selection**:

- Use automated identification as starting point
- Visually inspect topographies and activations
- When in doubt, err on side of caution (don't remove)
- Document which components were removed

**Validation**:

- Compare before/after data quality
- Check that brain components preserved
- Verify artifact reduction in problematic channels/trials

### When to Use ICA

**Recommended**:

- Eye movement artifacts in visual tasks
- Muscle artifacts in speech/motor tasks
- Heartbeat artifacts in source localization
- Isolating specific brain rhythms

**Not recommended**:

- Very small datasets (<30 channels, <1 minute)
- Electrode bridging or poor contact (interpolate instead)
- When simple rejection suffices (e.g., few contaminated trials)

### Continuous vs. Epoched ICA

**Continuous data**:

- Use full dataset for decomposition

**Epoched data**:

- Can concatenate epochs for ICA

Both approaches work - choose based on your workflow and data characteristics.

## Workflow Summary

This demo shows the complete ICA pipeline:

### 1. Prepare Data

- Load and preprocess continuous data
- Create EOG channels
- Detect extreme values
- Apply high-pass filter

### 2. Run ICA

- Standard infomax on full dataset
- Extended infomax on subset (20%)
- Compare algorithms

### 3. Visualize Components

- Topographic maps of components
- Component activations over time
- Frequency spectra
- Interactive databrowser

### 4. Identify Artifacts

- Automated identification (all types)
- Individual identification methods
- Plot component features

### 5. Remove and Validate

- Remove identified artifact components
- Reconstruct to verify correctness
- Compare original vs. cleaned data

### 6. Apply to Epochs

- Run ICA on epoched data
- Same identification and removal workflow
- Validate removal quality
