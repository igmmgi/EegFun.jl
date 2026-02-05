This demo demonstrates artifact detection functionality, including identification, repair, and rejection strategies.

### What are EEG Artifacts?

EEG artifacts are unwanted signals in the recording that originate from sources other than brain activity. They can severely compromise data quality. Common artifact sources include:

**Physiological Artifacts:**

- **Eye movements (EOG)** - Horizontal and vertical eye movements create large amplitude deflections, particularly in frontal electrodes
- **Eye blinks** - Produce characteristic sharp, high-amplitude signals in frontal channels
- **Muscle activity (EMG)** - Jaw clenching, facial movements, and neck tension create high-frequency noise
- **Cardiac signals (ECG)** - Heart electrical activity can contaminate EEG

**External Artifacts:**

- **Electrode issues** - Poor contact, bridging, or movement can create extreme values
- **Movement artifacts** - Head or body movements produce low-frequency drifts and sharp transients
- **Environmental interference** - Line noise (50/60 Hz), electromagnetic interference from equipment

It is common in EEG experiments to lose trials due to artifacts, and the percentage of data remaining after artifact rejection is often reported in method sections. In some cases, entire datasets may be considered too noisy, resulting in too few trials remaining for a participant to be included in the final analysis. Again, such details are often (should be) reported in the methods section of published papers.

### Detection Strategies

Available approaches include:

1. **Threshold-based detection** - Identify samples/epochs exceeding absolute amplitude criteria (e.g., ±100 μV)
2. **Statistical detection** - Use z-scores to flag outliers based on variance, range, kurtosis, and other metrics
3. **EOG-specific detection** - Compute differential channels and detect characteristic eye movement patterns
4. **Visual inspection** - Interactive tools for manual review and refinement of automatic detection
5. **Epoch windows** - Mark time periods around triggers for focused artifact assessment

### Repair vs. Rejection

After detecting artifacts, there are two main options:

**Artifact Repair (Interpolation):**

- **When to use**: Bad channels in otherwise good epochs, isolated channel artifacts
- **Methods**:
  - `:neighbor_interpolation`
  - `:spherical_spline`
- **Advantage**: Retains maximum data, maintains epoch count
- **Limitation**: Cannot fix all artifacts (e.g., extreme eye movements affecting multiple channels)

**Artifact Rejection:**

- **When to use**: Entire epochs, or a large number of individual channels contaminated, movement artifacts, or when interpolation is inadequate
- **Advantage**: Removes contribution of noisy epoch to analysis
- **Limitation**: Reduces trial count, may bias remaining data if rejection is systematic

**Best/Recommended Practice**: Combine both approaches - repair isolated channel issues, reject severely contaminated epochs.

## IMPORTANT
While artifact detection and correction techniques (including ICA, interpolation, and rejection) are valuable tools, no computational method can rescue extremely poor data quality. **Garbage in, garbage out** applies to EEG analysis. The most effective approach is careful data collection: proper electrode preparation, clear participant instructions, comfortable seating, and monitoring during recording. Time invested in quality data acquisition will far exceed any benefit from post-hoc correction of noisy data.

## Workflow Summary

This demo demonstrates some features available in EegFun for artifact detection and correction:

### 1. Data Preparation

- Load raw BioSemi data and electrode layout
- Apply average reference and high-pass filtering (0.5 Hz)
- Visualize data in interactive databrowser

### 2. EOG Detection (Continuous Data)

- Compute vertical EOG channel (Fp1/Fp2 - IO1/IO2)
- Compute horizontal EOG channel (F9 - F10)
- Automatically detect eye movement onsets with threshold criteria

### 3. Extreme Value Detection (Continuous Data)

- Detect samples exceeding absolute voltage thresholds
- Apply detection globally or to specific channel/sample subsets
- Use predicate functions for flexible channel selection

### 4. Epoch Window Marking

- Mark time windows around experimental triggers
- Apply artifact detection specifically within epoch windows
- Visualize epoch boundaries and detected artifacts in databrowser

### 5. Bad Epoch Detection (Epoched Data)

- Extract epochs around experimental events
- **Automatic detection**: Combine absolute and z-score criteria across multiple metrics
  - Default metrics: variance, max, min, abs, range, kurtosis
  - Customizable thresholds and metric selection
- **Visual inspection**: Interactive grid view for manual verification
- **Hybrid approach**: Start with automatic detection, refine manually

### 6. Artifact Information

- Query unique rejection reasons, affected channels, and bad epochs
- Retrieve comprehensive artifact summaries per condition

### 7. Artifact Repair

- Interpolate bad channels using neighbor interpolation or spherical spline
- Visualize before/after comparison with repair plots
- Compare repair quality across different methods

### 8. Artifact Rejection

- Remove contaminated epochs entirely
- Visualize rejection impact on epoch count and data distribution
