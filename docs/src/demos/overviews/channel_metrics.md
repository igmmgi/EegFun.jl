This demo demonstrates how to calculate and visualize channel quality metrics for identifying bad channels and detecting artifacts.

### What are Channel Metrics?

Channel quality metrics quantify various characteristics of the EEG signal at each electrode. These metrics help identify problematic channels that may need repair or rejection before analysis.

**Variance**:

- Measures overall signal power/amplitude variability at each channel
- High variance may indicate artifacts (muscle, movement)
- Very low variance suggests poor contact or disconnected electrode

**Kurtosis**:

- Measures the distribution shape (tailedness) of the signal
- High kurtosis indicates spiky signals with extreme values (artifacts)

**Correlation**:

- Average correlation between a channel and its spatial neighbors
- Low correlation suggests the channel is behaving differently from nearby electrodes
- Detects bad contacts, bridging, or electrode-specific noise

**Joint Probability**:

- Statistical outlier detection across channels simultaneously
- Identifies channels that deviate from the typical multi-channel pattern
- Combines information from variance and correlation metrics

### Use Cases

**Quality control**:

- Identify bad channels before epoching or averaging
- Screen data quality during or after acquisition
- Track bad channels across participants - if a specific electrode is consistently bad across multiple participants, this may indicate a broken electrode or faulty electrode set that needs replacement

## Workflow Summary

This demo shows channel quality assessment workflows:

### 1. Basic Channel Metrics

- Load and preprocess raw data (average reference, high-pass filter)
- Calculate channel joint probability metrics
- Identify channels with extreme values

### 2. EOG Correlation Analysis

- Compute vertical and horizontal EOG channels
- Detect EOG onsets automatically
- Calculate correlation between EEG channels and EOG
- Partition bad channels into EOG-related vs. non-EOG artifacts

### 3. Metric Interpretation

- Add z-score columns for standardized thresholds
- Identify channels that are artifact-related vs. bad contacts
- Guide decisions on repair (interpolation) vs. rejection
