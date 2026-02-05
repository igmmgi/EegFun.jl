This demo demonstrates how to generate summary statistics across all channels for data quality assessment.

### What is Channel Summary?

Channel summary provides aggregate statistics across all EEG channels, offering a quick overview of data characteristics. This complements channel-specific metrics by revealing patterns across the entire montage.

**Statistical measures**:

- Mean amplitude per channel
- Variance/standard deviation
- Minimum and maximum values
- Sample counts

### Use Cases

**Data quality overview**:

- Quick assessment of recording quality across all channels
- Identify systematic patterns or issues across the montage
- Compare data quality before and after preprocessing

**Channel comparison**:

- Identify channels with systematically different properties
- Detect asymmetries or systematic biases in the recording
- Verify that preprocessing affected channels as expected

## Workflow Summary

This demo shows channel summary analysis:

### 1. Generate Initial Summary

- Load and preprocess data (average reference, high-pass filter)
- Calculate summary statistics across all channels
- Display results using formatted table

### 2. Summary with Sample Selection

- Mark extreme values for exclusion
- Recalculate summary excluding artifacts
- Compare clean vs. raw summary statistics
