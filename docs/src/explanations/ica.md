# Independent Component Analysis (ICA)

Understanding ICA for EEG artifact removal.

## What is ICA?

Independent Component Analysis is a blind source separation technique that decomposes EEG signals into statistically independent components.

## Why Use ICA for EEG?

ICA is particularly effective because:

1. **Artifact sources are independent**: Eye movements, heart beats, and brain activity have distinct patterns
2. **Non-invasive separation**: Recovers sources without prior knowledge
3. **Preserves brain activity**: Removes artifacts while keeping neural signals

## How ICA Works

> **TODO**: Add mathematical explanation
> - Independence assumption
> - Mixing model
> - Unmixing matrix

## Identifying Artifact Components

### EOG (Eye Movement) Components

Characteristics:
- Strong frontal topography
- Large amplitude
- Stereotyped temporal patterns

> **TODO**: Add example images

### ECG (Cardiac) Components

Characteristics:
- Regular temporal pattern (~1 Hz)
- Specific topography  
- High correlation with ECG channel

### Muscle Artifacts

Characteristics:
- High-frequency content
- Peripheral scalp distribution
- Bursting temporal pattern

## ICA Algorithms

> **TODO**: Discuss different algorithms
> - FastICA
> - Infomax
> - JADE
> - When to use each

## Best Practices

### Data Requirements

- **Length**: At least 100,000 samples (better: 500,000+)
- **Channels**: More data points than channels
- **Preprocessing**: High-pass filter first (removes slow drifts)

### Number of Components

> **TODO**: Guidance on component count selection

### When NOT to Use ICA

- Very short recordings
- Few channels (<20)
- Data with discontinuities

## Component Rejection vs. Artifact Correction

> **TODO**: Discuss trade-offs

## See Also

- [ICA workflow tutorial](../tutorials/ica-workflow.md)
- [ICA functions reference](../reference/ica.md)
