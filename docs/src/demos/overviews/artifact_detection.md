Artifact management is a critical step in EEG preprocessing to ensure that eye blinks, muscle activity, and technical noise do not bias your results. EegFun.jl provides a suite of tools for identifying and handling these artifacts.

### Key Detection Methods

- **Extreme Value Detection**: Identifies periods where the voltage exceeds reasonable physiological limits (e.g., ±100 μV).
- **Z-Score Rejection**: Uses statistical measures (variance, max, range, kurtosis) across epochs to identify outliers that deviate significantly from the mean signal.
- **EOG Onset Detection**: Specifically targets eye movement artifacts (blinks and saccades) using peak detection on EOG channels.

### Workflow

1. **Continuous Data**: Run extreme value detection or EOG onset detection.
2. **Epoch Data**: Apply statistical rejection using `detect_bad_epochs_automatic()`.
3. **Review**: Use the `plot_databrowser()` to manually verify flagged artifacts.

This demo demonstrates the core functions for each of these steps.
