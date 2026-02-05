## Overview

This demo demonstrates the complete ICA workflow for artifact identification and removal.

### What is ICA?

Independent Component Analysis (ICA) decomposes EEG into maximally independent sources:
- Separates brain activity from artifacts (eye movements, muscle, heartbeat)
- Identifies components representing distinct sources
- Enables selective artifact removal without discarding trials

### ICA Workflow

1. **Decomposition**: Run ICA to separate independent components
2. **Component labeling**: Identify artifact vs. brain components
   - Eye movements: Strong frontal loading, correlates with EOG
   - Muscle artifacts: High-frequency, temporal regions
   - Brain components: Posterior alpha, central mu rhythm
3. **Artifact removal**: Remove bad components and reconstruct data
4. **Validation**: Verify artifact removal quality

### When to Use ICA

- **Eye movement correction**: Better than simple rejection for paradigms requiring visual fixation
- **Muscle artifact removal**: Particularly useful for speech/motor tasks
- **Source separation**: Isolate specific brain rhythms

### Considerations

- Requires sufficient data (recommended: >30 channels, several minutes of data)
- Component interpretation requires expertise
- Not appropriate for all artifacts (e.g., electrode bridging)
