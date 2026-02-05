## Overview

This demo demonstrates multivariate pattern analysis (MVPA) for decoding experimental conditions from EEG data.

### What is MVPA Decoding?

MVPA uses machine learning to classify trials based on spatial patterns of activity:
- **Single-timepoint decoding**: Classify conditions from spatial patterns at each time point
- **Time-generalization**: Test if patterns learned at one time generalize to other times
- **Cross-validation**: Rigorous testing using train/test splits

### Workflow

1. **Extract epochs** around events of interest
2. **Prepare data** for classification (channel × time patterns)
3. **Train classifier** (e.g., SVM) on labeled trials
4. **Cross-validate** to estimate classification accuracy
5. **Statistical testing** (permutation tests for significance)

### Applications

- Decode stimulus categories from neural patterns
- Test information content at different processing stages
- Compare representational structure across conditions
