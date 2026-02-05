## Overview

This demo demonstrates statistical analysis methods for ERP data.

### Statistical Testing

Compare experimental conditions using appropriate statistics:

**Parametric Tests:**
- **T-tests**: Compare two conditions
- **ANOVA**: Multiple conditions or factors
- **Assumptions**: Normality, equal variance

**Non-parametric Tests:**
- **Cluster-based permutation**: Control for multiple comparisons
- **No distribution assumptions**: More robust
- **Spatial-temporal clustering**: Accounts for dependencies

### Cluster-Based Permutation Tests

Powerful method for ERP analysis:
1. Compute test statistic at each time point/channel
2. Find clusters of contiguous significance
3. Permute condition labels and repeat
4. Compare observed clusters to permutation distribution

### Advantages

- Controls family-wise error rate
- Sensitive to spatiotemporal effects
- No stringent parametric assumptions
- Accounts for multiple comparisons elegantly

### Applications

- Condition comparisons in ERP studies
- Group differences
- Time-frequency power comparisons
