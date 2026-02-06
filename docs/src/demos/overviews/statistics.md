This demo demonstrates statistical approaches for comparing ERP conditions and controlling for multiple comparisons.

### Statistical Testing Approaches

**Parametric Tests:**

- **T-tests**: Two-condition comparisons
- **ANOVA**: Multiple conditions or factorial designs

**Non-parametric Tests:**

- **Cluster-based permutation**: Controls for multiple comparisons across space and time
- **Distribution-free**: No normality assumptions
- **Spatiotemporal sensitivity**: Leverages natural clustering in EEG data

### Cluster-Based Permutation Testing

This method addresses the multiple comparisons problem inherent in ERP analysis:

1. Compute test statistic at each time point/channel
2. Identify clusters of contiguous (spatial/temporal) significant effects
3. Permute condition labels and repeat many times
4. Compare observed cluster mass to permutation distribution

**Key advantages:**

- Controls family-wise error rate without being overly conservative (e.g., Bonferroni)
- Sensitive to spatially and temporally distributed effects

### Workflow Summary

This demo shows:

1. **Loading group data**: Multiple participants with condition labels
2. **Statistical comparison**: T-tests and cluster permutation
3. **Visualization**: Plotting significant time windows and topographies
