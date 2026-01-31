
# Statistical Methods {#Statistical-Methods}

Understanding statistical analysis approaches in EegFun.jl.

## Overview {#Overview}

EegFun.jl provides several statistical testing approaches:
1. **Permutation Testing** - Non-parametric hypothesis testing
  
2. **Cluster-Based Statistics** - Control for multiple comparisons
  
3. **MVPA Decoding** - Multivariate pattern analysis
  
4. **RSA** - Representational similarity analysis
  

## Permutation Testing {#Permutation-Testing}

### Principle {#Principle}

Generate null distribution by randomly shuffling labels:
> 
> **TODO**: Add mathematical explanation and examples
> 


### When to Use {#When-to-Use}
- Non-normal data distributions
  
- Small sample sizes
  
- Robust inference
  

## Cluster-Based Permutation Tests {#Cluster-Based-Permutation-Tests}

### The Multiple Comparisons Problem {#The-Multiple-Comparisons-Problem}

Testing thousands of time points/channels inflates Type I error.

### Cluster-Based Solution {#Cluster-Based-Solution}
> 
> **TODO**: Explain clustering approach
> - Cluster formation
>   
> - Cluster-level statistics
>   
> - Permutation distribution
>   
> 


### Implementation {#Implementation}

```julia
# TODO: Add cluster test example
```


## MVPA Decoding {#MVPA-Decoding}

Classify experimental conditions from EEG patterns:
> 
> **TODO**: Add decoding workflow explanation
> - Cross-validation
>   
> - Classification algorithms
>   
> - Temporal generalization
>   
> 


## Representational Similarity Analysis (RSA) {#Representational-Similarity-Analysis-RSA}

Compare representational geometries:
> 
> **TODO**: Add RSA explanation
> - Dissimilarity matrices
>   
> - Model comparison
>   
> - Searchlight analysis
>   
> 


## Multiple Comparison Corrections {#Multiple-Comparison-Corrections}

### Bonferroni Correction {#Bonferroni-Correction}

Conservative but simple:
> 
> **TODO**: Add formulae and examples
> 


### FDR (False Discovery Rate) {#FDR-False-Discovery-Rate}

Less conservative:
> 
> **TODO**: Add explanation
> 


### Cluster-Based Correction {#Cluster-Based-Correction}

Already covered above.

## Effect Sizes {#Effect-Sizes}
> 
> **TODO**: Add effect size calculations
> - Cohen&#39;s d
>   
> - Partial eta-squared
>   
> - Interpretation guidelines
>   
> 


## See Also {#See-Also}
- [Statistical functions reference](../reference/statistics.md)
  
- Statistical testing tutorials (TODO)
  
