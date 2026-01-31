# ICA Workflow Tutorial

Independent Component Analysis (ICA) is a powerful technique for removing artifacts from EEG data.

## When to Use ICA

ICA is particularly effective for removing:
- Eye blinks and movements
- Cardiac artifacts (ECG)
- Muscle artifacts
- Channel noise

> **TODO**: Add guidance on when NOT to use ICA

## Running ICA

```julia
using EegFun

# Run ICA decomposition
# TODO: Add actual ICA function call
# ica_result = EegFun.run_ica(dat)
```

## Identifying Artifact Components

> **TODO**: Add component identification workflow
> - EOG correlation
> - ECG correlation  
> - Visual inspection of components
> - Spatial patterns

## Removing Artifacts

```julia
# TODO: Add component removal example
# clean_dat = EegFun.remove_ica_components(dat, ica_result, [1, 3, 5])
```

## Visualization

```julia
# TODO: Add ICA visualization examples
# - Component topographies
# - Component activations
# - Power spectra
```

## Best Practices

> **TODO**: Add best practices section
> - How many components to compute
> - Data requirements (continuous vs epoched)
> - When to run ICA in the pipeline

## Complete Example

```julia
# TODO: Complete end-to-end ICA workflow
```

## Next Steps

- Learn more about [ICA concepts](../explanations/ica.md)
- Explore [ICA plotting functions](../reference/ica.md)
