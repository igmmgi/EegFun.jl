# Global Field Power

This demo shows how to compute Global Field Power (GFP) and Global Dissimilarity from ERP data.

### What are GFP and Global Dissimilarity?

- **Global Field Power (GFP)** is the standard deviation across all channels at each time point. It provides a reference-independent measure of overall response strength — high GFP means strong, synchronised activity across the scalp.

- **Global Dissimilarity** measures the rate of topographic change over time. Peaks in dissimilarity indicate moments when the scalp distribution is changing rapidly, which may mark transitions between ERP components or brain states.

### Key Functions

| Function | Returns |
| --- | --- |
| `gfp` | DataFrame with `:time` and `:gfp` columns |
| `global_dissimilarity` | DataFrame with `:time` and `:dissimilarity` columns |
| `gfp_and_dissimilarity` | DataFrame with both `:gfp` and `:dissimilarity` |

### Options

- `channel_selection` — restrict to specific channels
- `normalize` — scale to 0–100% for cross-dataset comparison
- `condition_selection` — compute for specific conditions

## Workflow Summary

### GFP Calculation

- All channels or a subset; raw or normalised

### Global Dissimilarity

- Identify topographic transition points

### Combined Calculation

- Compute both metrics efficiently in one call


## Code Examples

::: details Show Code

```julia
# Demo: Global Field Power (GFP) and Global Dissimilarity
# Shows how to compute reference-independent measures of response strength
# and topographic change from ERP data.

using EegFun

#######################################################################
# LOAD ERP DATA
#######################################################################

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")


#######################################################################
# GLOBAL FIELD POWER (GFP)
#######################################################################

# GFP = standard deviation across all channels at each time point
# High GFP = strong, synchronized activity

# All channels
gfp_result = EegFun.gfp(dat)

# Specific channels
gfp_result = EegFun.gfp(dat, channel_selection = EegFun.channels([:Cz, :Pz, :Fz]))

# Normalized to 0-100% (for comparing across datasets)
gfp_result = EegFun.gfp(dat, normalize = true)


#######################################################################
# GLOBAL DISSIMILARITY
#######################################################################

# Global dissimilarity = rate of topographic change over time
# Peaks indicate transitions between brain states or ERP components

gd_result = EegFun.global_dissimilarity(dat)

# With normalization
gd_result = EegFun.global_dissimilarity(dat, normalize = true)


#######################################################################
# COMBINED CALCULATION
#######################################################################

# Calculate both GFP and global dissimilarity in one call
result = EegFun.gfp_and_dissimilarity(dat)

# Access the values (one result per condition)
result[1].time
result[1].gfp
result[1].dissimilarity


#######################################################################
# MULTIPLE CONDITIONS
#######################################################################

# GFP for all conditions
gfp_all = EegFun.gfp(dat)

# GFP for a specific condition
gfp_cond1 = EegFun.gfp(dat, condition_selection = EegFun.conditions([1]))
```

:::

## See Also

- [API Reference](../../reference/index.md)
