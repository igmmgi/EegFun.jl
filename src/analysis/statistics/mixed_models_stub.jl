"""
    fit_mass_lmm(args...; kwargs...)

Fit a mass-univariate linear mixed model across all channels and timepoints.

**Requires `MixedModels.jl` and `StatsModels.jl` to be loaded.**

To use this function, you must first import the required packages:
```julia
using EegFun
using MixedModels, StatsModels
```
"""
function fit_mass_lmm(args...; kwargs...)
    error("To use Mass-Univariate Mixed Models, you must first load the packages: `using MixedModels, StatsModels`")
end
