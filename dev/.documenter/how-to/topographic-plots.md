
# How to Create Topographic Plots {#How-to-Create-Topographic-Plots}

Create publication-quality topographic maps of EEG scalp distributions.

## Basic Topographic Map {#Basic-Topographic-Map}

Plot a single time point:

```julia
using EegFun
using GLMakie

# TODO: Add basic topography plotting example
# fig = EegFun.plot_topography(erps, time=0.3)
```


## Topographic Map Series {#Topographic-Map-Series}

Plot multiple time points:

```julia
# TODO: Add time series example
# times = [0.1, 0.2, 0.3, 0.4, 0.5]
# fig = EegFun.plot_topography_series(erps, times=times)
```


## Customization {#Customization}

### Custom Time Windows {#Custom-Time-Windows}

```julia
# TODO: Average over time window
# fig = EegFun.plot_topography(erps, time=(0.3, 0.5))
```


### Color Scale {#Color-Scale}

```julia
# TODO: Add colormap customization
```


### Electrode Markers {#Electrode-Markers}

```julia
# TODO: Show/hide electrodes, customize markers
```


## Exporting for Publications {#Exporting-for-Publications}

Save in vector format for publications:

```julia
# TODO: Add save example
# save("topography.pdf", fig)
# save("topography.svg", fig)
```


## Common Patterns {#Common-Patterns}

### Pattern 1: P300 Component {#Pattern-1:-P300-Component}

```julia
# TODO: Show typical P300 topography (300-500ms)
```


### Pattern 2: N170 Component {#Pattern-2:-N170-Component}

```julia
# TODO: Show typical N170 topography
```


## Troubleshooting {#Troubleshooting}
> 
> **TODO**: Add common issues:
> - Missing electrodes
>   
> - Layout problems
>   
> - Color scale issues
>   
> 


## See Also {#See-Also}
- [Visualization architecture](../explanations/visualization.md)
  
- [Layout management](../reference/layouts.md)
  
