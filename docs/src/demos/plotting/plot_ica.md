# Plot ICA

This demo shows how to visualise ICA component topographies as scalp maps.

### Why Plot ICA Topographies?

- **Identify artefact components** — blink and eye movement components have characteristic frontal patterns
- **Verify decomposition** — well-separated components should have distinct, interpretable scalp maps
- **Select components** — focus on a subset of components for inspection or removal

### Key Functions

| Function | Purpose | Typical Use |
| --- | --- | --- |
| `plot_topography(ica)` | Plot all components | Full overview |
| `plot_topography(ica, component_selection=...)` | Plot selected components | Focused inspection |

### Key Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `component_selection` | `components()` | Which components to display |
| `dims` | auto | Grid size as `(rows, cols)` |
| `use_global_scale` | `false` | Share colour scale across all maps |
| `colormap` | `:RdBu` | Colour map for the topography |
| `colorbar_plot` | `true` | Show colorbars |

### What You'll Learn

1. Displaying all ICA component topographies in a grid
2. Selecting specific components by index or range
3. Adjusting grid layout and colour scaling
4. Using keyboard shortcuts to scale the colour range


## Code Examples

::: details Show Code

```julia
# Demo: Plotting ICA Component Topographies
# Shows how to visualise ICA component scalp maps as topographic plots,
# select specific components, and adjust display options.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

#######################################################################
# LOAD DATA AND RUN ICA
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# highpass at 1 Hz for ICA (critical for good decomposition)
EegFun.highpass_filter!(dat, 1.0)

# average reference
EegFun.rereference!(dat, :avg)

# create additional Bool colum for v. extreme values
EegFun.is_extreme_value!(dat, 200)

# run ICA (excludinng the v. extreme samples)
ica = EegFun.run_ica(dat, sample_selection = EegFun.samples_not(:is_extreme_value_200))

#######################################################################
# BASIC TOPOGRAPHY GRID
#######################################################################

# plot all components in a grid
EegFun.plot_topography(ica, label_plot = false, point_plot = false)

#######################################################################
# SELECT SPECIFIC COMPONENTS
#######################################################################

# plot only the first 10 components
EegFun.plot_topography(ica, component_selection = EegFun.components(1:10))

#######################################################################
# CUSTOM GRID SIZE
#######################################################################

# control rows × columns layout
EegFun.plot_topography(ica, component_selection = EegFun.components(1:12), dims = (4, 3))

#######################################################################
# DISPLAY OPTIONS
#######################################################################

# shared colour scale across all components
EegFun.plot_topography(ica, component_selection = EegFun.components(1:10), use_global_scale = true)

# adjust colour map
EegFun.plot_topography(ica, component_selection = EegFun.components(1:6), colormap = :RdBu)
```

:::

## See Also

- [API Reference](../../reference/index.md)
