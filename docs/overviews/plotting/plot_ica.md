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
