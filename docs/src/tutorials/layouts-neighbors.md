# Layouts and Neighbours

Spatial analysis, topographic plotting, and artifact repair in `EegFun.jl` all depend on knowing where electrodes are physically located. A **layout** stores electrode positions; **neighbours** are the set of nearby electrodes computed from those positions.

## Reading a Layout

A layout CSV contains electrode labels and their polar coordinates (incidence and azimuth angles, in degrees). EegFun ships with layouts for BioSemi and EasyCap caps:

```julia
layout = EegFun.read_layout("resources/layouts/biosemi/biosemi64.csv")

# Check the electrode labels
EegFun.channel_labels(layout)
```

> [!TIP]
> When you load data with `create_eegfun_data`, the layout is automatically attached to your data object. You can access it via `dat.layout`.

## Coordinate Systems

Layout files store positions in **polar coordinates** (incidence from the vertex, azimuth from the nose). Different operations require different Cartesian representations.

### 2D Coordinates (for Topographic Plots)

Topographic maps need flat x/y positions. The conversion projects electrodes onto a 2D plane, with the vertex (Cz) at the centre:

```julia
EegFun.polar_to_cartesian_xy!(layout)

# Check positions
EegFun.positions_2D(layout)
```

By default, radial distance is preserved — electrodes beyond the scalp edge (inc > 90°, e.g. eye channels from EEGLAB imports) will appear outside the head circle, matching EEGLAB's topoplot rendering. To normalise all electrodes to fit within a fixed radius instead:

```julia
EegFun.polar_to_cartesian_xy!(layout, preserve_radial_distance = false)
```

### 3D Coordinates (for Spherical Spline Repair)

Spherical spline interpolation and source localisation require 3D positions on the unit sphere:

```julia
EegFun.polar_to_cartesian_xyz!(layout)

# Check positions
EegFun.positions_3D(layout)
```

> [!NOTE]
> Both coordinate conversions modify the layout in-place and clear any previously calculated neighbour information, since the distance metric changes.

## Visualising the Layout

Use `plot_layout_2d` or `plot_layout_3d` to verify that electrode positions look correct before proceeding with analysis:

```julia
# 2D view with head outline
EegFun.plot_layout_2d(layout)

# 3D view
EegFun.plot_layout_3d(layout)
```

## Computing Neighbours

Many algorithms — neighbour interpolation for channel repair, spatial clustering, correlation-based diagnostics — require knowing which electrodes are near each other. EegFun finds neighbours based on Euclidean distance within a specified criterion.

### 2D Neighbours

Computed from the 2D projected coordinates. Suitable for topographic analysis:

```julia
EegFun.get_neighbours_xy!(layout, 0.5)
```

### 3D Neighbours

Computed from the 3D spherical coordinates. Preferred for channel repair since distances on the sphere are more accurate than projected distances:

```julia
EegFun.get_neighbours_xyz!(layout, 0.5)
```

> [!IMPORTANT]
> The **distance criterion** is in normalised coordinate units where 1.0 = the scalp equator (ear level). Think of it as a fraction of the head radius: 0.5 means "neighbours within 50% of the head radius." After computing neighbours, check the average count to ensure the value makes sense for your cap density:
>
> ```julia
> avg = EegFun.average_number_of_neighbours(layout.neighbours)
> ```
>
> A good target is typically 4–6 neighbours per electrode, although this can vary depending on the cap density. If the average is much lower, increase the criterion; if much higher, decrease it. But visual inspection is always recommended.

### Inspecting Neighbours

For transparency and reproducibility, you can export the full neighbour structure — including distances and interpolation weights — to a TOML file:

```julia
EegFun.print_layout_neighbours(layout, "neighbours.toml")
```

### Visualising Neighbours

Pass `neighbours = true` to the layout plot to interactively inspect connections. Hovering over an electrode highlights its neighbours:

```julia
EegFun.plot_layout_2d(layout, neighbours = true)
EegFun.plot_layout_3d(layout, neighbours = true)
```

## Checking and Querying

```julia
# Does the layout have computed neighbours?
EegFun.has_neighbours(layout)

# What distance criterion was used?
EegFun.neighbour_criterion(layout)

# Does the layout have valid (non-zero) coordinates?
EegFun.has_valid_coordinates(layout)
```

> [!NOTE]
> Auto-generated layouts (created when loading data without a layout file) have all-zero coordinates. These work for basic data browsing but cannot be used for topographic plots or spatial operations. `has_valid_coordinates` returns `false` for such layouts.

## Summary

A typical layout setup for analysis:

```julia
# Read the layout
layout = EegFun.read_layout("resources/layouts/biosemi/biosemi64.csv")

# Convert coordinates
EegFun.polar_to_cartesian_xy!(layout)   # for topoplots
EegFun.polar_to_cartesian_xyz!(layout)  # for spherical spline repair

# Compute neighbours
EegFun.get_neighbours_xyz!(layout, 0.5)

# Verify visually
EegFun.plot_layout_2d(layout, neighbours = true)
```

## Next Steps

- **Channel repair** — [Artifact Handling](artifact-handling.md) uses the neighbour structure computed here for `repair_channels!` and `repair_artifacts!`
- **Manual preprocessing** — [Manual Preprocessing](manual-preprocessing.md) shows how layout setup fits into the full preprocessing workflow
- **Batch pipelines** — [Batch Processing](batch-processing.md) configures the layout and neighbour criterion via TOML
