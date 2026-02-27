This demo demonstrates how to load, inspect, modify, and manage electrode layouts and their spatial neighbour relationships.

### What are Layouts?

Layouts define the spatial positions of EEG electrodes on the scalp. They are essential for:

- Topographic plotting (scalp maps)
- Spherical spline interpolation (channel repair)
- Cluster-based permutation tests (spatial adjacency)
- Neighbour-based artifact repair

### Coordinate Systems

EegFun supports multiple coordinate representations:

**Polar coordinates**: The raw format stored in layout CSV files (theta/radius pairs)

**2D Cartesian (x, y)**: Used for topographic plots and 2D neighbour calculations (`polar_to_cartesian_xy!`)

**3D Cartesian (x, y, z)**: Used for spherical spline interpolation and 3D neighbour calculations (`polar_to_cartesian_xyz!`)

### Neighbours

Spatial neighbours define which electrodes are "close" to each other. This adjacency information is used by:

- Cluster-based permutation tests (defining spatial clusters)
- Neighbour-based artifact repair (interpolating bad channels from nearby good channels)

Neighbours can be calculated in 2D (`get_neighbours_xy!`) or 3D (`get_neighbours_xyz!`) space with a configurable distance threshold.

## Workflow Summary

This demo covers:

- Loading layout files and converting coordinates
- Calculating and inspecting neighbour relationships
- Subsetting layouts to specific channels
- Renaming channels in a layout
- Creating custom layouts from scratch
- Validating layout structure
