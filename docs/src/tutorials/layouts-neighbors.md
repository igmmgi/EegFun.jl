# Working with Layouts and Neighbors

Spatial analysis and artifact repair in `EegFun.jl` depend on knowing where electrodes are and how they relate to each other. This tutorial explains how to manage layouts and discover electrode neighbors.

## 1. Importing Layouts

A layout is typically stored in a CSV file containing electrode labels and their polar coordinates (incidence and azimuth).

```julia
# Import a standard layout
layout = read_layout("biosemi64.csv")

# Access channel labels
labels = channel_labels(layout)
```

## 2. Coordinate Transformations

EEG systems often provide data in polar coordinates (angles), but plotting and spatial statistics often require Cartesian coordinates ($x, y, z$).

### 2D Transformations (for Topoplots)

To generate a 2D topographic map, you need flat Cartesian coordinates.

```julia
# Convert polar to 2D Cartesian (x, y)
# By default, preserves radial distance (inc=90° maps to radius=1.0)
polar_to_cartesian_xy!(layout)
```

### 3D Transformations (for Repair & Source Localization)

Advanced artifact repair methods (like Spherical Splines) require 3D positions.

```julia
# Convert polar to 3D Cartesian (x, y, z)
polar_to_cartesian_xyz!(layout)
```

## 3. Neighbor Discovery

Many algorithms (spatial clustering, channel repair) require a list of "neighboring" sensors for each electrode. `EegFun.jl` calculates these based on Euclidean distance.

```julia
# Define a distance criterion (in mm or normalized units)
# and find all neighbors within that radius
get_neighbours_xy!(layout, 40.0)

# Check the average number of neighbors per electrode
avg = average_number_of_neighbours(layout)
println("Average neighbors: $avg")
```

The results are stored in an `OrderedDict` within the layout, mapping each electrode to its neighbors, their distances, and interpolation weights.

## 4. Inspecting Neighbors

It is good practice to verify your neighbor structure to ensure the distance criterion makes sense for your electrode density.

```julia
# Print the neighbor structure to a TOML file for review
print_layout_neighbours(layout, "my_neighbors.toml")
```

The generated TOML will show exactly which channels are considered neighbors for every electrode, which is essential for transparent and reproducible analysis.

## Summary

```julia
# Typical layout setup
layout = read_layout("cap.csv")
polar_to_cartesian_xy!(layout)
get_neighbours_xy!(layout, 45.0)

# The layout is now ready for spatial statistics or artifact repair
```
