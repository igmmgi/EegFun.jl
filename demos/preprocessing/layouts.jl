# Demo: Layout and Neighbour Management
# Shows how to load, inspect, modify, and manage electrode layouts
# and their spatial neighbour relationships.

using EegFun

#######################################################################
# LOADING AND INSPECTING LAYOUTS
#######################################################################

# Load a layout file (CSV with polar coordinates)
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
layout  # show summary

# Convert polar to Cartesian coordinates for 2D plotting
EegFun.polar_to_cartesian_xy!(layout)

# Convert to 3D Cartesian coordinates (for spherical spline interpolation)
EegFun.polar_to_cartesian_xyz!(layout)

# Check coordinate validity
EegFun.has_valid_coordinates(layout)


#######################################################################
# NEIGHBOUR CALCULATION
#######################################################################

# Calculate 2D neighbours using a distance criterion (in arbitrary units)
EegFun.get_neighbours_xy!(layout, 0.5)
layout  # now shows neighbours

# Check how many neighbours each channel has
EegFun.average_number_of_neighbours(layout.neighbours)

# Print neighbours to screen
EegFun.print_layout_neighbours(layout)

# Save neighbours to a TOML file
# EegFun.print_layout_neighbours(layout, "my_neighbours.toml")

# Check if a specific channel has enough neighbours
EegFun.check_channel_neighbors(layout, :Cz)

# Clear and recalculate with a different criterion
EegFun.clear_neighbours!(layout)
EegFun.get_neighbours_xyz!(layout, 0.5)  # 3D neighbours


#######################################################################
# LAYOUT SUBSETTING
#######################################################################

# Subset layout to only include specific channels
layout_subset = EegFun.subset_layout(layout, [:Fp1, :Fp2, :F3, :F4, :Fz, :Cz, :Pz])
layout_subset  # 7 channels

# In-place subsetting (modifies the original)
layout_copy = copy(layout)
EegFun.subset_layout!(layout_copy, [:Cz, :Pz, :Oz])
layout_copy  # 3 channels


#######################################################################
# RENAMING CHANNELS
#######################################################################

# Rename a single channel (non-mutating)
layout_renamed = EegFun.rename_channel(layout, :Fp1, :FP1_new)

# Rename in-place
layout_copy2 = copy(layout)
EegFun.rename_channel!(layout_copy2, :Fp2, :FP2_new)


#######################################################################
# CUSTOM LAYOUTS
#######################################################################

# Create a custom layout from scratch (e.g., for non-standard montages)
custom_layout = EegFun.create_custom_layout(labels = [:Ch1, :Ch2, :Ch3, :Ch4], x = [0.0, 1.0, 0.0, -1.0], y = [1.0, 0.0, -1.0, 0.0])
custom_layout

# Validate a layout has required structure
EegFun.validate_layout(layout)
