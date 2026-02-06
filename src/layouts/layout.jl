# Layout type is defined in types.jl
"""
    channel_labels(layout::Layout) -> Vector{Symbol}

Get electrode label column names from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: Electrode label column names
"""
channel_labels(layout::Layout)::Vector{Symbol} = layout.data.label

"""
    positions_polar(layout::Layout) -> DataFrame

Get polar coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing polar coordinate columns (inc, azi)
"""
positions_polar(layout::Layout)::DataFrame = layout.data[:, [:inc, :azi]]

"""
    positions_2D(layout::Layout) -> DataFrame

Get 2D Cartesian coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing 2D Cartesian coordinate columns (x2, y2) if available
"""
positions_2D(layout::Layout)::DataFrame = layout.data[:, [:x2, :y2]]

"""
    positions_3D(layout::Layout) -> DataFrame

Get 3D Cartesian coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing 3D Cartesian coordinate columns (x3, y3, z3) if available
"""
positions_3D(layout::Layout)::DataFrame = layout.data[:, [:x3, :y3, :z3]]




# === LAYOUT VALIDATION AND MANIPULATION ===
"""
    validate_layout(layout::Layout) -> Layout

Validate that a layout contains all required columns and data.

This function checks that the layout DataFrame contains the minimum required
columns for EEG layout operations. It validates both column presence and
data type requirements.

# Arguments
- `layout::Layout`: The layout object to validate

# Returns
- `Layout`: The validated layout object

# Throws
- `ArgumentError`: If required columns are missing

# Required Columns
- `:label`: Electrode labels (Symbol type)
- `:inc`: Incidence angles in degrees
- `:azi`: Azimuth angles in degrees

# Examples
```julia
layout = read_layout("biosemi64.csv")
validated_layout = validate_layout(layout)
```
"""
function validate_layout(layout::Layout)
    required_cols = [:label, :inc, :azi]
    missing_cols = setdiff(required_cols, propertynames(layout.data))
    if !isempty(missing_cols)
        throw(ArgumentError("Layout missing required columns: $missing_cols"))
    end
    return layout
end

"""
    has_2d_coords(layout::Layout) -> Bool

Check if the layout has 2D Cartesian coordinates.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if x2 and y2 columns exist
"""
has_2d_coords(layout::Layout) = all(col -> col in propertynames(layout.data), [:x2, :y2])

"""
    has_3d_coords(layout::Layout) -> Bool

Check if the layout has 3D Cartesian coordinates.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if x3, y3, and z3 columns exist
"""
has_3d_coords(layout::Layout) = all(col -> col in propertynames(layout.data), [:x3, :y3, :z3])
"""
    has_valid_coordinates(layout::Layout) -> Bool

Check if the layout has non-zero coordinates suitable for topographic plotting.
Auto-generated layouts have all zeros which cannot be used for spatial visualization.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if layout has non-zero inc/azi coordinates

# Examples
```julia
if !has_valid_coordinates(layout)
    error("Cannot create topographic plot without a proper electrode layout")
end
```
"""
function has_valid_coordinates(layout::Layout)::Bool
    # Check if layout has inc and azi columns
    if !all(col -> col in propertynames(layout.data), [:inc, :azi])
        return false
    end
    # Check if any coordinates are non-zero
    return any(layout.data.inc .!= 0) || any(layout.data.azi .!= 0)
end

# === LAYOUT I/O FUNCTIONS ===
"""
    read_layout(layout_file_name::String)

Reads a layout file in CSV format and returns a DataFrame containing the layout information.

# Arguments
- `layout_file_name::String`: The path to the CSV file containing the layout data.

# Returns
- `DataFrame`: A DataFrame containing the layout data, with columns for electrode labels, incidence angles, azimuth angles, and calculated Cartesian coordinates (if applicable).

# Throws
- `SystemError`: If the file does not exist or cannot be accessed.
"""
function read_layout(file)
    if !isfile(file)
        @minimal_error_throw "Cannot open file: $file"
    end
    @info "Reading layout from $file"
    df = DataFrame(CSV.File(file, types = Dict(:label => Symbol)))
    rename!(df, Symbol.(names(df)))

    return Layout(df, nothing, nothing)
end
# === COORDINATE CONVERSIONS ===
"""
    _ensure_coordinates_2d!(layout::Layout)

Helper function to ensure 2D coordinates exist in the layout DataFrame.
Converts polar coordinates to Cartesian if needed.

# Arguments
- `layout::Layout`: The layout to check and potentially convert

# Modifies
- `layout`: Adds 2D Cartesian coordinates if they don't exist

# Examples
```julia
_ensure_coordinates_2d!(layout)
```
"""
function _ensure_coordinates_2d!(layout::Layout)
    if !has_2d_coords(layout)
        @info "Converting polar coordinates to 2D Cartesian coordinates"
        polar_to_cartesian_xy!(layout)
    end
end

"""
    _ensure_coordinates_3d!(layout::Layout)

Helper function to ensure 3D coordinates exist in the layout DataFrame.
Converts polar coordinates to Cartesian if needed.

# Arguments
- `layout::Layout`: The layout to check and potentially convert

# Modifies
- `layout`: Adds 3D Cartesian coordinates if they don't exist

# Examples
```julia
_ensure_coordinates_3d!(layout)
```
"""
function _ensure_coordinates_3d!(layout::Layout)
    if !has_3d_coords(layout)
        @info "Converting polar coordinates to 3D Cartesian coordinates"
        polar_to_cartesian_xyz!(layout)
    end
end

"""
    polar_to_cartesian_xy!(layout::Layout; normalization_radius::Float64=1.0, preserve_radial_distance::Bool=true)

Converts polar coordinates (incidence and azimuth angles) from a layout into Cartesian coordinates (x, y).

# Arguments
- `layout::Layout`: A Layout containing the layout information with columns for incidence angles (`:inc`) and azimuth angles (`:azi`).
- `normalization_radius::Float64`: Maximum radius for electrode positions (default: 1.0). Only used when `preserve_radial_distance=false`.
- `preserve_radial_distance::Bool`: If true, preserves true radial distances from polar coordinates, allowing electrodes with inc>90° to appear outside the unit circle (default: true). When false, normalizes all electrodes to fit within normalization_radius.

# Modifies
- The input `layout` is modified in place to include new columns `x2` and `y2`, which represent the Cartesian coordinates calculated from the polar coordinates.
- Clears any existing neighbour information since coordinates have changed.

# Returns
- Nothing. The function modifies the `layout` directly.

# Notes
- When `preserve_radial_distance=true`, electrodes imported from EEGLAB with inc>90° (e.g., eye electrodes beyond the scalp) will appear outside the standard head circle, matching EEGLAB's topoplot rendering.
"""
function polar_to_cartesian_xy!(layout::Layout; normalization_radius::Float64 = 1.0, preserve_radial_distance::Bool = true)
    df = layout.data

    # Check for required columns and numeric types
    if !all(col -> col in propertynames(df), [:inc, :azi])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end
    if !(eltype(df.inc) <: Number && eltype(df.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    # Vectorized conversion (degrees to radians + polar to cartesian)
    rad_factor = pi / 180
    inc = df.inc .* rad_factor
    azi = df.azi .* rad_factor

    x2 = inc .* cos.(azi)
    y2 = inc .* sin.(azi)

    # Note: We do NOT center the coordinates here!
    # Polar coordinates are naturally centered at the vertex (inc=0°).
    # Centering would shift vertices like CZ away from (0,0), especially when 
    # eye electrodes or other non-scalp sensors are present.

    if !preserve_radial_distance
        # Normalize to fit all electrodes within normalization_radius (traditional behavior)
        max_r = maximum(sqrt.(x2 .^ 2 .+ y2 .^ 2))
        if max_r > 0
            scale = normalization_radius / max_r
            df[!, :x2] = x2 .* scale
            df[!, :y2] = y2 .* scale
        else
            df[!, :x2] = x2
            df[!, :y2] = y2
        end
    else
        # Preserve true radial distances - electrodes with inc>90° will be outside unit circle
        # Map inc in degrees to radius: 90° = 1.0 (equator/unit circle)
        scale_factor = 1.0 / (π / 2)  # Scale so that inc=90° (π/2 radians) → radius=1.0
        df[!, :x2] = x2 .* scale_factor
        df[!, :y2] = y2 .* scale_factor
    end

    clear_neighbours!(layout)
    return nothing
end

"""
    polar_to_cartesian_xyz!(layout::Layout)

Converts polar coordinates (incidence and azimuth angles) from a layout into Cartesian coordinates (x, y, z).

# Arguments
- `layout::Layout`: A Layout containing the layout information with columns for incidence angles (`:inc`) and azimuth angles (`:azi`).

# Modifies
- The input `layout` is modified in place to include new columns `x3`, `y3`, and `z3`, which represent the Cartesian coordinates calculated from the polar coordinates.
- Clears any existing neighbour information since coordinates have changed.

# Returns
- Nothing. The function modifies the `layout` directly.
"""
function polar_to_cartesian_xyz!(layout::Layout)
    df = layout.data

    if !all(col -> col in propertynames(df), [:inc, :azi])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end
    if !(eltype(df.inc) <: Number && eltype(df.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    rad_factor = pi / 180
    inc = df.inc .* rad_factor
    azi = df.azi .* rad_factor

    x3 = sin.(inc) .* cos.(azi)
    y3 = sin.(inc) .* sin.(azi)
    z3 = cos.(inc)

    # Optimize range and centering
    min_x, max_x = extrema(x3)
    min_y, max_y = extrema(y3)
    min_z, max_z = extrema(z3)

    max_range = max(max_x - min_x, max_y - min_y, max_z - min_z)

    if max_range > 0
        half_range = max_range / 2
        df[!, :x3] = (x3 .- (max_x + min_x) / 2) ./ half_range
        df[!, :y3] = (y3 .- (max_y + min_y) / 2) ./ half_range
        df[!, :z3] = (z3 .- (max_z + min_z) / 2) ./ half_range
    else
        df[!, :x3] = x3
        df[!, :y3] = y3
        df[!, :z3] = z3
    end

    clear_neighbours!(layout)
    return nothing
end

# === DISTANCE CALCULATIONS ===
"""
    calculate_distance_xy(x1::Real, y1::Real, x2::Real, y2::Real)

Calculates the Euclidean distance between two points in 2D space.

# Arguments
- `x1, y1, x2, y2`: The coordinates of the two points.

# Returns
- `Float64`: The distance between the two points.
"""
@inline function distance_xy(x1::Real, y1::Real, x2::Real, y2::Real)::Float64
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

"""
    squared_distance_xy(x1, y1, x2, y2)

Calculates the squared distance between two points in 2D space.

# Arguments
- `x1, y1, x2, y2`: The coordinates of the two points.

# Returns
- `Float64`: The squared distance between the two points.
"""
@inline function squared_distance_xy(x1, y1, x2, y2)
    return (x1 - x2)^2 + (y1 - y2)^2
end

"""
    calculate_distance_xyz(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)

Calculates the Euclidean distance between two points in 3D space.

# Arguments
- `x1, y1, z1, x2, y2, z2`: The coordinates of the two points.

# Returns
- `Float64`: The distance between the two points.
"""
@inline function distance_xyz(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)::Float64
    return sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
end

"""
    squared_distance_xyz(x1, y1, z1, x2, y2, z2)

Calculates the squared distance between two points in 3D space.

# Arguments
- `x1, y1, z1, x2, y2, z2`: The coordinates of the two points.

# Returns
- `Float64`: The squared distance between the two points.
"""
@inline function squared_distance_xyz(x1, y1, z1, x2, y2, z2)
    return (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2
end

# === NEIGHBOR CALCULATIONS ===
"""
    get_neighbours_xy!(layout::Layout, distance_criterion::Real)

Identifies the neighbours of each electrode based on their normalized 2D Cartesian coordinates.

# Arguments
- `layout::Layout`: A Layout containing the layout information.
- `distance_criterion::Real`: The maximum distance in normalized units to consider two electrodes as neighbours.

# Modifies
- The input `layout` is modified in place to update its `neighbours` field and `criterion`.
"""

"""
    _find_neighbours(labels, coords, distance_criterion) -> OrderedDict

Internal helper to identify neighbours utilizing distance symmetry.
"""
function _find_neighbours(labels, coords, distance_criterion)
    n = length(labels)
    crit_sq = distance_criterion^2

    # Initialize adjacency-like structure
    # electrode_index => (neighbour_indices, distances)
    adj = [(Int[], Float64[]) for _ = 1:n]

    # Nested loop exploiting symmetry
    for i = 1:n
        p1 = coords[i]
        for j = (i+1):n
            p2 = coords[j]
            # dist_sq using GeometryBasics points is fast
            d_sq = sum((p1 .- p2) .^ 2)

            if d_sq <= crit_sq
                d = sqrt(d_sq)
                # Store for both i and j
                push!(adj[i][1], j)
                push!(adj[i][2], d)
                push!(adj[j][1], i)
                push!(adj[j][2], d)
            end
        end
    end

    # Transform to OrderedDict with weights
    neighbour_dict = OrderedDict{Symbol,Neighbours}()

    for i = 1:n
        label = Symbol(labels[i])
        nb_indices, dists = adj[i]

        if isempty(nb_indices)
            neighbour_dict[label] = Neighbours([], [], [])
            continue
        end

        nb_labels = [Symbol(labels[idx]) for idx in nb_indices]

        # Calculate inverse distance weights
        # Handle potential zero distances to avoid Inf
        inv_d = 1.0 ./ max.(dists, 1e-6)
        total_inv = sum(inv_d)
        weights = inv_d ./ total_inv

        neighbour_dict[label] = Neighbours(nb_labels, dists, weights)
    end

    return neighbour_dict
end

# === NEIGHBOR UTILITIES ===
"""
    _format_neighbours_toml(neighbours_dict::OrderedDict{Symbol, Neighbours}, nneighbours::Real)

Helper function to format the neighbours dictionary as TOML structure.
Preserves the order of electrodes from the original OrderedDict.
"""
function _format_neighbours_toml(neighbours_dict::OrderedDict{Symbol,Neighbours}, nneighbours::Real)
    toml_dict = OrderedDict{String,Any}()

    # Add metadata first
    toml_dict["metadata"] = OrderedDict(
        "average_neighbours_per_electrode" => round(nneighbours, digits = 2),
        "total_electrodes" => length(neighbours_dict),
        "generated_at" => string(now()),
    )

    # Add electrode data in original order
    toml_dict["electrodes"] = OrderedDict{String,Any}()

    for (electrode, neighbours) in neighbours_dict
        electrode_str = string(electrode)
        toml_dict["electrodes"][electrode_str] = OrderedDict(
            "neighbours" => [string(n) for n in neighbours.channels],
            "distances" => [round(d, digits = 4) for d in neighbours.distances],
            "weights" => [round(w, digits = 6) for w in neighbours.weights],
            "neighbour_count" => length(neighbours.channels),
        )
    end

    return toml_dict
end

"""
    print_neighbours_dict(neighbours_dict::OrderedDict{Symbol, Neighbours}, filename::String)

Write the neighbors dictionary to a TOML file in a structured format, showing for each electrode:
- Its neighbors
- The distances to each neighbor
- The weights used for interpolation
- The average number of neighbors per electrode

The electrode order from the original OrderedDict is preserved.

# Arguments
- `neighbours_dict::OrderedDict{Symbol, Neighbours}`: Dictionary returned by get_electrode_neighbours_xy/xyz
- `filename::String`: Path to the output TOML file

# Example
```julia
layout = read_layout("./layouts/biosemi64.csv")
neighbours = get_neighbours_xy!(layout, 40)

# Write to TOML file
print_neighbours_dict(neighbours, "neighbours.toml")
```
"""
function print_layout_neighbours(neighbours_dict::OrderedDict{Symbol,Neighbours}, filename::String)
    nneighbours = average_number_of_neighbours(neighbours_dict)
    toml_data = _format_neighbours_toml(neighbours_dict, nneighbours)
    @info "Printing neighbours to $filename"
    open(filename, "w") do io
        TOML.print(io, toml_data)
    end
end

function print_layout_neighbours(layout::Layout, filename::String)
    if isnothing(layout.neighbours)
        @minimal_error "No neighbours to print"
    end
    print_layout_neighbours(layout.neighbours, filename)
end

"""
    average_number_of_neighbours(neighbours_dict::OrderedDict{Symbol, Neighbours})

Calculate the average number of neighbours per electrode from a neighbours dictionary.

# Arguments
- `neighbours_dict::OrderedDict{Symbol, Neighbours}`: Dictionary containing neighbour information for each electrode

# Returns
- `Float64`: The average number of neighbours per electrode

# Example
```julia
layout = read_layout("./layouts/biosemi64.csv")
neighbours, _ = get_neighbours_xy!(layout, 40)
avg_neighbours = average_neighbours_per_electrode(neighbours)
```
"""
function average_number_of_neighbours(neighbours_dict::OrderedDict{Symbol,Neighbours})
    if isempty(neighbours_dict)
        return 0.0
    end
    total_neighbours = sum(length(neighbours.channels) for neighbours in values(neighbours_dict))
    return total_neighbours / length(neighbours_dict)
end

# Helper methods for neighbour management
"""
    has_neighbours(layout::Layout) -> Bool

Check if the layout has calculated neighbor information.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if neighbors have been calculated and stored

# Examples
```julia
if has_neighbours(layout)
    # Use neighbor information
end
```
"""
has_neighbours(layout::Layout) = !isnothing(layout.neighbours)

"""
    neighbour_criterion(layout::Layout) -> Union{Nothing, Float64}

Get the distance criterion used for neighbor calculations.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Union{Nothing, Float64}`: The distance criterion in mm, or nothing if not calculated

# Examples
```julia
criterion = neighbour_criterion(layout)
```
"""
neighbour_criterion(layout::Layout) = layout.criterion

"""
    neighbours(layout::Layout) -> Union{Nothing, OrderedDict{Symbol, Neighbours}}

Get the neighbor information for all electrodes.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Union{Nothing, OrderedDict{Symbol, Neighbours}}`: Dictionary of neighbors for each electrode, or nothing if not calculated

# Examples
```julia
neighbor_dict = neighbours(layout)
```
"""
neighbours(layout::Layout) = layout.neighbours

"""
    clear_neighbours!(layout::Layout)

Clear all neighbor information from the layout.

This function removes all calculated neighbor data and resets the distance
criterion. This is useful when layout coordinates change and neighbor
calculations need to be redone.

# Arguments
- `layout::Layout`: The layout object to clear neighbors from

# Modifies
- `layout`: Removes neighbor data and criterion

# Examples
```julia
clear_neighbours!(layout)
```
"""
function clear_neighbours!(layout::Layout)
    layout.neighbours = nothing
    layout.criterion = nothing
end

"""
    _format_electrode(io, electrode, neighbours)

Helper function to format a single electrode entry for display.

This internal function formats the neighbor information for a single
electrode in a human-readable format, including neighbor count,
average distance, and detailed neighbor information.

# Arguments
- `io`: Output stream for writing
- `electrode`: Electrode symbol
- `neighbours`: Neighbours struct containing neighbor information

# Examples
```julia
_format_electrode(io, :Fp1, neighbours_struct)
```
"""
function _format_electrode(io, electrode, neighbours)
    n_neighbours = length(neighbours.channels)
    avg_distance = round(mean(neighbours.distances), digits = 1)

    println(io, "$(rpad(string(electrode), 6)): $(n_neighbours) neighbours (avg dist: $(avg_distance))")
    if n_neighbours > 0
        neighbour_details = []
        for j = 1:n_neighbours
            neighbour = string(neighbours.channels[j])
            distance = round(neighbours.distances[j], digits = 1)
            weight = round(neighbours.weights[j], digits = 3)
            push!(neighbour_details, "$(neighbour) ($(distance),$(weight))")
        end
        neighbour_list = join(neighbour_details, ", ")
        println(io, "    Neighbours: $neighbour_list")
    end
    println(io)
end

# Get neighbours with automatic computation if needed (mutating)
"""
    get_neighbours_xy!(layout::Layout, distance_criterion::Real)

Get neighbors with automatic computation if needed (mutating version).

This function automatically calculates 2D neighbors if they don't exist or
if the distance criterion has changed. It caches the results for efficiency.

# Arguments
- `layout::Layout`: The layout object
- `distance_criterion::Real`: The distance criterion for neighbors in mm

# Modifies
- `layout`: Updates neighbor information and criterion

# Examples
```julia
get_neighbours_xy!(layout, 40.0)
```
"""
function get_neighbours_xy!(layout::Layout, distance_criterion::Real)
    if !has_neighbours(layout) || layout.criterion != distance_criterion
        if distance_criterion <= 0
            @minimal_error "Distance criterion must be positive"
        end
        _ensure_coordinates_2d!(layout)

        @info "Calculating 2D neighbours with distance criterion $distance_criterion"

        # Precompute coordinates as Point2f for efficiency
        coords = [Point2f(x, y) for (x, y) in zip(layout.data.x2, layout.data.y2)]

        layout.neighbours = _find_neighbours(layout.data.label, coords, distance_criterion)
        layout.criterion = Float64(distance_criterion)
    end
    return nothing
end

"""
    get_neighbours_xyz!(layout::Layout, distance_criterion::Real)

Get neighbors with automatic computation if needed (mutating version).

This function automatically calculates 3D neighbors if they don't exist or
if the distance criterion has changed. It caches the results for efficiency.

# Arguments
- `layout::Layout`: The layout object
- `distance_criterion::Real`: The distance criterion for neighbors in mm
"""
function get_neighbours_xyz!(layout::Layout, distance_criterion::Real)
    if !has_neighbours(layout) || layout.criterion != distance_criterion
        if distance_criterion <= 0
            @minimal_error "Distance criterion must be positive"
        end
        _ensure_coordinates_3d!(layout)

        @info "Calculating 3D neighbours with distance criterion $distance_criterion"

        # Precompute coordinates as Point3f for efficiency
        coords = [Point3f(x, y, z) for (x, y, z) in zip(layout.data.x3, layout.data.y3, layout.data.z3)]

        layout.neighbours = _find_neighbours(layout.data.label, coords, distance_criterion)
        layout.criterion = Float64(distance_criterion)
    end
    return nothing
end

# === LAYOUT VALIDATION AND MANIPULATION ===
"""
    rename!(layout::Layout, rename_dict::Dict{Symbol, Symbol})

Rename channels in a Layout object using a dictionary mapping old names to new names.
Modifies the layout in place.

# Arguments
- `layout::Layout`: The layout object to modify
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `nothing` (modifies the layout in place)

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
rename!(layout, rename_dict)

# Rename a single channel
rename!(layout, Dict(:Cz => :Cz_new))
```

# Notes
- Only channels that exist in the layout will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- Clears any cached neighbour information since channel names have changed
"""
function rename_channel!(layout::Layout, rename_dict::Dict{Symbol,Symbol})
    # Check if any channels in the rename_dict exist in the layout
    existing_channels = Set(layout.data.label)
    channels_to_rename = keys(rename_dict)
    channels_found = intersect(existing_channels, channels_to_rename)

    if isempty(channels_found)
        @info "rename!: No channels found to rename"
        return nothing
    end

    # Check for potential duplicate names before applying any renames
    final_names = Symbol[]
    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels
            push!(final_names, new_name)
        end
    end

    # Check for duplicates in final names
    if length(final_names) != length(unique(final_names))
        duplicate_names = filter(x -> count(==(x), final_names) > 1, unique(final_names))
        @minimal_error_throw "Cannot rename channels to duplicate names: $(join(duplicate_names, ", "))"
    end

    # Apply the renaming with proper swap handling
    # First, collect all the final rename mappings to avoid interference
    final_renames = Dict{Int,Symbol}()  # row_index => final_name

    for (old_name, new_name) in rename_dict
        if old_name ∈ existing_channels
            # Find the row index for this channel
            row_idx = findfirst(==(old_name), layout.data.label)
            if !isnothing(row_idx)
                final_renames[row_idx] = new_name
            end
        end
    end

    # Now apply all renames simultaneously
    for (row_idx, final_name) in final_renames
        layout.data[row_idx, :label] = final_name
    end

    # Clear any cached neighbour information since channel names have changed
    if has_neighbours(layout)
        @info "rename!: Clearing neighbours since channel names have changed"
        clear_neighbours!(layout)
    end

    @info "rename!: Renamed $(length(channels_found)) channels"
    return nothing
end

"""
    rename_channel(layout::Layout, rename_dict::Dict{Symbol, Symbol})

Create a renamed copy of a Layout object using a dictionary mapping old names to new names.

# Arguments
- `layout::Layout`: The layout object to rename
- `rename_dict::Dict{Symbol, Symbol}`: Dictionary mapping old channel names to new names

# Returns
- `Layout`: A new layout object with renamed channels

# Examples
```julia
# Rename Fp1 to Fpz and Fp2 to Fpz
rename_dict = Dict(:Fp1 => :Fpz, :Fp2 => :Fpz)
new_layout = rename(layout, rename_dict)
```

# Notes
- Only channels that exist in the layout will be renamed
- If multiple channels would be renamed to the same name, an error is thrown to prevent duplicates
- The original layout is not modified
"""
function rename_channel(layout::Layout, rename_dict::Dict{Symbol,Symbol})
    # Create a copy of the layout and apply renaming
    renamed_layout = Layout(copy(layout.data), nothing, nothing)
    rename_channel!(renamed_layout, rename_dict)
    return renamed_layout
end

# === LAYOUT SUBSETTING ===
"""
    get_selected_channels(layout::Layout, channel_selection::Function) -> Vector{Symbol}

Get channel labels that match the selection criteria.

This helper function applies a channel selection predicate to the layout's
electrode labels and returns the matching channels.

# Arguments
- `layout::Layout`: The layout object
- `channel_selection::Function`: Channel selection predicate function

# Returns
- `Vector{Symbol}`: Vector of selected channel labels

# Examples
```julia
selected = get_selected_channels(layout, channels([:Fp1, :Fp2]))
```
"""
function get_selected_channels(layout::Layout, channel_selection::Function)
    # Get all channel labels from the layout and apply the channel selection predicate
    all_channels = layout.data.label
    return all_channels[channel_selection(all_channels)]
end

"""
    subset_layout!(layout::Layout; channel_selection = channels())

Subset a Layout object to include only channels that match the channel selection predicate.
Modifies the layout in place.

# Arguments
- `layout::Layout`: The layout to subset
- `channel_selection::Function`: Channel selection predicate function (default: all channels)

# Returns
- `nothing` (modifies the layout in place)
"""
function subset_layout!(layout::Layout; channel_selection = channels())

    n_channels_before = length(layout.data.label)
    selected_channels = get_selected_channels(layout, channel_selection)

    # Filter the layout data to keep only selected channels
    layout.data = DataFrames.filter(:label => in(selected_channels), layout.data)

    # Clear any cached neighbour information since channels have changed
    if has_neighbours(layout)
        @info "subset_layout!: Clearing neighbours since channels have changed"
        clear_neighbours!(layout)
    end

    @info "subset_layout!: Subset layout from $(n_channels_before) to $(length(selected_channels)) channels"
end

"""
    subset_layout(layout::Layout; channel_selection = channels())

Create a subset copy of a Layout object using channel selection predicates.

# Arguments
- `layout::Layout`: The layout to subset
- `channel_selection::Function`: Channel selection predicate function (default: all channels)

# Returns
- `Layout`: A new subset layout object
"""
function subset_layout(layout::Layout; channel_selection = channels())
    # Create a copy of the layout and apply subsetting
    subset_layout = Layout(copy(layout.data), nothing, nothing)
    subset_layout!(subset_layout, channel_selection = channel_selection)
    return subset_layout
end
