# Layout type is defined in types.jl

# === LAYOUT METADATA ACCESSORS ===
# Layout metadata group accessors
"""
    channel_column_labels(layout::Layout) -> Vector{Symbol}

Get electrode label column names from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: Electrode label column names

# Examples
```julia
labels = channel_column_labels(layout)
```
"""
channel_column_labels(layout::Layout) = _get_cols_by_group(layout.data, :label)

"""
    channel_column_data(layout::Layout) -> DataFrame

Get electrode label data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing electrode label data

# Examples
```julia
data = channel_column_data(layout)
```
"""
channel_column_data(layout::Layout) = layout.data[:, _get_cols_by_group(layout.data, :label)]

"""
    position_polar_labels(layout::Layout) -> Vector{Symbol}

Get polar coordinate column names from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: Polar coordinate column names (inc, azi)

# Examples
```julia
polar_cols = position_polar_labels(layout)
```
"""
position_polar_labels(layout::Layout) = _get_cols_by_group(layout.data, :polar_coords)

"""
    position_2D_labels(layout::Layout) -> Vector{Symbol}

Get 2D Cartesian coordinate column names from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: 2D Cartesian coordinate column names (x2, y2) if available

# Examples
```julia
cartesian_2d_cols = position_2D_labels(layout)
```
"""
position_2D_labels(layout::Layout) = _get_cols_by_group(layout.data, :cartesian_2d)

"""
    position_3D_labels(layout::Layout) -> Vector{Symbol}

Get 3D Cartesian coordinate column names from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: 3D Cartesian coordinate column names (x3, y3, z3) if available

# Examples
```julia
cartesian_3d_cols = position_3D_labels(layout)
```
"""
position_3D_labels(layout::Layout) = _get_cols_by_group(layout.data, :cartesian_3d)

# Layout coordinate data functions
"""
    positions_polar_data(layout::Layout) -> DataFrame

Get polar coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing polar coordinate columns (inc, azi)

# Examples
```julia
polar_data = positions_polar_data(layout)
```
"""
positions_polar_data(layout::Layout) = layout.data[:, position_polar_labels(layout)]

"""
    positions_2D_data(layout::Layout) -> DataFrame

Get 2D Cartesian coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing 2D Cartesian coordinate columns (x2, y2) if available

# Examples
```julia
cartesian_2d_data = positions_2D_data(layout)
```
"""
positions_2D_data(layout::Layout) = layout.data[:, position_2D_labels(layout)]

"""
    positions_3D_data(layout::Layout) -> DataFrame

Get 3D Cartesian coordinate data from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `DataFrame`: DataFrame containing 3D Cartesian coordinate columns (x3, y3, z3) if available

# Examples
```julia
cartesian_3d_data = positions_3D_data(layout)
```
"""
positions_3D_data(layout::Layout) = layout.data[:, position_3D_labels(layout)]

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
    _add_metadata!(df::DataFrame, columns::Vector{Symbol}, group::Symbol)

Add metadata to DataFrame columns while preserving existing metadata.

This internal function adds metadata group tags to specified columns while
preserving any existing metadata on those columns. It handles missing
columns gracefully and provides warnings for non-existent columns.

# Arguments
- `df::DataFrame`: The DataFrame to add metadata to
- `columns::Vector{Symbol}`: Column names to tag with metadata
- `group::Symbol`: The metadata group to assign (e.g., :label, :channels, :polar_coords)

# Modifies
- `df`: Adds metadata to existing columns

# Examples
```julia
_add_metadata!(df, [:label], :label)
_add_metadata!(df, [:inc, :azi], :polar_coords)
```
"""
function _add_metadata!(df::DataFrame, columns::Vector{Symbol}, group::Symbol)

    # Filter to only existing columns
    existing_cols = [col for col in columns if hasproperty(df, col)]
    if isempty(existing_cols)
        @minimal_error "No existing columns found for group $group"
        return
    end
    
    # Report any missing columns
    missing_cols = setdiff(columns, existing_cols)
    if !isempty(missing_cols)
        @minimal_error "Columns not found in DataFrame: $missing_cols"
    end

    # Set metadata for each existing column
    for col in existing_cols
        col_str = string(col)
        metadata!(df, col_str, "group" => group)
    end

end

# Accessor methods
"""
    get_labels(layout::Layout) -> Vector{Symbol}

Get all electrode labels from the layout.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Vector{Symbol}`: Vector of electrode labels

# Examples
```julia
labels = get_labels(layout)
```
"""
get_labels(layout::Layout) = layout.data.label

"""
    has_2d_coords(layout::Layout) -> Bool

Check if the layout has 2D Cartesian coordinates.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if x2 and y2 columns exist

# Examples
```julia
if has_2d_coords(layout)
    # Use 2D coordinates
end
```
"""
has_2d_coords(layout::Layout) = all(col -> col in propertynames(layout.data), [:x2, :y2])

"""
    has_3d_coords(layout::Layout) -> Bool

Check if the layout has 3D Cartesian coordinates.

# Arguments
- `layout::Layout`: The layout object

# Returns
- `Bool`: True if x3, y3, and z3 columns exist

# Examples
```julia
if has_3d_coords(layout)
    # Use 3D coordinates
end
```
"""
has_3d_coords(layout::Layout) = all(col -> col in propertynames(layout.data), [:x3, :y3, :z3])

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
        @minimal_error "Cannot open file: $file"
    end
    @info "Reading layout from $file"
    df = DataFrame(CSV.File(file, types = Dict(:label => Symbol)))
    rename!(df, Symbol.(names(df)))
    
    # Add metadata for column groups
    _add_metadata!(df, [:label], :label)
    _add_metadata!(df, [:inc, :azi], :polar_coords)
    
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
    polar_to_cartesian_xy!(layout::Layout)

Converts polar coordinates (incidence and azimuth angles) from a layout into Cartesian coordinates (x, y).

# Arguments
- `layout::Layout`: A Layout containing the layout information with columns for incidence angles (`:inc`) and azimuth angles (`:azi`).

# Modifies
- The input `layout` is modified in place to include new columns `x2` and `y2`, which represent the Cartesian coordinates calculated from the polar coordinates.
- Clears any existing neighbour information since coordinates have changed.

# Returns
- Nothing. The function modifies the `layout` directly.
"""
function polar_to_cartesian_xy!(layout::Layout)
    # Get the DataFrame from Layout
    df = layout.data
    
    # Check for required columns
    if !all([col in propertynames(df) for col in [:inc, :azi]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end

    # Validate data types
    if !(eltype(df.inc) <: Number && eltype(df.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    radius = 88 # mm
    inc = df[!, :inc] .* (pi / 180)
    azi = df[!, :azi] .* (pi / 180)

    # Store existing metadata before adding columns
    existing_metadata = copy(DataFrames.metadata(df))
    
    df[!, :x2] = inc .* cos.(azi) .* radius
    df[!, :y2] = inc .* sin.(azi) .* radius

    # Restore existing metadata
    for (col, meta) in existing_metadata
        DataFrames.metadata!(df, col, meta)
    end
    
    # Add metadata for the new 2D Cartesian coordinates
    _add_metadata!(df, [:x2, :y2], :cartesian_2d)

    # Clear neighbours since coordinates have changed
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
    # Get the DataFrame from Layout
    df = layout.data
    
    # Check for required columns
    if !all([col in propertynames(df) for col in [:inc, :azi]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end

    # Validate data types
    if !(eltype(df.inc) <: Number && eltype(df.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    radius = 88.0  # mm
    inc = df[!, :inc] .* (pi / 180)  # Convert to radians
    azi = df[!, :azi] .* (pi / 180)  # Convert to radians

    # Store existing metadata before adding columns
    existing_metadata = copy(DataFrames.metadata(df))
    
    # Standard spherical to Cartesian conversion
    df[!, :x3] = radius .* sin.(inc) .* cos.(azi)
    df[!, :y3] = radius .* sin.(inc) .* sin.(azi)
    df[!, :z3] = radius .* cos.(inc)

    # Restore existing metadata
    for (col, meta) in existing_metadata
        DataFrames.metadata!(df, col, meta)
    end

    # Add metadata for the new 3D Cartesian coordinates
    _add_metadata!(df, [:x3, :y3, :z3], :cartesian_3d)

    # Clear neighbours since coordinates have changed
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
    get_electrode_neighbours_xy!(layout::Layout, distance_criterion::Real)

Identifies the neighbours of each electrode based on their Cartesian coordinates.

# Arguments
- `layout::Layout`: A Layout containing the layout information with columns for electrode labels, Cartesian coordinates (`x2`, `y2`).
- `distance_criterion::Real`: The maximum distance to consider two electrodes as neighbours.

# Returns
- `OrderedDict{Symbol,Neighbours}`: A dictionary where each key is an electrode label, and the value is a Neighbours struct containing neighbour information.

# Throws
- `ArgumentError`: If the layout does not contain the required columns.
"""
function get_layout_neighbours_xy!(layout::Layout, distance_criterion::Real)

    if distance_criterion <= 0
        @minimal_error "Distance criterion must be positive"
    end
    _ensure_coordinates_2d!(layout)

    @info "Calculating neighbours with distance criterion $distance_criterion mm"
    # Precompute coordinates
    coords = Matrix{Float64}(undef, size(layout.data, 1), 2)
    coords[:, 1] = layout.data.x2
    coords[:, 2] = layout.data.y2

    neighbour_dict = OrderedDict{Symbol,Neighbours}()

    for (idx1, label1) in enumerate(layout.data.label)

        neighbour_dict[Symbol(label1)] = Neighbours([], [], [])

        for (idx2, label2) in enumerate(layout.data.label)
            if idx1 == idx2
                continue
            end

            # Compute squared distance
            distance_sq = squared_distance_xy(coords[idx1, 1], coords[idx1, 2], coords[idx2, 1], coords[idx2, 2])

            if distance_sq <= distance_criterion^2
                distance = sqrt(distance_sq)
                push!(neighbour_dict[Symbol(label1)].electrodes, label2)
                push!(neighbour_dict[Symbol(label1)].distances, distance)
            end
        end

        # Compute weights (inverse distance weighting)
        distances = neighbour_dict[Symbol(label1)].distances
        inv_distances = 1 ./ distances  # Inverse of distances
        total_inv_distance = sum(inv_distances)
        for idx in eachindex(neighbour_dict[Symbol(label1)].electrodes)
            push!(neighbour_dict[Symbol(label1)].weights, inv_distances[idx] / total_inv_distance)
        end

    end

    layout.neighbours = neighbour_dict
    layout.criterion = distance_criterion

    return nothing

end

"""
    get_electrode_neighbours_xyz!(layout::Layout, distance_criterion::Real)

Identifies the neighbours of each electrode based on their Cartesian coordinates.

# Arguments
- `layout::Layout`: A Layout containing the layout information with columns for electrode labels, Cartesian coordinates (`x3`, `y3`, `z3`).
- `distance_criterion::Real`: The maximum distance to consider two electrodes as neighbours.

# Returns
- `OrderedDict{Symbol,Neighbours}`: A dictionary where each key is an electrode label, and the value is a Neighbours struct containing neighbour information.

# Throws
- `ArgumentError`: If the layout does not contain the required columns.
"""
function get_layout_neighbours_xyz!(layout::Layout, distance_criterion::Real)

    if distance_criterion <= 0
        @minimal_error "Distance criterion must be positive"
    end
    _ensure_coordinates_3d!(layout)

    @info "Calculating neighbours with distance criterion $distance_criterion mm"
    # Precompute coordinates
    coords = Matrix{Float64}(undef, size(layout.data, 1), 3)
    coords[:, 1] = layout.data.x3
    coords[:, 2] = layout.data.y3
    coords[:, 3] = layout.data.z3

    neighbour_dict = OrderedDict{Symbol, Neighbours}()

    for (idx1, label1) in enumerate(layout.data.label)

        neighbour_dict[Symbol(label1)] = Neighbours([], [], [])

        for (idx2, label2) in enumerate(layout.data.label)

            if idx1 == idx2
                continue
            end

            # Compute squared distance
            distance_sq = squared_distance_xyz(
                coords[idx1, 1],
                coords[idx1, 2],
                coords[idx1, 3],
                coords[idx2, 1],
                coords[idx2, 2],
                coords[idx2, 3],
            )

            if distance_sq <= distance_criterion^2
                distance = sqrt(distance_sq)
                push!(neighbour_dict[Symbol(label1)].electrodes, label2)
                push!(neighbour_dict[Symbol(label1)].distances, distance)
            end

        end

        # Compute weights (inverse distance weighting)
        distances = neighbour_dict[Symbol(label1)].distances
        inv_distances = 1 ./ distances  # Inverse of distances
        total_inv_distance = sum(inv_distances)
        for idx in eachindex(neighbour_dict[Symbol(label1)].electrodes)
            push!(neighbour_dict[Symbol(label1)].weights, inv_distances[idx] / total_inv_distance)
        end

    end

    layout.neighbours = neighbour_dict
    layout.criterion = distance_criterion

    return nothing

end

# === NEIGHBOR UTILITIES ===
"""
    _format_neighbours_toml(neighbours_dict::OrderedDict{Symbol, Neighbours}, nneighbours::Real)

Helper function to format the neighbours dictionary as TOML structure.
Preserves the order of electrodes from the original OrderedDict.
"""
function _format_neighbours_toml(neighbours_dict::OrderedDict{Symbol, Neighbours}, nneighbours::Real)
    toml_dict = OrderedDict{String, Any}()
    
    # Add metadata first
    toml_dict["metadata"] = OrderedDict(
        "average_neighbours_per_electrode" => round(nneighbours, digits=2),
        "total_electrodes" => length(neighbours_dict),
        "generated_at" => string(now())
    )
    
    # Add electrode data in original order
    toml_dict["electrodes"] = OrderedDict{String, Any}()
    
    for (electrode, neighbours) in neighbours_dict
        electrode_str = string(electrode)
        toml_dict["electrodes"][electrode_str] = OrderedDict(
            "neighbours" => [string(n) for n in neighbours.electrodes],
            "distances" => [round(d, digits=4) for d in neighbours.distances],
            "weights" => [round(w, digits=6) for w in neighbours.weights],
            "neighbour_count" => length(neighbours.electrodes)
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
neighbours = get_layout_neighbours_xy!(layout, 40)

# Write to TOML file
print_neighbours_dict(neighbours, "neighbours.toml")
```
"""
function print_layout_neighbours(neighbours_dict::OrderedDict{Symbol, Neighbours}, filename::String)
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
neighbours, _ = get_layout_neighbours_xy!(layout, 40)
avg_neighbours = average_neighbours_per_electrode(neighbours)
```
"""
function average_number_of_neighbours(neighbours_dict::OrderedDict{Symbol, Neighbours})
    if isempty(neighbours_dict)
        return 0.0
    end
    total_neighbours = sum(length(neighbours.electrodes) for neighbours in values(neighbours_dict))
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
    n_neighbours = length(neighbours.electrodes)
    avg_distance = round(mean(neighbours.distances), digits=1)
    
    println(io, "$(rpad(string(electrode), 6)): $(n_neighbours) neighbours (avg dist: $(avg_distance)mm)")
    if n_neighbours > 0
        neighbour_details = []
        for j in 1:n_neighbours
            neighbour = string(neighbours.electrodes[j])
            distance = round(neighbours.distances[j], digits=1)
            weight = round(neighbours.weights[j], digits=3)
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
        layout.neighbours = get_layout_neighbours_xy!(layout, distance_criterion)
        layout.criterion = distance_criterion
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

# Modifies
- `layout`: Updates neighbor information and criterion

# Examples
```julia
get_neighbours_xyz!(layout, 40.0)
```
"""
function get_neighbours_xyz!(layout::Layout, distance_criterion::Real)
    if !has_neighbours(layout) || layout.criterion != distance_criterion
        layout.neighbours = get_layout_neighbours_xyz!(layout, distance_criterion)
        layout.criterion = distance_criterion
    end
    return nothing
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
    layout.data = filter(:label => in(selected_channels), layout.data)
    
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