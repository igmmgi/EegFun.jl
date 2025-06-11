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
        throw(SystemError("Cannot open file: $file"))
    end
    df = DataFrame(CSV.File(file, types = Dict(:label => Symbol)))
    rename!(df, Symbol.(names(df)))
    return df
end


"""
    polar_to_cartesian_xy!(layout::DataFrame)

Converts polar coordinates (incidence and azimuth angles) from a layout DataFrame into Cartesian coordinates (x, y).

# Arguments
- `layout::DataFrame`: A DataFrame containing the layout information with columns for incidence angles (`:inc`) and azimuth angles (`:azi`).

# Modifies
- The input `layout` DataFrame is modified in place to include new columns `x2` and `y2`, which represent the Cartesian coordinates calculated from the polar coordinates.

# Returns
- Nothing. The function modifies the `layout` DataFrame directly.
"""
function polar_to_cartesian_xy!(layout::DataFrame)
    # Check for required columns
    if !all([col in propertynames(layout) for col in [:inc, :azi]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end

    # Validate data types
    if !(eltype(layout.inc) <: Number && eltype(layout.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    radius = 88 # mm
    inc = layout[!, :inc] .* (pi / 180)
    azi = layout[!, :azi] .* (pi / 180)

    layout[!, :x2] = inc .* cos.(azi) .* radius
    layout[!, :y2] = inc .* sin.(azi) .* radius

    return nothing
end

"""
    polar_to_cartesian_xyz!(layout::DataFrame)

Converts polar coordinates (incidence and azimuth angles) from a layout DataFrame into Cartesian coordinates (x, y, z).

# Arguments
- `layout::DataFrame`: A DataFrame containing the layout information with columns for incidence angles (`:inc`) and azimuth angles (`:azi`).

# Modifies
- The input `layout` DataFrame is modified in place to include new columns `x3`, `y3`, and `z3`, which represent the Cartesian coordinates calculated from the polar coordinates.

# Returns
- Nothing. The function modifies the `layout` DataFrame directly.
"""
function polar_to_cartesian_xyz!(layout::DataFrame)
    # Check for required columns
    if !all([col in propertynames(layout) for col in [:inc, :azi]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end

    # Validate data types
    if !(eltype(layout.inc) <: Number && eltype(layout.azi) <: Number)
        throw(ArgumentError(":inc and :azi columns must contain numeric values"))
    end

    radius = 88.0  # mm
    inc = layout[!, :inc] .* (pi / 180)  # Convert to radians
    azi = layout[!, :azi] .* (pi / 180)  # Convert to radians

    layout[!, :x3] = radius .* sin.(inc) .* cos.(azi)
    layout[!, :y3] = radius .* sin.(inc) .* sin.(azi)
    layout[!, :z3] = radius .* cos.(inc)

    return nothing
end

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




"""
    get_electrode_neighbours_xy(layout::DataFrame, distance_criterion::Real)

Identifies the neighbours of each electrode based on their Cartesian coordinates.

# Arguments
- `layout::DataFrame`: A DataFrame containing the layout information with columns for electrode labels, Cartesian coordinates (`x2`, `y2`).
- `distance_criterion::Real`: The maximum distance to consider two electrodes as neighbours.

# Returns
- `OrderedDict{Symbol,Neighbours}`: A dictionary where each key is an electrode label, and the value is a Neighbours struct containing neighbour information.

# Throws
- `ArgumentError`: If the layout DataFrame does not contain the required columns.
"""
function get_electrode_neighbours_xy(layout::DataFrame, distance_criterion::Real)

    if !all([col in propertynames(layout) for col in [:x2, :y2, :label]])
        throw(ArgumentError("Layout must contain x2, y2, and :label columns"))
    end

    if distance_criterion <= 0
        throw(ArgumentError("Distance criterion must be positive"))
    end

    # Precompute coordinates
    coords = Matrix{Float64}(undef, size(layout, 1), 2)
    coords[:, 1] = layout.x2
    coords[:, 2] = layout.y2

    neighbour_dict = OrderedDict{Symbol,Neighbours}()

    for (idx1, label1) in enumerate(layout.label)

        neighbour_dict[Symbol(label1)] = Neighbours([], [], [])

        for (idx2, label2) in enumerate(layout.label)
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

    return neighbour_dict

end

"""
    get_electrode_neighbours_xyz(layout::DataFrame, distance_criterion::Real)

Identifies the neighbours of each electrode based on their Cartesian coordinates.

# Arguments
- `layout::DataFrame`: A DataFrame containing the layout information with columns for electrode labels, Cartesian coordinates (`x3`, `y3`, `z3`).
- `distance_criterion::Real`: The maximum distance to consider two electrodes as neighbours.

# Returns
- `OrderedDict{Symbol,Neighbours}`: A dictionary where each key is an electrode label, and the value is a Neighbours struct containing neighbour information.

# Throws
- `ArgumentError`: If the layout DataFrame does not contain the required columns.
"""
function get_electrode_neighbours_xyz(layout::DataFrame, distance_criterion::Real)

    if !all([col in propertynames(layout) for col in [:x3, :y3, :z3, :label]])
        throw(ArgumentError("Layout must contain :x3, :y3, :z3, and :label columns"))
    end

    if distance_criterion <= 0
        throw(ArgumentError("Distance criterion must be positive"))
    end

    # Precompute coordinates
    coords = Matrix{Float64}(undef, size(layout, 1), 3)
    coords[:, 1] = layout.x3
    coords[:, 2] = layout.y3
    coords[:, 3] = layout.z3

    neighbour_dict = OrderedDict{Symbol, Neighbours}()

    for (idx1, label1) in enumerate(layout.label)

        neighbour_dict[Symbol(label1)] = Neighbours([], [], [])

        for (idx2, label2) in enumerate(layout.label)

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

    return neighbour_dict

end


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
            "neighbors" => [string(n) for n in neighbours.electrodes],
            "distances" => [round(d, digits=4) for d in neighbours.distances],
            "weights" => [round(w, digits=6) for w in neighbours.weights],
            "neighbor_count" => length(neighbours.electrodes)
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
neighbours = get_electrode_neighbours_xy(layout, 40)

# Write to TOML file
print_neighbours_dict(neighbours, "neighbours.toml")
```
"""
function print_neighbours_dict(neighbours_dict::OrderedDict{Symbol, Neighbours}, filename::String)
    nneighbours = average_number_of_neighbours(neighbours_dict)
    toml_data = _format_neighbours_toml(neighbours_dict, nneighbours)
    open(filename, "w") do io
        TOML.print(io, toml_data)
    end
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
neighbours, _ = get_electrode_neighbours_xy(layout, 40)
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


