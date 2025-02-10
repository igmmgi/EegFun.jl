
"""
    read_layout(layout_file_name::String)

Reads a layout file in CSV format and returns a DataFrame containing the layout information.

# Arguments
- `layout_file_name::String`: The path to the CSV file containing the layout data.

# Returns
- `DataFrame`: A DataFrame containing the layout data, with columns for electrode labels, incidence angles, azimuth angles, and calculated Cartesian coordinates (if applicable).
"""
read_layout(layout_file_name) = DataFrame(CSV.File(layout_file_name))::DataFrame

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
    # Validate input columns
    if !all([col in names(layout) for col in ["inc", "azi"]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end
    radius = 88 # mm
    inc = layout[!, :inc] .* (pi / 180)
    azi = layout[!, :azi] .* (pi / 180)
    layout[!, "x2"] = inc .* cos.(azi) .* radius
    layout[!, "y2"] = inc .* sin.(azi) .* radius
    return
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
    # Validate input columns
    if !all([col in names(layout) for col in ["inc", "azi"]])
        throw(ArgumentError("Layout must contain :inc and :azi columns"))
    end

    radius = 88.0  # mm
    inc = layout[!, :inc] .* (pi / 180)  # Convert to radians
    azi = layout[!, :azi] .* (pi / 180)  # Convert to radians

    layout[!, :x3] = radius .* sin.(inc) .* cos.(azi)
    layout[!, :y3] = radius .* sin.(inc) .* sin.(azi)
    layout[!, :z3] = radius .* cos.(inc)
end

# distances
"""
    calculate_distance_xy(x1::Real, y1::Real, x2::Real, y2::Real)

Calculates the Euclidean distance between two points in 2D space.

# Arguments
- `x1, y1, x2, y2`: The coordinates of the two points.

# Returns
- `Float64`: The distance between the two points.
"""
function calculate_distance_xy(x1::Real, y1::Real, x2::Real, y2::Real)::Float64
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

"""
    calculate_distance_xyz(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)

Calculates the Euclidean distance between two points in 3D space.

# Arguments
- `x1, y1, z1, x2, y2, z2`: The coordinates of the two points.

# Returns
- `Float64`: The distance between the two points.
"""
function calculate_distance_xyz(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)::Float64
    return sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
end

"""
    get_electrode_neighbours_xy(layout::DataFrame, distance_criterion::Real)

Identifies the neighbours of each electrode based on their Cartesian coordinates.

# Arguments
- `layout::DataFrame`: A DataFrame containing the layout information with columns for electrode labels, Cartesian coordinates (`x2`, `y2`).
- `distance_criterion::Real`: The maximum distance to consider two electrodes as neighbours.

# Returns
- `OrderedDict{Symbol,Vector{Symbol}}`: A dictionary where each key is an electrode label, and the value is a list of labels of its neighbours.

# Throws
- `ArgumentError`: If the layout DataFrame does not contain the required columns.
"""
function get_electrode_neighbours_xy(layout::DataFrame, distance_criterion::Real)
    if !all(in.([:x2, :y2, :label], names(layout)))
        throw(ArgumentError("Layout must contain :x2, :y2, and :label columns"))
    end

    if distance_criterion <= 0
        throw(ArgumentError("Distance criterion must be positive"))
    end

    neighbour_dict = OrderedDict{Symbol,Vector{Symbol}}()

    for (idx1, label1) in enumerate(layout.label)
        neighbour_dict[Symbol(label1)] = Symbol[]
        for (idx2, label2) in enumerate(layout.label)
            if idx1 == idx2
                continue
            end

            distance = calculate_distance_xy(layout.x2[idx1], layout.y2[idx1], layout.x2[idx2], layout.y2[idx2])

            if distance <= distance_criterion
                push!(neighbour_dict[Symbol(label1)], Symbol(label2))
            end
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
- `OrderedDict{Symbol,Vector{Symbol}}`: A dictionary where each key is an electrode label, and the value is a list of labels of its neighbours.

# Throws
- `ArgumentError`: If the layout DataFrame does not contain the required columns.
"""
function get_electrode_neighbours_xyz(layout::DataFrame, distance_criterion::Real)
    if !all(in.([:x3, :y3, :z3, :label], names(layout)))
        throw(ArgumentError("Layout must contain :x3, :y3, :z3, and :label columns"))
    end

    if distance_criterion <= 0
        throw(ArgumentError("Distance criterion must be positive"))
    end

    neighbour_dict = OrderedDict{Symbol,Vector{Symbol}}()

    for (idx1, label1) in enumerate(layout.label)
        neighbour_dict[Symbol(label1)] = Symbol[]
        for (idx2, label2) in enumerate(layout.label)
            if idx1 == idx2
                continue
            end

            distance = calculate_distance_xyz(
                layout.x3[idx1],
                layout.y3[idx1],
                layout.z3[idx1],
                layout.x3[idx2],
                layout.y3[idx2],
                layout.z3[idx2],
            )

            if distance <= distance_criterion
                push!(neighbour_dict[Symbol(label1)], Symbol(label2))
            end
        end
    end

    return neighbour_dict
end
