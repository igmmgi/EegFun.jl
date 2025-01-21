mutable struct CoordXY
    x::Any
    y::Any
end

mutable struct Coord
    coord::Array{CoordXY}
end

read_layout(layout_file_name) = DataFrame(CSV.File(layout_file_name))


function polar_to_cartesian_xy!(layout::DataFrame)
    radius = 88 # mm
    inc = layout[!, :inc] .* (pi / 180)
    azi = layout[!, :azi] .* (pi / 180)
    layout[!, "x2"] = inc .* cos.(azi) .* radius
    layout[!, "y2"] = inc .* sin.(azi) .* radius
    return
end

function polar_to_cartesian_xyz!(layout::DataFrame)

    radius = 88 # mm
    inc = layout[!, :inc] .* (pi / 180)
    azi = layout[!, :azi] .* (pi / 180)
    layout[!, "x3"] = radius .* sin.(inc) .* cos.(azi)
    layout[!, "y3"] = radius .* sin.(inc) .* sin.(azi)
    layout[!, "z3"] = radius .* cos.(inc)
    return
end

# distances
calculate_distance_xy(x1, y1, x2, y2) =  sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))
calculate_distance_xyz(x1, y1, z1, x2, y2, z2) = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)) + ((z1 - z2) * (z1 - z2)))

function get_electrode_neighbours_xy(layout, distance_criterion)
    neighbour_dict = OrderedDict()
    for (idx_electrode1, label_electrode1) in enumerate(layout.label)
        neighbour_dict[Symbol(label_electrode1)] = []
        for (idx_electrode2, label_electrode2) in enumerate(layout.label)
            if (idx_electrode1 == idx_electrode2)
                continue
            end
            distance = calculate_distance_xy(
                layout.x2[idx_electrode1],
                layout.y2[idx_electrode1],
                layout.x2[idx_electrode2],
                layout.y2[idx_electrode2],
            )
            if distance <= distance_criterion
                push!(neighbour_dict[Symbol(label_electrode1)], Symbol(label_electrode2))
            end
        end
    end
    return neighbour_dict
end

function get_electrode_neighbours_xyz(layout, distance_criterion)
    neighbour_dict = OrderedDict()
    for (idx_electrode1, label_electrode1) in enumerate(layout.label)
        neighbour_dict[Symbol(label_electrode1)] = []
        for (idx_electrode2, label_electrode2) in enumerate(layout.label)
            if (idx_electrode1 == idx_electrode2)
                continue
            end
            distance = calculate_distance_xyz(
                layout.x3[idx_electrode1],
                layout.y3[idx_electrode1],
                layout.z3[idx_electrode1],
                layout.x2[idx_electrode2],
                layout.y2[idx_electrode2],
                layout.z3[idx_electrode2],
            )
            if distance <= distance_criterion
                push!(neighbour_dict[Symbol(label_electrode1)], Symbol(label_electrode2))
            end
        end
    end
    return neighbour_dict
end
