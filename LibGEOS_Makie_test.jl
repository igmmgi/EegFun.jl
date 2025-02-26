using GLMakie

# Create points in a + configuration
points = [
    Point2f(0.5, 0.3),  # Bottom of vertical line
    Point2f(0.5, 0.7),  # Top of vertical line
    Point2f(0.3, 0.5),  # Left of horizontal line
    Point2f(0.7, 0.5),  # Right of horizontal line
    Point2f(0.5, 0.5)   # Center point
]

xs = [p[1] for p in points]
ys = [p[2] for p in points]

function graham_scan(points::Vector{Point2f})
    # Need at least 3 points for a hull
    length(points) < 3 && return points
    
    # Find leftmost point
    start = 1
    for i in 2:length(points)
        if points[i][1] < points[start][1]
            start = i
        end
    end
    p0 = points[start]
    
    # Sort by angle from p0
    other_points = vcat(points[1:start-1], points[start+1:end])
    sorted = sort(other_points, 
        by = p -> atan(p[2] - p0[2], p[1] - p0[1]))
    
    # Initialize hull with first point
    hull = [p0]
    
    # Build hull
    for p in sorted
        while length(hull) > 1
            v1 = hull[end] - hull[end-1]
            v2 = p - hull[end]
            cross = v1[1] * v2[2] - v1[2] * v2[1]
            if cross > 0  # Left turn
                break
            end
            pop!(hull)
        end
        push!(hull, p)
    end
    
    return hull
end

function knn_concave_hull(points::Vector{Point2f}, k::Int)
    # Start with convex hull
    hull = graham_scan(points)
    push!(hull, hull[1])  # Close the polygon
    
    # Iteratively refine the hull
    max_iterations = 100
    iteration = 0
    while iteration < max_iterations
        # Find point pairs that are too far apart
        max_dist = 0.0
        max_idx = 0
        for i in 1:length(hull)-1
            dist = sqrt((hull[i+1][1] - hull[i][1])^2 + (hull[i+1][2] - hull[i][2])^2)
            if dist > max_dist
                max_dist = dist
                max_idx = i
            end
        end
        
        # If all distances are acceptable, stop
        if max_dist < 0.3
            break
        end
        
        # Find nearest point to midpoint of longest edge
        mid_x = (hull[max_idx+1][1] + hull[max_idx][1]) / 2
        mid_y = (hull[max_idx+1][2] + hull[max_idx][2]) / 2
        
        min_dist = Inf
        nearest_point = nothing
        for p in points
            if p ∉ hull
                dist = sqrt((p[1] - mid_x)^2 + (p[2] - mid_y)^2)
                if dist < min_dist
                    min_dist = dist
                    nearest_point = p
                end
            end
        end
        
        # If found a point, insert it
        if nearest_point !== nothing
            insert!(hull, max_idx + 1, nearest_point)
        end
        
        iteration += 1
    end
    
    return hull
end

# Create circle points around each point
border_size = 0.05
circle_points = 0:2π/36:2π
border_points = Point2f[]
for p in points
    for θ in circle_points
        push!(border_points, Point2f(
            p[1] + border_size * sin(θ),
            p[2] + border_size * cos(θ)
        ))
    end
end

# Plot both hulls
fig = Figure(resolution=(800, 400))

# Convex Hull (left plot)
ax1 = Axis(fig[1, 1], title="Convex Hull")
convex_hull = graham_scan(border_points)
push!(convex_hull, convex_hull[1])  # Close the polygon
convex_xs = [p[1] for p in convex_hull]
convex_ys = [p[2] for p in convex_hull]
poly!(ax1, convex_xs, convex_ys, color=(:blue, 0.2))
lines!(ax1, convex_xs, convex_ys, color=:blue)
scatter!(ax1, xs, ys, color=:red, markersize=8)
xlims!(ax1, 0.0, 1.0)
ylims!(ax1, 0.0, 1.0)

# Concave Hull (right plot)
ax2 = Axis(fig[1, 2], title="Concave Hull")
concave_hull = knn_concave_hull(border_points, 3)
concave_xs = [p[1] for p in concave_hull]
concave_ys = [p[2] for p in concave_hull]
poly!(ax2, concave_xs, concave_ys, color=(:blue, 0.2))
lines!(ax2, concave_xs, concave_ys, color=:blue)
scatter!(ax2, xs, ys, color=:red, markersize=8)
xlims!(ax2, 0.0, 1.0)
ylims!(ax2, 0.0, 1.0)

display(fig)