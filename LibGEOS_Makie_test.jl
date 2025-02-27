using Meshes

points = [
    (0.5, 0.0),    # Center
    (1.0, 0.0),    # Right
    (1.0, 1.0),    # Top right
    (0.0, 1.0),    # Top
]

points2 = [Meshes.Point(p...) for p in points]
chul = convexhull(points2)
coords=[[vert.coords.x.val,vert.coords.y.val] for vert in chul.rings[1].vertices]
fig = Figure()
ax = Axis(fig[1,1])

# draw the points
scatter!(ax, Point2f.(points))

# draw the convex hull
push!(coords,coords[1])
lines!(ax, Point2f.(coords))
display(fig)
