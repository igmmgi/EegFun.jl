using LibGEOS, GLMakie

# Create a simple polygon using LibGEOS
coords = [
    (0.0, 0.0),    # Center
    (1.0, 0.0),    # Right
    (1.0, 1.0),    # Top right
    (0.0, 1.0),    # Top
    (0.0, 0.0),    # Back to start
]

# Convert to LibGEOS geometry
shell = map(x -> LibGEOS.Coordinate(x[1], x[2]), coords)
geom = LibGEOS.Polygon(shell)

# Create a buffer around the polygon
buffered = buffer(geom, 0.2)

# Extract coordinates for plotting
coords_out = coordinates(buffered)
x = [c[1] for c in coords_out]
y = [c[2] for c in coords_out]

# Plot using Makie
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, x, y)
scatter!(ax, x, y)

# Make axes equal and display
ax.aspect = DataAspect()
display(fig)