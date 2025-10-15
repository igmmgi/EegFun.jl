fig = Figure()
Axis(fig[1, 1])
# vertical colorbars
Colorbar(fig[1, 2], limits = (0, 10), colormap = :viridis, tellheight = true, tellwidth = true)
Axis(fig[2, 1])
Axis(fig[3, 1])
Axis(fig[1, 3])
Axis(fig[2, 3])
Axis(fig[3, 3])
Colorbar(fig[3, 4], limits = (0, 10), colormap = :viridis, tellheight = true, tellwidth = true)
fig
