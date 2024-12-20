###########################################################
mapcols(x -> filtfilt(filter, x), dat.data[:, 3:end])
mapcols!(x -> filtfilt(filter, x), dat.data[:, 3:end])

###########################################################
using CoordinateTransformations

layout = read_layout("/home/ian/Documents/Julia/EEGfun/layouts/biosemi72.csv")
polar_to_cartesian_xy!(layout)
polar_to_cartesian_xyz!(layout)

s = Spherical.(88, layout.inc, layout.azi)
c = CartesianFromSpherical().(s)
scatter(c)


###########################################################
# potential start point for 3D head topo
brain = load(assetpath("/home/ian/Downloads/OBJ/Super Average Head.obj"))
f = Figure()
ax = Axis3(f[1, 1])
hidedecorations!(ax)  # hides ticks, grid and lables
hidespines!(ax)  # hid
clip_planes = [Plane3f(Point3f(0), Vec3f(0, 1, 0))]
mesh!(ax, brain, color = :grey, clip_planes = clip_planes)
wireframe!(ax, brain, clip_planes = clip_planes, color = :grey)



# is = IntervalSlider(fig[2, 1], range=@lift(data[$xrange, :time]), startvalues=(1, 2))
# vl = lift(x -> x[1], is.interval)
# vr = lift(x -> x[2], is.interval)
# vspan!(fig[1, 1], vl, vr)
