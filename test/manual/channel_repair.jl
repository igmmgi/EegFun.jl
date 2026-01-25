using EegFun
using GLMakie

# Load some test data
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)
EegFun.polar_to_cartesian_xyz!(layout_file)

dat = EegFun.read_bdf(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)
EegFun.rereference!(dat, :avg)
EegFun.filter_data!(dat, "hp", 1)

# Select a few channels to repair
test_channels = [:Fp1]
available_channels = dat.layout.data.label
channels_to_repair = filter(ch -> ch in available_channels, test_channels)

# Calculate neighbors first
EegFun.get_layout_neighbours_xyz!(dat.layout, 0.5)

# Store original data for comparison
original_data = copy(dat.data)

# Try neighbor interpolation
EegFun.repair_channels!(dat, channels_to_repair, method = :neighbor_interpolation)

# Try spherical spline
EegFun.repair_channels!(dat, channels_to_repair, method = :spherical_spline)

# Check if data changed (using isapprox to handle floating point precision)
data_changed = any(any(.!isapprox.(dat.data[!, ch], original_data[!, ch], rtol = 1e-10)) for ch in channels_to_repair)
