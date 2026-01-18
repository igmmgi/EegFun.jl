using eegfun
using GLMakie

# Load some test data
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv")
eegfun.polar_to_cartesian_xy!(layout_file)
eegfun.polar_to_cartesian_xyz!(layout_file)

dat = eegfun.read_bdf(data_file)
dat = eegfun.create_eeg_dataframe(dat, layout_file)
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# Select a few channels to repair
test_channels = [:Fp1]
available_channels = dat.layout.data.label
channels_to_repair = filter(ch -> ch in available_channels, test_channels)

# Calculate neighbors first
eegfun.get_layout_neighbours_xyz!(dat.layout, 0.5)

# Store original data for comparison
original_data = copy(dat.data)

# Try neighbor interpolation
eegfun.repair_channels!(dat, channels_to_repair, method = :neighbor_interpolation)

# Try spherical spline
eegfun.repair_channels!(dat, channels_to_repair, method = :spherical_spline)

# Check if data changed (using isapprox to handle floating point precision)
data_changed = any(any(.!isapprox.(dat.data[!, ch], original_data[!, ch], rtol = 1e-10)) for ch in channels_to_repair)
