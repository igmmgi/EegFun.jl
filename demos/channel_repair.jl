using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)

EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Select a few channels to repair
test_channels = [:Fp1]
available_channels = dat.layout.data.label
channels_to_repair = filter(ch -> ch in available_channels, test_channels)

# Calculate neighbors first
EegFun.get_neighbours_xyz!(dat.layout, 0.5)

# Store original data for comparison
original_data = copy(dat.data)

# Try neighbor interpolation
EegFun.repair_channels!(dat, channels_to_repair, method = :neighbor_interpolation)

# Try spherical spline
EegFun.repair_channels!(dat, channels_to_repair, method = :spherical_spline)

# Check if data changed (using isapprox to handle floating point precision)
data_changed = any(any(.!isapprox.(dat.data[!, ch], original_data[!, ch], rtol = 1e-10)) for ch in channels_to_repair)
