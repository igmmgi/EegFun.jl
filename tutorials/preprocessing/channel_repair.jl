# Demo: Channel Repair
# Shows neighbor interpolation and spherical spline for repairing bad channels.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

dat = EegFun.create_eegfun_data(dat, layout)

# minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

EegFun.plot_databrowser(dat)

# Select a channel to repair and make this channel noisy!
channel_to_repair = :Cz
dat.data[!, channel_to_repair] .+= randn(size(dat.data[:, channel_to_repair])) * 200 # v. noisy!

# We can now see this noisy channel in the databrowser
# NB. we can actually press "R" and select Cz and apply the repair in the browser
EegFun.plot_databrowser(dat)

# Try neighbor interpolation
EegFun.repair_channels!(dat, [channel_to_repair], method = :neighbor_interpolation)

# Cz is now repaired
EegFun.plot_databrowser(dat)

# Try neighbor interpolation
EegFun.repair_channels!(dat, [channel_to_repair], method = :spherical_spline)

# Cz is now repaired
EegFun.plot_databrowser(dat)


