using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file)

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


