# Demo: Joint Probability Plot
# Shows joint probability distributions for artifact detection.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# test plot_channel_summary
jp = EegFun.channel_joint_probability(dat)
EegFun.plot_joint_probability(jp)
