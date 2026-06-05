# Demo: Channel Summary
# Shows how to create summary statistics for channels across epochs.

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
EegFun.highpass_filter!(dat, 0.1)

# summary statistics across all channels
summary = EegFun.channel_summary(dat)
EegFun.log_pretty_table(summary; title = "Initial Channel Summary")

# summary statistics across all channels excluding v. extreme values
EegFun.is_extreme_value!(dat, 200);
summary = EegFun.channel_summary(dat, sample_selection = EegFun.samples_not(:is_extreme_value_200))
EegFun.log_pretty_table(summary; title = "Channel Summary (excluding extreme values)")

# summary statistics across all Midline channels via predicate selection
summary = EegFun.channel_summary(dat, channel_selection = EegFun.channels(x -> endswith.(string.(x), "z")))
EegFun.log_pretty_table(summary; title = "Channel Summary (Midline)")
