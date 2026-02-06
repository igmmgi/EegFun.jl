using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eegfun_data(dat, layout_file)

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
