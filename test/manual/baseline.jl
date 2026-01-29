using EegFun

# read raw data
dat = EegFun.read_raw_data("./data/raw_files/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Baseline stuff
EegFun.baseline!(dat)

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 1));
dat.data

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 2));
dat.data

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 10));
dat.data

EegFun.baseline!(dat, EegFun.IntervalTime(start = 0.0, stop = 0.0));
dat.data

EegFun.baseline!(dat, EegFun.IntervalTime(start = 0.0, stop = 0.004));
dat.data
