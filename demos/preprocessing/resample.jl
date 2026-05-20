# Demo: Resampling
# Shows how to resample data to different sampling rates.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

EegFun.sample_rate(dat)      # current sample rate
EegFun.trigger_count(dat)    # current triggers in file
EegFun.plot_databrowser(dat) # view current data

dat_new = EegFun.resample(dat, 2) # downsample by a factor of 2
EegFun.sample_rate(dat_new)       # should = original ÷ 2
EegFun.trigger_count(dat_new)     # triggers should be preserved
EegFun.plot_databrowser(dat_new)  # view current data


dat_new = EegFun.resample(dat, 4) # downsample by a factor of 4
EegFun.sample_rate(dat_new)       # should = original ÷ 4
EegFun.trigger_count(dat_new)     # triggers should be preserved
EegFun.plot_databrowser(dat_new)  # view current data
