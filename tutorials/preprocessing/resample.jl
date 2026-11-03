# Demo: Resampling
# Shows how to resample EEG data to a target sampling rate (Hz).

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

EegFun.sample_rate(dat)      # current sample rate (Hz)
EegFun.trigger_count(dat)    # current triggers in file
EegFun.plot_databrowser(dat) # view current data

dat_512 = EegFun.resample(dat, 512) # resample to 512 Hz (polyphase anti-aliasing)
EegFun.sample_rate(dat_512)         # 512 Hz
EegFun.trigger_count(dat_512)       # triggers are preserved
EegFun.plot_databrowser(dat_512)    # view resampled data

dat_256 = EegFun.resample(dat, 256) # resample to 256 Hz
EegFun.sample_rate(dat_256)         # 256 Hz
EegFun.trigger_count(dat_256)       # triggers are preserved
EegFun.plot_databrowser(dat_256)    # view resampled data
