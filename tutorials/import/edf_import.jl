# Demo: Loading and Processing EDF Files
# Shows how to load European Data Format (.edf) files,
# create EegFun data structures, apply basic preprocessing,
# and visualize the data.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.edf")

using EegFun
using GLMakie

# Load raw EDF data
raw_edf = EegFun.read_raw_data(EegFun.example_path("data/edf/test.edf"))

# Using the fallback auto-generated layout (if you don't have a specific .csv layout)
dat = EegFun.create_eegfun_data(raw_edf)

EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.all_labels(dat)
EegFun.channel_labels(dat)
EegFun.meta_labels(dat)
EegFun.extra_labels(dat) # empty


# Basic preprocessing
EegFun.highpass_filter!(dat, 1.0) # 1 Hz high-pass filter
EegFun.rereference!(dat, :avg)    # Average reference

println("\nTrigger summary:")
EegFun.trigger_count(dat)

# Open interactive databrowser
EegFun.plot_databrowser(dat);
