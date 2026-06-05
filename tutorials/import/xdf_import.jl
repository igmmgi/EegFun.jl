# Demo: Loading and Processing XDF Files
# Shows how to load Extensible Data Format (.xdf) files,
# create EegFun data structures, apply basic preprocessing,
# and visualize the data.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.xdf")

using EegFun

# Load raw XDF data
raw_xdf = EegFun.read_raw_data(EegFun.example_path("data/xdf/test1.xdf"))

# Using the fallback auto-generated layout
dat = EegFun.create_eegfun_data(raw_xdf)

# Basic preprocessing
EegFun.highpass_filter!(dat, 1.0) # 1 Hz high-pass filter
EegFun.rereference!(dat, :avg)    # Average reference

println("\nTrigger summary:")
EegFun.trigger_count(dat)

# Open interactive databrowser
# EegFun.plot_databrowser(dat);
