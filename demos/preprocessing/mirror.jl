# Demo: Mirror Padding
# Shows how to pad epoch and ERP data with mirrored edges (:pre, :post, :both)
# to reduce edge artifacts when filtering.

using EegFun
# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

# read raw data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"));

# read and prepare layout file
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"));
EegFun.polar_to_cartesian_xy!(layout)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Create some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :pre)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :post)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :both)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :pre)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :post)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :both)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))


