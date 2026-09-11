# Demo: Correlation Heatmap
# Shows correlation matrices between channels or trials.

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

# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# Correlation Matrix
cm = EegFun.correlation_matrix(dat)
EegFun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")

cm = EegFun.correlation_matrix(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :C3, :C4]))
EegFun.plot_correlation_heatmap(cm, title = "Selected Correlations")

# Calculate EOG signals
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # horizontal EOG = F9 - F10

# Calculate correlations between all channels and EOG channels
cm_eog = EegFun.correlation_matrix_dual_selection(
    dat,
    sample_selection = EegFun.samples(),  # All samples
    channel_selection1 = EegFun.channels(),  # All EEG channels
    channel_selection2 = EegFun.channels([:vEOG, :hEOG]),  # EOG channels
)

# Visualize which EEG channels are correlated with EOG
EegFun.plot_correlation_heatmap(cm_eog, title = "EEG-EOG Correlations", xlabel = "EOG Channels", ylabel = "EEG Channels")

# Display channel correlations in layout
cm = EegFun.correlation_matrix(dat)
EegFun.plot_layout_2d(layout, correlation_matrix = cm)
