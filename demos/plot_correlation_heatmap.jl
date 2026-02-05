using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

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
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)

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
EegFun.plot_layout_2d(layout_file, correlation_matrix = cm)
