using EegFun
using GLMakie
using BenchmarkTools
# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..","..", "AttentionExp", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);
EegFun.rereference!(dat, :mastoid)
EegFun.filter_data!(dat, "hp", 1)

# Correlation Matrix
cm = EegFun.correlation_matrix(dat)
fig, ax = EegFun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")

cm = EegFun.correlation_matrix(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :C3, :C4]))
fig, ax = EegFun.plot_correlation_heatmap(cm, title = "Selected Correlations")

# More custom figure
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
cm = EegFun.correlation_matrix(dat)
EegFun.plot_correlation_heatmap!(fig, ax1, cm, title = "Full Correlation Matrix")
cm = EegFun.correlation_matrix(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :C3, :C4]))
EegFun.plot_correlation_heatmap!(fig, ax2, cm, title = "Selected Correlation Matrix", colorbar_position = (2, 2))
fig


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
EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# Calculate correlations between all channels and EOG channels
cm = EegFun.correlation_matrix_dual_selection(
    dat,
    sample_selection = EegFun.samples(),  # All samples
    channel_selection1 = EegFun.channels(),  # All EEG channels
    channel_selection2 = EegFun.channels([:vEOG, :hEOG]),  # EOG channels
)
EegFun.add_zscore_columns!(cm)


cm = EegFun.correlation_matrix(dat)
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

EegFun.plot_layout_2d(layout_file, correlation_matrix = cm)
