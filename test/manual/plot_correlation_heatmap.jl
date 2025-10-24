using eegfun
using GLMakie
using BenchmarkTools
# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

# Correlation Matrix
cm = eegfun.correlation_matrix(dat)
fig, ax = eegfun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")

cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :C3, :C4]))
fig, ax = eegfun.plot_correlation_heatmap(cm, title = "Selected Correlations")

# More custom figure
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
cm = eegfun.correlation_matrix(dat)
eegfun.plot_correlation_heatmap!(fig, ax1, cm, title = "Full Correlation Matrix")
cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :C3, :C4]))
eegfun.plot_correlation_heatmap!(fig, ax2, cm, title = "Selected Correlation Matrix", colorbar_position = (2, 2))
fig


# Calculate EOG signals
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:Fp1, :Fp2]),
    channel_selection2 = eegfun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.channel_difference!(
    dat,
    channel_selection1 = eegfun.channels([:F9]),
    channel_selection2 = eegfun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
eegfun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
eegfun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)



# Calculate correlations between all channels and EOG channels
cm = eegfun.correlation_matrix_dual_selection(dat, 
    sample_selection = eegfun.samples(),  # All samples
    channel_selection1 = eegfun.channels(),  # All EEG channels
    channel_selection2 = eegfun.channels([:vEOG, :hEOG]),  # EOG channels
)
eegfun.add_zscore_columns!(cm)
