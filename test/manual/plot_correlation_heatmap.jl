using eegfun
using GLMakie
using BenchmarkTools
# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
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

# M;ore custom figure
fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
cm = eegfun.correlation_matrix(dat)
eegfun.plot_correlation_heatmap!(fig, ax1, cm, title = "Full Correlation Matrix")
cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :C3, :C4]))
eegfun.plot_correlation_heatmap!(fig, ax2, cm, title = "Selected Correlation Matrix", colorbar_position = (2, 2))
fig





