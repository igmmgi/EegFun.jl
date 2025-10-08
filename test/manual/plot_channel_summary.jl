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

cs = eegfun.channel_summary(dat)
eegfun.plot_channel_summary(cs, :range)

# TODO: make this a 
fig = Figure()
ax1 = Axis(fig[1, 1])
eegfun.plot_channel_summary!(fig, ax1, cs, :min)
ax2 = Axis(fig[1, 2])
eegfun.plot_channel_summary!(fig, ax2, cs, :max)
ax3 = Axis(fig[1, 3])
eegfun.plot_channel_summary!(fig, ax3, cs, :std)
ax4 = Axis(fig[2, 1])
eegfun.plot_channel_summary!(fig, ax4, cs, :range)
ax5 = Axis(fig[2, 2])
eegfun.plot_channel_summary!(fig, ax5, cs, :var)
ax6 = Axis(fig[2, 3])
eegfun.plot_channel_summary!(fig, ax6, cs, :zvar)
fig

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [ eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

cs = eegfun.channel_summary(epochs[1])
eegfun.plot_channel_summary(cs, :range, average_over = :epoch)

fig = Figure()
ax1 = Axis(fig[1, 1])
eegfun.plot_channel_summary!(fig, ax1, cs, :min, average_over = :epoch)
ax2 = Axis(fig[1, 2])
eegfun.plot_channel_summary!(fig, ax2, cs, :max, average_over = :epoch)
ax3 = Axis(fig[1, 3])
eegfun.plot_channel_summary!(fig, ax3, cs, :std, average_over = :epoch)
ax4 = Axis(fig[2, 1])
eegfun.plot_channel_summary!(fig, ax4, cs, :range, average_over = :epoch)
ax5 = Axis(fig[2, 2])
eegfun.plot_channel_summary!(fig, ax5, cs, :var, average_over = :epoch)
ax6 = Axis(fig[2, 3])
eegfun.plot_channel_summary!(fig, ax6, cs, :zvar, average_over = :epoch)
fig






