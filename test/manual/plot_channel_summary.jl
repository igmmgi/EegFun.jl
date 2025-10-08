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
eegfun.plot_channel_summary(cs, :min, grid_visible = false)
eegfun.plot_channel_summary(cs, :min, bar_color = :red)

eegfun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar])
eegfun.plot_channel_summary(cs, [:min, :max, :std], dims = (1, 3), grid_visible = false)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [ eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

cs = eegfun.channel_summary(epochs[1])
eegfun.plot_channel_summary(cs, :range, average_over = :epoch)

eegfun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)
eegfun.plot_channel_summary(cs, [:min, :max, :std], dims = (1, 3), grid_visible = false, average_over = :epoch)


