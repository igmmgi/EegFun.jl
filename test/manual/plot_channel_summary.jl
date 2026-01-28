using EegFun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

cs = EegFun.channel_summary(dat)

# # basic channel summary statistics
# EegFun.channel_summary(dat)
# EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
# EegFun.channel_summary(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]), sample_selection = x -> x.sample .< 2000)
# EegFun.channel_summary(dat, channel_selection = x -> endswith.(string.(x), "z")) # all midline channels 
# EegFun.channel_summary(dat, channel_selection = x -> .!(endswith.(string.(x), "z"))) # all non-midline channels 
# EegFun.channel_summary(dat, include_extra = true) # include additional channels (e.g. vEOG, hEOG)


EegFun.plot_channel_summary(cs, :range)
EegFun.plot_channel_summary(cs, :min, grid_visible = false)
EegFun.plot_channel_summary(cs, :min, bar_color = :red)

EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar])
EegFun.plot_channel_summary(cs, [:min, :max, :std], dims = (1, 3), grid_visible = false)

#################################
# Epoched DataFrameEeg
#################################
# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

cs = EegFun.channel_summary(epochs[1])
EegFun.plot_channel_summary(cs, :range, average_over = :epoch)

EegFun.plot_channel_summary(cs, [:min, :max, :std, :range, :var, :zvar], average_over = :epoch)
EegFun.plot_channel_summary(cs, [:min, :max, :std], dims = (1, 3), grid_visible = false, average_over = :epoch)
