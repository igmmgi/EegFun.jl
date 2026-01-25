using EegFun
using GLMakie

# Load some test data
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)
EegFun.polar_to_cartesian_xyz!(layout_file)

dat = EegFun.read_bdf(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

summary = EegFun.channel_summary(dat)
EegFun.log_pretty_table(summary; title = "Initial Channel Summary")

EegFun.is_extreme_value!(dat, 100);
summary = EegFun.channel_summary(dat, sample_selection = EegFun.samples_not(:is_extreme_value_100))
