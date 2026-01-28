# package
using EegFun
using GLMakie

dat = EegFun.read_raw_data("../Flank_C_3.bdf");
layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = EegFun.create_eeg_dataframe(dat, layout);
EegFun.highpass_filter!(dat, 1)
EegFun.rereference!(dat, :avg)

# test plot_channel_summary
jp = EegFun.channel_joint_probability(dat)
EegFun.plot_joint_probability(jp)
