# package
using EegFun
using GLMakie

dat = EegFun.read_bdf("../Flank_C_3.bdf");
layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = EegFun.create_eeg_dataframe(dat, layout);
EegFun.filter_data!(dat, "hp", 1)
EegFun.rereference!(dat, :avg)

# test plot_channel_summary
jp = EegFun.channel_joint_probability(dat)
EegFun.plot_joint_probability(jp)
