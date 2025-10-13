# package
using eegfun
using GLMakie

dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);
eegfun.filter_data!(dat, "hp", 1)
eegfun.rereference!(dat, :avg)

# test plot_channel_summary
jp = eegfun.channel_joint_probability(dat)
eegfun.plot_joint_probability(jp)

