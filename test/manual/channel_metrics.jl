using eegfun
using GLMakie

# Load some test data
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv")
eegfun.polar_to_cartesian_xy!(layout_file)
eegfun.polar_to_cartesian_xyz!(layout_file)

dat = eegfun.read_bdf(data_file)
dat = eegfun.create_eeg_dataframe(dat, layout_file)
eegfun.rereference!(dat, :avg)
eegfun.filter_data!(dat, "hp", 1)

channel_joint_probability = eegfun.channel_joint_probability(dat)
