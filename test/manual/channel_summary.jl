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
eegfun.filter_data!(dat, "hp", 0.1)

summary = eegfun.channel_summary(dat)
eegfun.log_pretty_table(summary; title = "Initial Channel Summary")

eegfun.is_extreme_value!(dat, 100);
summary = eegfun.channel_summary(dat, sample_selection = eegfun.samples_not(:is_extreme_value_100))
