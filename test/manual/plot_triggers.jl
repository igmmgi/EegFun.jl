using eegfun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)
dat = eegfun.read_bdf(data_file);

count = eegfun.trigger_count(dat)

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "exp2_1002")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)

dat = eegfun.read_bdf(data_file);
count = eegfun.trigger_count(dat)

dat = eegfun.create_eeg_dataframe(dat, layout_file);
count = eegfun.trigger_count(dat)

count = eegfun.trigger_count(dat)

eegfun.plot_trigger_overview(dat)
eegfun.plot_trigger_overview(dat; ignore_triggers = [3])

eegfun.plot_trigger_timing(dat)
eegfun.plot_trigger_timing(dat; ignore_triggers = [3])
