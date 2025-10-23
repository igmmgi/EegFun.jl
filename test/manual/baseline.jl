using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

eegfun.baseline!(dat)

eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 1));
dat.data

eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 2));
dat.data

eegfun.baseline!(dat, eegfun.IntervalIndex(start = 1, stop = 10));
dat.data

eegfun.baseline!(dat, eegfun.IntervalTime(start = 0.0, stop = 0.0));
dat.data

eegfun.baseline!(dat, eegfun.IntervalTime(start = 0.0, stop = 0.004));
dat.data
