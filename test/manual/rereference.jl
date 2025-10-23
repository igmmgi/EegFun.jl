using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

dat.data

eegfun.rereference!(dat, :Fp1)
dat.data

eegfun.rereference!(dat, :AF7)
dat.data

eegfun.rereference!(dat, :Fp1)
dat.data
