using EegFun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = EegFun.read_bdf(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);

dat.data

EegFun.rereference!(dat, :Fp1)
dat.data

EegFun.rereference!(dat, :AF7)
dat.data

EegFun.rereference!(dat, :Fp1)
dat.data
