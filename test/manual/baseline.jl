using EegFun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");

dat = EegFun.read_raw_data(data_file)
dat = EegFun.create_eeg_dataframe(dat, layout_file)

EegFun.baseline!(dat)

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 1));
dat.data

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 2));
dat.data

EegFun.baseline!(dat, EegFun.IntervalIndex(start = 1, stop = 10));
dat.data

EegFun.baseline!(dat, EegFun.IntervalTime(start = 0.0, stop = 0.0));
dat.data

EegFun.baseline!(dat, EegFun.IntervalTime(start = 0.0, stop = 0.004));
dat.data
