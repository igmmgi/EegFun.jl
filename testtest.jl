# TODO: what should we put in here?
# TODO: how much does it readlly help?
using EegFun


# Reading *.bdf files and opening databrowser
# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf")
# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)
# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eegfun_data(dat, layout_file)
# Some minimal preprocessing (average reference and highpass filter)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)
EegFun.lowpass_filter!(dat, 50)
# Basic databrowser
EegFun.plot_databrowser(dat)