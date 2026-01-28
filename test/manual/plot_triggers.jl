using EegFun
using GLMakie
using BenchmarkTools

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
# data_file = joinpath(@__DIR__, "..", "..", "..", "Flank_C_3.bdf")
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded", "Flank_C_12.bdf")
# data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_12.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file);
count = EegFun.trigger_count(dat)

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "exp2_1002")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.read_raw_data(data_file);
count = EegFun.trigger_count(dat)

dat = EegFun.create_eeg_dataframe(dat, layout_file);
count = EegFun.trigger_count(dat)

count = EegFun.trigger_count(dat)

EegFun.plot_trigger_overview(dat)
EegFun.plot_trigger_overview(dat; ignore_triggers = [3])

EegFun.plot_trigger_timing(dat)
EegFun.plot_trigger_timing(dat; ignore_triggers = [3])
