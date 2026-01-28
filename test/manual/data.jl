using EegFun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "Flank_C_3.bdf")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.read_raw_data(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# data/channel access functions
EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.channel_data(dat)
EegFun.extra_data(dat)

EegFun.all_labels(dat)
EegFun.channel_labels(dat)
EegFun.meta_labels(dat)
EegFun.extra_labels(dat)

# some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[3]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

EegFun.all_data(epochs)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.meta_data(epochs)
EegFun.meta_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.channel_data(epochs)
EegFun.channel_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.extra_data(epochs)
EegFun.extra_data(epochs, epoch_selection = EegFun.epochs(1:2))

# test subsetting
dat_subset = EegFun.subset(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
dat_subset = EegFun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = EegFun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds

# test subsetting
epoch_subset = EegFun.subset(epochs, channel_selection = EegFun.channels([:Fp1, :Fp2]))
epoch_subset = EegFun.subset(epochs, epoch_selection = EegFun.epochs(1:2))
