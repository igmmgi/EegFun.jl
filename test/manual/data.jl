using EegFun

# read raw data
dat = EegFun.read_raw_data("./data/files/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# data/channel access functions
EegFun.all_data(dat)
EegFun.meta_data(dat)
EegFun.channel_data(dat)
EegFun.extra_data(dat)

# some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

EegFun.all_data(epochs)
EegFun.all_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.meta_data(epochs)
EegFun.meta_data(epochs, epoch_selection = EegFun.epochs(1:2))
EegFun.channel_data(epochs)
EegFun.channel_data(epochs, epoch_selection = EegFun.epochs(1:2))

# test subsetting
dat_subset = EegFun.subset(dat, channel_selection = EegFun.channels([:Fp1, :Fp2]))
dat_subset = EegFun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = EegFun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds
