using eegfun
using GLMakie

# Get some basic data with initial preprocessing steps (high-pass filter, epoch)
data_file = joinpath(@__DIR__, "..", "..", "..",  "Flank_C_3.bdf")
layout_file = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
eegfun.polar_to_cartesian_xy!(layout_file)

dat = eegfun.read_bdf(data_file);
dat = eegfun.create_eeg_dataframe(dat, layout_file);

# data/channel access functions
eegfun.all_data(dat)
eegfun.meta_data(dat)
eegfun.channel_data(dat)
eegfun.extra_data(dat)

eegfun.all_labels(dat)
eegfun.channel_labels(dat)
eegfun.meta_labels(dat)
eegfun.extra_labels(dat)

# some epoched data
epoch_cfg = [
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[3]]),
]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

eegfun.all_data(epochs)
eegfun.all_data(epochs, epoch_selection = eegfun.epochs(1:2))
eegfun.meta_data(epochs)
eegfun.meta_data(epochs, epoch_selection = eegfun.epochs(1:2))
eegfun.channel_data(epochs)
eegfun.channel_data(epochs, epoch_selection = eegfun.epochs(1:2))
eegfun.extra_data(epochs)
eegfun.extra_data(epochs, epoch_selection = eegfun.epochs(1:2))

# test subsetting
dat_subset = eegfun.subset(dat, channel_selection = eegfun.channels([:Fp1, :Fp2]))
dat_subset = eegfun.subset(dat, sample_selection = x -> x.sample .<= 10_000) # first 10000 samples
dat_subset = eegfun.subset(dat, sample_selection = x -> x.time .<= 10) # first 10 seconds

# test subsetting
epoch_subset = eegfun.subset(epochs, channel_selection = eegfun.channels([:Fp1, :Fp2]))
epoch_subset = eegfun.subset(epochs, epoch_selection = eegfun.epochs(1:2))

