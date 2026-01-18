using eegfun
# using BenchmarkTools

dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);

epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

epochs_new = eegfun.mirror(epochs[1], :both)
epochs[1].data
epochs_new.data

epochs_new = eegfun.mirror(epochs, :both)
epochs[1].data
epochs_new[1].data
