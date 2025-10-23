using eegfun
# using BenchmarkTools

dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);
dat.sample_rate
eegfun.trigger_count(dat);

dat_new = eegfun.resample(dat, 2)
dat_new.sample_rate
eegfun.trigger_count(dat_new);  # triggers should be preserved

dat_new = eegfun.resample(dat, 4)
dat_new.sample_rate
eegfun.trigger_count(dat_new);  # triggers should be preserved

# mutating version
eegfun.resample!(dat, 2)
dat.sample_rate
eegfun.trigger_count(dat);  # triggers should be preserved


#################################
# Epoched DataFrameEeg
#################################
dat = eegfun.read_bdf("../Flank_C_3.bdf");
layout = eegfun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = eegfun.create_eeg_dataframe(dat, layout);

epoch_cfg = [eegfun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = eegfun.extract_epochs(dat, epoch_cfg, -2, 4)

epochs_new = eegfun.resample(epochs, 2)
epochs_new[1].sample_rate

eegfun.resample!(epochs, 2)
epochs[1].sample_rate
