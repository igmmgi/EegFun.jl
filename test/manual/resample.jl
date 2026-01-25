using EegFun
# using BenchmarkTools

dat = EegFun.read_bdf("../Flank_C_3.bdf");
layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = EegFun.create_eeg_dataframe(dat, layout);
dat.sample_rate
EegFun.trigger_count(dat);

dat_new = EegFun.resample(dat, 2)
dat_new.sample_rate
EegFun.trigger_count(dat_new);  # triggers should be preserved

dat_new = EegFun.resample(dat, 4)
dat_new.sample_rate
EegFun.trigger_count(dat_new);  # triggers should be preserved

# mutating version
EegFun.resample!(dat, 2)
dat.sample_rate
EegFun.trigger_count(dat);  # triggers should be preserved


#################################
# Epoched DataFrameEeg
#################################
dat = EegFun.read_bdf("../Flank_C_3.bdf");
layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = EegFun.create_eeg_dataframe(dat, layout);

epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

epochs_new = EegFun.resample(epochs, 2)
epochs_new[1].sample_rate

EegFun.resample!(epochs, 2)
epochs[1].sample_rate
