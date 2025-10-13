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


