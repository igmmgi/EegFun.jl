using EegFun
# using BenchmarkTools

dat = EegFun.read_raw_data("../Flank_C_3.bdf");
layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv");
dat = EegFun.create_eeg_dataframe(dat, layout);

epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -2, 4)

epochs_new = EegFun.mirror(epochs[1], :both)
epochs[1].data
epochs_new.data

epochs_new = EegFun.mirror(epochs, :both)
epochs[1].data
epochs_new[1].data
