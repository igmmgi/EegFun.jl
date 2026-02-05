using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/bdf/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.highpass_filter!(dat, 1)

# Rereference to :Fp1
EegFun.rereference!(dat, :Fp1)
dat.data

# Rereference to :AF7
EegFun.rereference!(dat, :AF7)
dat.data

# We can also see the influence of reference using the databrowser
EegFun.plot_databrowser(dat)


################################
# Epoching
################################

# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))  # -200 to 1000 ms

EegFun.rereference!(epochs, :Fp1)
EegFun.all_data(epochs[1]) # first condition
EegFun.all_data(epochs[2]) # second condition


################################
# ERPs
################################
# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.rereference!(erps, :Fp1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])

EegFun.rereference!(erps, :F1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])