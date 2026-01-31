
# Rereferencing Demo {#Rereferencing-Demo}

Demonstrates different rereferencing schemes for EEG data.

## Overview {#Overview}

Rereferencing changes the voltage reference for all channels. Common schemes include average reference, single channel reference, and mastoid reference.

## Source Code {#Source-Code}

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference and highpass filter)
EegFun.highpass_filter!(dat, 1)


EegFun.rereference!(dat, :Fp1)
dat.data

EegFun.rereference!(dat, :AF7)
dat.data

EegFun.rereference!(dat, :Fp1)
dat.data

EegFun.plot_databrowser(dat)


# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.2, 1.0)  # -200 to 1000 ms

EegFun.rereference!(epochs, :Fp1)
EegFun.all_data(epochs[1])
EegFun.all_data(epochs[2])

EegFun.rereference!(epochs, :F1)
EegFun.all_data(epochs[1])
EegFun.all_data(epochs[2])


# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.rereference!(erps, :Fp1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])

EegFun.rereference!(erps, :F1)
EegFun.all_data(erps[1])
EegFun.all_data(erps[2])
```


## See Also {#See-Also}
- [Preprocessing Reference](../reference/preprocessing.md)
  
