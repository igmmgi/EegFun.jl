
# Mirror Padding Demo {#Mirror-Padding-Demo}

Demonstrates mirror padding to reduce edge effects in filtering.

## Overview {#Overview}

Mirror padding replicates data at epoch boundaries in reverse order to mitigate filter edge artifacts. Supports `:pre`, `:post`, or `:both` padding modes.

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
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.1)

# Create some epoched data
epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.2,1.0)  # -200 to 1000 ms

EegFun.plot_epochs(epochs, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :pre)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :post)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

epochs_new = EegFun.mirror(epochs, :both)
EegFun.plot_epochs(epochs_new, channel_selection = EegFun.channels([:Fp1]))

# ERPs
erps = EegFun.average_epochs(epochs)

EegFun.plot_erp(erps, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :pre)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :post)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))

erps_new = EegFun.mirror(erps, :both)
EegFun.plot_erp(erps_new, channel_selection = EegFun.channels([:Fp1]))
```


## See Also {#See-Also}
- [Preprocessing Reference](../reference/preprocessing.md)
  
