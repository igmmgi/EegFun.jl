
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igmmgi.github.io/EegFun.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# EegFun.jl

<img src="images/EegFunLogo.png" alt="EegFun Logo" width="150"/>

A Julia package for EEG/ERP data analysis and visualization. Currently under active development (Alpha 0.1).

## Documentation

**[View the full documentation →](https://igmmgi.github.io/EegFun.jl/dev/)**

The documentation includes:

* Getting started guide and installation instructions
* Interactive demos and tutorials
* Complete API reference

## Features

* **EEG/ERP Analysis**
* **EEG/ERP Interactive Plots (via Makie.jl)**
* **Time-Frequency Analysis**
* **Raw data to full ERP batch preprocessing pipelines**

## Some Features

### Data Browser

<details>
<summary>Show Code</summary>

```julia
using EegFun

# raw data file and channel coordinates
dat = EegFun.read_raw_data("my_raw_file.bdf");
layout_file = EegFun.read_layout("my_layout.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# Julia EegFun type
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.plot_databrowser(dat);
```

</details>

<img src="images/data_browser.png" alt="Data Browser" width="800"/>

### ICA Data Browser

<details>
<summary>Show Code</summary>

```julia
using EegFun

# raw data file and channel coordinates
dat = EegFun.read_raw_data("my_raw_file.bdf");

layout_file = EegFun.read_layout("my_layout.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

dat = EegFun.create_eeg_dataframe(dat, layout_file);

# rereference data and apply 1Hz high-pass filter for ICA
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)

# calculate EOG channels
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # horizontal EOG = F9 - F10

# detect some extreme values
EegFun.is_extreme_value!(dat, 200);

# ICA on continuous data without extreme values
ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_200)) 

EegFun.plot_ica_component_activation(dat, ica_result)
```

</details>

<img src="images/data_browser_ica.png" alt="Data Browser ICA" width="800"/>

### Plot Examples

#### Artifact Detection

<img src="images/artifact_detection.png" alt="Artifact Detection" width="600"/>

<details>
<summary>More Plot Examples</summary>

#### Epoch Plots (Grid Layout)

<img src="images/epochs_grid_layout.png" alt="Epochs Grid Layout" width="800"/>

#### ERP (Topo Layout)

<img src="images/erp_topo_layout.png" alt="ERP Topo Layout" width="600"/>

#### ERP Image (Topo Layout)

<img src="images/erp_image_topo_layout.png" alt="ERP Image Topo Layout" width="600"/>

</details>

## TODO

* Add additional file formats to read_raw_data (currently only Biosemi BDF and BrainVision) []
* Tutorial examples []
* Lots more ....
