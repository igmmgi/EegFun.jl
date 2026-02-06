
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igmmgi.github.io/EegFun.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Julia 1.10+](https://img.shields.io/badge/julia-v1.10+-blue.svg)
![Status: Alpha](https://img.shields.io/badge/Status-Alpha-orange.svg)


# EegFun.jl

<img src="images/EegFunLogo.png" alt="EegFun Logo" width="150"/>

A Julia package for EEG/ERP data analysis and visualization. Currently under active development.

## Documentation

**[View the full documentation →](https://igmmgi.github.io/EegFun.jl/dev/)**

The documentation includes:

* Installation instructions and getting started guide
* Tutorials and Demos
* Complete API reference

## Installation

EegFun.jl is currently unregistered. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/igmmgi/EegFun.jl")
```

## Features

* **EEG/ERP Analysis**
* **EEG/ERP Interactive Plots (via Makie.jl)**
* **Time-Frequency Analysis**
* **Raw data to full ERP batch preprocessing pipelines**

## Showcase

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

<details>
<summary><b>Plot Examples</b></summary>

<details>
<summary>Artifact Detection</summary>

<img src="images/artifact_detection.png" alt="Artifact Detection" width="600"/>

</details>

<details>
<summary>Epoch Plots (Grid Layout)</summary>

<img src="images/epochs_grid_layout.png" alt="Epochs Grid Layout" width="800"/>

</details>

<details>
<summary>ERP (Topo Layout)</summary>

<img src="images/erp_topo_layout.png" alt="ERP Topo Layout" width="600"/>

</details>

<details>
<summary>ERP Image (Topo Layout)</summary>

<img src="images/erp_image_topo_layout.png" alt="ERP Image Topo Layout" width="600"/>

</details>

</details>

<details>
<summary><i>TODO</i></summary>

* Add additional file formats to read_raw_data (currently only Biosemi BDF and BrainVision)
* Lots more ....

</details>
