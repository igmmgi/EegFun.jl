# EegFun.jl

![EegFun Logo](images/EegFunLogo.png)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igmmgi.github.io/EegFun.jl/stable)
[![Build Status](https://github.com/igmmgi/EegFun.jl/workflows/Documentation/badge.svg)](https://github.com/igmmgi/EegFun.jl/actions)
[![CI](https://github.com/igmmgi/EegFun.jl/workflows/Tests/badge.svg)](https://github.com/igmmgi/EegFun.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia package for EEG/ERP data analysis and visualization. Currently under active development (Alpha 0.1).

## Quick Start

```julia
using EegFun

data_file = "my_raw_file.bdf"
layout_file = EegFun.read_layout("my_layout.csv");
EegFun.polar_to_cartesian_xy!(layout_file)
dat = EegFun.read_raw_data(data_file);
dat = EegFun.create_eeg_dataframe(dat, layout_file);

EegFun.plot_databrowser(dat);
```

## Features & Visualizations

### Interactive Data Browser

Inspect raw EEG data, mark artifacts, and apply filters interactively.

![Data Browser](images/data_browser.png)

### Automated Artifact Detection

Detect various types of artifacts using customizable criteria.

![Artifact Detection](images/artifact_detection.png)

### Topographical Plots

Visualize ERP distribution across the scalp.

![ERP Topography](images/erp_topo_layout.png)
