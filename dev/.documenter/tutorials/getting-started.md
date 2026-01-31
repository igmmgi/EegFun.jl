
# Getting Started with EegFun.jl {#Getting-Started-with-EegFun.jl}

Welcome to EegFun.jl! This tutorial will guide you through your first steps with EEG data analysis in Julia.

## Requirements {#Requirements}
- Julia 1.9 or later
  
- An internet connection for package installation
  
- Optional: Sample EEG data in BDF or similar format
  

## Installation {#Installation}

First, install EegFun.jl and a plotting backend:

```julia
using Pkg
Pkg.add("EegFun")
Pkg.add("GLMakie")  # For interactive plotting
```


Load the packages:

```julia
using EegFun
using GLMakie
```


## Loading Your First Dataset {#Loading-Your-First-Dataset}

EegFun.jl supports various EEG formats. Here&#39;s how to load Biosemi BDF data:

```julia
# Load raw EEG data
dat = EegFun.read_raw_data("path/to/your/data.bdf")

# Load electrode layout
layout = EegFun.read_layout("biosemi64.csv")

# Create an EEG dataframe with spatial information
dat = EegFun.create_eeg_dataframe(dat, layout)
```

> 
> **TODO**: Add information about where to find sample datasets
> 


## Basic Data Exploration {#Basic-Data-Exploration}

Inspect your data structure:

```julia
# View column names
names(dat.data)

# Check sampling rate
dat.sampling_rate

# Number of channels
length(dat.layout.labels)

# Number of samples
nrow(dat.data)
```


## Simple Visualization {#Simple-Visualization}

Plot a few seconds of raw data:

```julia
# TODO: Add example of plotting raw data
# fig, ax = EegFun.plot_raw(dat, duration=5.0)
```


## Next Steps {#Next-Steps}

Now that you have data loaded, you can:
- [Preprocess your data](basic-preprocessing.md) with filtering and artifact detection
  
- [Extract epochs](erp-analysis.md) around experimental events
  
- [Apply ICA](ica-workflow.md) for artifact removal
  
> 
> **TODO**: Complete this tutorial with more detailed examples
> - Add screenshots/plots
>   
> - Include troubleshooting section
>   
> - Add links to sample datasets
>   
> 

