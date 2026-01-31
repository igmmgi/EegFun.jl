
# Available Electrode Layouts {#Available-Electrode-Layouts}

EegFun.jl provides electrode layout files for various EEG cap manufacturers and configurations.

## Usage {#Usage}

```julia
using EegFun

# Load layout
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")

# Convert polar to Cartesian coordinates for plotting
EegFun.polar_to_cartesian_xy!(layout)

# Use with your data
dat = EegFun.create_eeg_dataframe(raw_data, layout)
```


## Biosemi ActiveTwo {#Biosemi-ActiveTwo}

### 16 Channels {#16-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi16.png" alt="Biosemi 16">


### 32 Channels {#32-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi32.png" alt="Biosemi 32">


### 64 Channels {#64-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi64.png" alt="Biosemi 64">


### 72 Channels (64 + 8 external) {#72-Channels-64-8-external}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi72.png" alt="Biosemi 72">


### 128 Channels {#128-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi128.png" alt="Biosemi 128">


### 160 Channels {#160-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi160.png" alt="Biosemi 160">


### 256 Channels {#256-Channels}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi256.png" alt="Biosemi 256">


## Easycap {#Easycap}

### M1 {#M1}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM1.png" alt="Easycap M1">


### M3 {#M3}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM3.png" alt="Easycap M3">


### M7 {#M7}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM7.png" alt="Easycap M7">


### M10 {#M10}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM10.png" alt="Easycap M10">


### M11 {#M11}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM11.png" alt="Easycap M11">


### M14 {#M14}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM14.png" alt="Easycap M14">


### M15 {#M15}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM15.png" alt="Easycap M15">


### M16 {#M16}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM16.png" alt="Easycap M16">


### M17 {#M17}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM17.png" alt="Easycap M17">


### M20 {#M20}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM20.png" alt="Easycap M20">


### M22 {#M22}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM22.png" alt="Easycap M22">


### M23 {#M23}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM23.png" alt="Easycap M23">


### M24 {#M24}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM24.png" alt="Easycap M24">


### M25 {#M25}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM25.png" alt="Easycap M25">


### M64 {#M64}

<img src="https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM64.png" alt="Easycap M64">


## Brain Products actiCap {#Brain-Products-actiCap}

Brain Products offers extensive layout options (60+ configurations). For the complete collection, browse the [acticap folder](https://github.com/igmmgi/EegFun.jl/tree/main/resources/layouts/acticap).

**Available systems:**
- **Standard BrainCap**: 22, 32, 64, 96, 128, 256 channels
  
- **BrainCap MR**: MRI-compatible (32, 64, 96, 128 channels)
  
- **BrainCap MR3**: MR3 variants
  
- **BrainCap Sleep**: Sleep monitoring (15 AASM, 32, 64 channels)
  
- **R-Net Systems**: Wireless caps for actiCHamp, BrainAmp, LiveAmp
  
- **LiveCap**: Portable (32, 64 channels)
  
- **Specialized**: TMS-compatible and MEG-compatible
  

## See Also {#See-Also}
- [Layouts API Reference](layouts.md)
  
- [Plot Layout Demo](../demos/plot-layout.md)
  
