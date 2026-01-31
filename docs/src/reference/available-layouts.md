# Available Electrode Layouts

EegFun.jl provides electrode layout files for various EEG cap manufacturers and configurations.

## Usage

```julia
using EegFun

# Load layout
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")

# Convert polar to Cartesian coordinates for plotting
EegFun.polar_to_cartesian_xy!(layout)

# Use with your data
dat = EegFun.create_eeg_dataframe(raw_data, layout)
```

## Biosemi ActiveTwo

### 16 Channels

![Biosemi 16](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi16.png)

### 32 Channels

![Biosemi 32](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi32.png)

### 64 Channels

![Biosemi 64](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi64.png)

### 72 Channels (64 + 8 external)

![Biosemi 72](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi72.png)

### 128 Channels

![Biosemi 128](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi128.png)

### 160 Channels

![Biosemi 160](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi160.png)

### 256 Channels

![Biosemi 256](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/biosemi/biosemi256.png)

## Easycap

### M1

![Easycap M1](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM1.png)

### M3

![Easycap M3](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM3.png)

### M7

![Easycap M7](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM7.png)

### M10

![Easycap M10](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM10.png)

### M11

![Easycap M11](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM11.png)

### M14

![Easycap M14](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM14.png)

### M15

![Easycap M15](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM15.png)

### M16

![Easycap M16](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM16.png)

### M17

![Easycap M17](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM17.png)

### M20

![Easycap M20](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM20.png)

### M22

![Easycap M22](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM22.png)

### M23

![Easycap M23](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM23.png)

### M24

![Easycap M24](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM24.png)

### M25

![Easycap M25](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM25.png)

### M64

![Easycap M64](https://raw.githubusercontent.com/igmmgi/EegFun.jl/main/resources/layouts/easycap/easycapM64.png)

## Brain Products actiCap

Brain Products offers extensive layout options (60+ configurations). For the complete collection, browse the [acticap folder](https://github.com/igmmgi/EegFun.jl/tree/main/resources/layouts/acticap).

**Available systems:**

- **Standard BrainCap**: 22, 32, 64, 96, 128, 256 channels
- **BrainCap MR**: MRI-compatible (32, 64, 96, 128 channels)
- **BrainCap MR3**: MR3 variants
- **BrainCap Sleep**: Sleep monitoring (15 AASM, 32, 64 channels)
- **R-Net Systems**: Wireless caps for actiCHamp, BrainAmp, LiveAmp
- **LiveCap**: Portable (32, 64 channels)
- **Specialized**: TMS-compatible and MEG-compatible

## See Also

- [Layouts API Reference](layouts.md)
- [Plot Layout Demo](../demos/plot-layout.md)
