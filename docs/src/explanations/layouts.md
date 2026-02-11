# Electrode Layouts

EegFun.jl includes a comprehensive collection of electrode layout files that define the spatial positions of EEG electrodes for various montage systems. These layouts are essential for topographic plotting, channel interpolation, and spatial analyses.

## Overview

Layout files are stored in CSV format with three required columns:

- `label`: Channel name (Symbol)
- `inc`: Incidence angle in degrees (vertical angle from vertex)
- `azi`: Azimuth angle in degrees (horizontal angle)

The polar coordinates (inc, azi) serve as the **Source of Truth** for electrode positions. Cartesian coordinates (2D and 3D) are computed on-demand from these polar coordinates.

## Sources

The electrode layout files included in EegFun.jl were obtained from the following sources:

- **BioSemi**: Downloaded from [BioSemi Head Cap Specifications](https://www.biosemi.com/headcap.htm)
- **BrainProducts actiCap**: Downloaded from [BrainProducts Cap Montages](https://www.brainproducts.com/downloads/cap-montages/)
- **EasyCap**: Generated with reference to:
  - [FieldTrip](https://www.fieldtriptoolbox.org/) template layouts (`fieldtrip/template/layout/`)
  - [MNE-Python](https://mne.tools/) channel montages (`mne-python/mne/channels/data/montages`)
  - [EasyCap electrode layout documentation (PDF)](https://www.easycap.de) - Converted to polar coordinates (label, inc, azi) for consistency with BioSemi format

## Usage

```julia
# Load a layout file
layout = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv")

# Convert to 2D Cartesian coordinates (preserving radial distances)
EegFun.polar_to_cartesian_xy!(layout)

# Calculate neighbors 
EegFun.get_neighbours_xy!(layout, 0.3)

# Plot the layout with neighbours; hover over electrodes to see neighbours
EegFun.plot_layout_2d(layout; neighbours=true)

EegFun.get_neighbours_xyz!(layout, 0.5)
EegFun.plot_layout_3d(layout; neighbours=true)
```

## Available Layouts

### BioSemi

Standard BioSemi ActiveTwo electrode caps and extended montages:

::: details biosemi16 (16 channels) - Standard 16-channel montage
![biosemi16 layout](/layouts/biosemi/biosemi16.png)
:::

::: details biosemi32 (32 channels) - Standard 32-channel montage
![biosemi32 layout](/layouts/biosemi/biosemi32.png)
:::

::: details biosemi64 (64 channels) - Standard 64-channel montage
![biosemi64 layout](/layouts/biosemi/biosemi64.png)
:::

::: details biosemi70 (70 channels) - 64 + 6 external channels
![biosemi70 layout](/layouts/biosemi/biosemi70.png)
:::

::: details biosemi72 (72 channels) - 64 + 8 external channels
![biosemi72 layout](/layouts/biosemi/biosemi72.png)
:::

::: details biosemi128 (128 channels) - Standard 128-channel montage
![biosemi128 layout](/layouts/biosemi/biosemi128.png)
:::

::: details biosemi160 (160 channels) - Extended high-density montage
![biosemi160 layout](/layouts/biosemi/biosemi160.png)
:::

::: details biosemi256 (256 channels) - Ultra high-density montage
![biosemi256 layout](/layouts/biosemi/biosemi256.png)
:::

### BrainProducts actiCap

BrainProducts electrode caps organized by product line:

#### Standard BrainCap

::: details BC-22 (22 channels) - 22-channel clinical montage
![BC-22 layout](/layouts/acticap/BrainCap/BrainCap_22_Channel/BC-22.png)
:::

::: details BC-32 (32 channels) - Standard 32-channel montage
![BC-32 layout](/layouts/acticap/BrainCap/BrainCap_32_Channel/BC-32.png)
:::

::: details BC-64 (64 channels) - Standard 64-channel montage
![BC-64 layout](/layouts/acticap/BrainCap/BrainCap_64_Channel/BC-64.png)
:::

::: details BC-96 (96 channels) - High-density 96-channel montage
![BC-96 layout](/layouts/acticap/BrainCap/BrainCap_96_Channel/BC-96.png)
:::

::: details BC-128 (128 channels) - High-density 128-channel montage
![BC-128 layout](/layouts/acticap/BrainCap/BrainCap_128_Channel/BC-128.png)
:::

::: details BC-256 (256 channels) - Ultra high-density montage
![BC-256 layout](/layouts/acticap/BrainCap/BrainCap_256_Channel/BC-256.png)
:::

#### BrainCap MR (MRI-compatible)

::: details BC-MR-32 (32 channels)
![BC-MR-32 layout](/layouts/acticap/BrainCap_MR/BrainCap_MR_32_Channel/BC-MR-32.png)
:::

::: details BC-MR-64 (64 channels)
![BC-MR-64 layout](/layouts/acticap/BrainCap_MR/BrainCap_MR_64_Channel/BC-MR-64.png)
:::

::: details BC-MR-96 (96 channels)
![BC-MR-96 layout](/layouts/acticap/BrainCap_MR/BrainCap_MR_96_Channel/BC-MR-96.png)
:::

::: details BC-MR-128 (128 channels)
![BC-MR-128 layout](/layouts/acticap/BrainCap_MR/BrainCap_MR_128_Channel/BC-MR-128.png)
:::

#### BrainCap MR3 (3T MRI-compatible)

::: details BC-MR3-32 (32 channels)
![BC-MR3-32 layout](/layouts/acticap/BrainCap_MR3/BrainCap_MR3_32_Channel/BC-MR3-32.png)
:::

::: details BC-MR3-64 (64 channels)
![BC-MR3-64 layout](/layouts/acticap/BrainCap_MR3/BrainCap_MR3_64_Channel/BC-MR3-64.png)
:::

::: details BC-MR3-96 (96 channels)
![BC-MR3-96 layout](/layouts/acticap/BrainCap_MR3/BrainCap_MR3_96_Channel/BC-MR3-96.png)
:::

::: details BC-MR3-128 (128 channels)
![BC-MR3-128 layout](/layouts/acticap/BrainCap_MR3/BrainCap_MR3_128_Channel/BC-MR3-128.png)
:::

#### BrainCap MEG (MEG-compatible)

::: details BC-MEG-32 (32 channels)
![BC-MEG-32 layout](/layouts/acticap/BrainCap_MEG/BrainCap_MEG_32/BC-MEG-32.png)
:::

::: details BC-MEG-64 (64 channels)
![BC-MEG-64 layout](/layouts/acticap/BrainCap_MEG/BrainCap_MEG_64/BC-MEG-64.png)
:::

::: details BC-MEG-128 (128 channels)
![BC-MEG-128 layout](/layouts/acticap/BrainCap_MEG/BrainCap_MEG_128/BC-MEG-128.png)
:::

#### BrainCap TMS (TMS-compatible)

::: details BC-TMS-32 (32 channels)
![BC-TMS-32 layout](/layouts/acticap/BrainCap_TMS/BrainCap_TMS_32/BC-TMS-32.png)
:::

::: details BC-TMS-64 (64 channels)
![BC-TMS-64 layout](/layouts/acticap/BrainCap_TMS/BrainCap_TMS_64/BC-TMS-64.png)
:::

::: details BC-TMS-96 (96 channels)
![BC-TMS-96 layout](/layouts/acticap/BrainCap_TMS/BrainCap_TMS_96/BC-TMS-96.png)
:::

::: details BC-TMS-128 (128 channels)
![BC-TMS-128 layout](/layouts/acticap/BrainCap_TMS/BrainCap_TMS_128/BC-TMS-128.png)
:::

#### BrainCap Sleep

::: details BC-MRS-15 (15 channels) - MR sleep - AASM standard
![BC-MRS-15 layout](/layouts/acticap/BrainCap_MR_Sleep/BrainCap_MR_Sleep_15_AASM/BC-MRS-15.png)
:::

::: details BC-MRS-32 (32 channels) - MR sleep montage
![BC-MRS-32 layout](/layouts/acticap/BrainCap_MR_Sleep/BrainCap_MR_Sleep_32/BC-MRS-32.png)
:::

::: details BC-MRS-64 (64 channels) - MR sleep montage
![BC-MRS-64 layout](/layouts/acticap/BrainCap_MR_Sleep/BrainCap_MR_sleep_64/BC-MRS-64.png)
:::

::: details BC-MRS3-15 (15 channels) - MR3 sleep - AASM standard
![BC-MRS3-15 layout](/layouts/acticap/BrainCap_MR3_Sleep/BrainCap_MR3_Sleep_15_AASM/BC-MRS3-15.png)
:::

::: details BC-MRS3-32 (32 channels) - MR3 sleep montage
![BC-MRS3-32 layout](/layouts/acticap/BrainCap_MR3_Sleep/BrainCap_MR3_Sleep_32/BC-MRS3-32.png)
:::

::: details BC-MRS3-64 (64 channels) - MR3 sleep montage
![BC-MRS3-64 layout](/layouts/acticap/BrainCap_MR3_Sleep/BrainCap_MR3_sleep_64/BC-MRS3-64.png)
:::

::: details BC-SL-32 (32 channels) - BrainAmp sleep
![BC-SL-32 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_BrainAmp/BC-SL-32.png)
:::

::: details BC-SL-64 (64 channels) - BrainAmp sleep
![BC-SL-64 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_BrainAmp/BC-SL-64.png)
:::

::: details BC-SL-128 (128 channels) - BrainAmp sleep
![BC-SL-128 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_BrainAmp/BC-SL-128.png)
:::

::: details BCA-SL-32 (32 channels) - actiCHamp sleep
![BCA-SL-32 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_actiCHamp/BCA-SL-32.png)
:::

::: details BCA-SL-64 (64 channels) - actiCHamp sleep
![BCA-SL-64 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_actiCHamp/BCA-SL-64.png)
:::

::: details BCA-SL-128 (128 channels) - actiCHamp sleep
![BCA-SL-128 layout](/layouts/acticap/BrainCap_Sleep_BCA-SL_BC-SL/BrainCap_Sleep_for_actiCHamp/BCA-SL-128.png)
:::

#### LiveCap

::: details LC-32 (32 channels) - LiveAmp 32-channel
![LC-32 layout](/layouts/acticap/LiveCap_LC/LiveCap_for_LiveAmp_32/LC-32.png)
:::

::: details LC-64 (64 channels) - LiveAmp 64-channel
![LC-64 layout](/layouts/acticap/LiveCap_LC/LiveCap_for_LiveAmp_64/LC-64.png)
:::

::: details LC-SL-32_Ch1-24 (32 channels) - LiveCap sleep (24 + 8)
![LC-SL-32_Ch1-24 layout](/layouts/acticap/LiveCap_LC-SL/LiveCAP_sleep_24_+_8/LC-SL-32_Ch1-24.png)
:::

#### R-Net

::: details RNP-BA-32 (32 channels) - R-Net for BrainAmp
![RNP-BA-32 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_BrainAmp_RNP-BA/RNP-BA-32.png)
:::

::: details RNP-BA-64 (64 channels) - R-Net for BrainAmp
![RNP-BA-64 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_BrainAmp_RNP-BA/RNP-BA-64.png)
:::

::: details RNP-BA-96 (96 channels) - R-Net for BrainAmp
![RNP-BA-96 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_BrainAmp_RNP-BA/RNP-BA-96.png)
:::

::: details RNP-BA-128 (128 channels) - R-Net for BrainAmp
![RNP-BA-128 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_BrainAmp_RNP-BA/RNP-BA-128.png)
:::

::: details RNP-LA-32 (32 channels) - R-Net for LiveAmp
![RNP-LA-32 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_LiveAmp_RNP-LA/RNP-LA-32.png)
:::

::: details RNP-LA-64 (64 channels) - R-Net for LiveAmp
![RNP-LA-64 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_LiveAmp_RNP-LA/RNP-LA-64.png)
:::

::: details RNP-AC-32 (32 channels) - R-Net for actiCHamp Plus
![RNP-AC-32 layout](/layouts/acticap/R-Net_RNP-AC_RNPA-AC_RNPA-AP/R-Net_for_actiCHamp_Plus_RNP-AC/RNP-AC-32/RNP-AC-32.png)
:::

::: details RNP-AC-64 (64 channels) - R-Net for actiCHamp Plus (with Iz)
![RNP-AC-64 layout](/layouts/acticap/R-Net_RNP-AC_RNPA-AC_RNPA-AP/R-Net_for_actiCHamp_Plus_RNP-AC/RNP-AC-64_with_Iz/RNP-AC-64.png)
:::

::: details RNP-AC-96 (96 channels) - R-Net for actiCHamp Plus (with Iz)
![RNP-AC-96 layout](/layouts/acticap/R-Net_RNP-AC_RNPA-AC_RNPA-AP/R-Net_for_actiCHamp_Plus_RNP-AC/RNP-AC-96_with_Iz/RNP-AC-96.png)
:::

::: details RNP-AC-128 (128 channels) - R-Net for actiCHamp Plus (with Iz)
![RNP-AC-128 layout](/layouts/acticap/R-Net_RNP-AC_RNPA-AC_RNPA-AP/R-Net_for_actiCHamp_Plus_RNP-AC/RNP-AC-128_with_Iz/RNP-AC-128.png)
:::

::: details RNP-BA-MR3-32 (32 channels) - R-Net MR3
![RNP-BA-MR3-32 layout](/layouts/acticap/R-Net_MR/RNP-BA-MR3-32.png)
:::

::: details RNP-BA-MR3-64 (64 channels) - R-Net MR3
![RNP-BA-MR3-64 layout](/layouts/acticap/R-Net_MR/RNP-BA-MR3-64.png)
:::

::: details RNP-BA-MR3-96 (96 channels) - R-Net MR3
![RNP-BA-MR3-96 layout](/layouts/acticap/R-Net_MR/RNP-BA-MR3-96.png)
:::

::: details RNP-BA-MR3-128 (128 channels) - R-Net MR3
![RNP-BA-MR3-128 layout](/layouts/acticap/R-Net_MR/RNP-BA-MR3-128.png)
:::

::: details RNP-AP-32 (32 channels) - R-Net for actiCHamp Plus
![RNP-AP-32 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_actiCHamp_Plus_RNP-AP/RNP-AP-32/RNP-AP-32.png)
:::

::: details RNP-AP-64 (64 channels) - R-Net for actiCHamp Plus (with FCz)
![RNP-AP-64 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_actiCHamp_Plus_RNP-AP/RNP-AP-64_with_FCz/RNP-AP-64.png)
:::

::: details RNP-AP-96 (96 channels) - R-Net for actiCHamp Plus (with FCz)
![RNP-AP-96 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_actiCHamp_Plus_RNP-AP/RNP-AP-96_with_FCz/RNP-AP-96.png)
:::

::: details RNP-AP-128 (128 channels) - R-Net for actiCHamp Plus (with FCz)
![RNP-AP-128 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_actiCHamp_Plus_RNP-AP/RNP-AP-128_with_FCz/RNP-AP-128.png)
:::

::: details RNP-AP-160 (160 channels) - R-Net for actiCHamp Plus (with FCz, Iz)
![RNP-AP-160 layout](/layouts/acticap/R-Net_RNP-AP_RNP-BA_RNP-LA/R-Net_for_actiCHamp_Plus_RNP-AP/RNP-AP-160_with_FCz_Iz/RNP-AP-160.png)
:::

#### actiCAP (Standard)

::: details AC-21 (21 channels)
![AC-21 layout](/layouts/acticap/actiCAP_CAC/actiCAP_21/AC-21.png)
:::

::: details AC-32 (32 channels)
![AC-32 layout](/layouts/acticap/actiCAP_CAC/actiCAP_32/AC-32.png)
:::

::: details AC-64 (64 channels)
![AC-64 layout](/layouts/acticap/actiCAP_CAC/actiCAP_64/AC-64.png)
:::

::: details AC-96 (96 channels)
![AC-96 layout](/layouts/acticap/actiCAP_CAC/actiCAP_96/AC-96.png)
:::

::: details AC-128 (128 channels)
![AC-128 layout](/layouts/acticap/actiCAP_CAC/actiCAP_128/AC-128.png)
:::

#### actiCAP Xpress Twist

::: details Twist-32 (32 channels) - for BrainAmp (with REF)
![Twist-32 REF layout](/layouts/acticap/actiCAP_Xpress_Twist_32/actiCAP_Xpress_Twist_for_BrainAmp/Twist-32_REF.png)
:::

::: details Twist-32 (32 channels) - for LiveAmp (with REF)
![Twist-32 REF layout](/layouts/acticap/actiCAP_Xpress_Twist_32/actiCAP_Xpress_Twist_for_LiveAmp/Twist-32_REF.png)
:::

::: details Twist-32 (32 channels) - for actiCHamp (no REF)
![Twist-32 NO REF layout](/layouts/acticap/actiCAP_Xpress_Twist_32/actiCAP_Xpress_Twist_for_actiCHamp/Twist-32_NO_REF.png)
:::

#### actiCHamp Electrode Cap (CAP)

::: details AP-32 (32 channels) - actiCHamp electrode cap
![AP-32 layout](/layouts/acticap/actiCAP_electrode_cap_CAP/actiCHamp_electrode_cap_32/AP-32.png)
:::

::: details AP-64 (64 channels) - actiCHamp electrode cap
![AP-64 layout](/layouts/acticap/actiCAP_electrode_cap_CAP/actiCHamp_electrode_cap_64/AP-64.png)
:::

::: details AP-96 (96 channels) - actiCHamp electrode cap
![AP-96 layout](/layouts/acticap/actiCAP_electrode_cap_CAP/actiCHamp_electrode_cap_96/AP-96.png)
:::

::: details AP-128 (128 channels) - actiCHamp electrode cap
![AP-128 layout](/layouts/acticap/actiCAP_electrode_cap_CAP/actiCHamp_electrode_cap_128/AP-128.png)
:::

::: details AP-160 (160 channels) - actiCHamp electrode cap
![AP-160 layout](/layouts/acticap/actiCAP_electrode_cap_CAP/actiCHamp_electrode_cap_160/AP-160.png)
:::

#### actiCAP UOL (CMA)

::: details CMA-21 (21 channels) - BrainAmp (with REF)
![CMA-21 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_21_Channel/BrainAmp/CMA-21_REF.png)
:::

::: details CMA-21 (21 channels) - actiCHamp (no REF)
![CMA-21 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_21_Channel/actiCHamp/CMA-21_NO_REF.png)
:::

::: details CMA-21 (21 channels) - actiCHamp Plus (no REF)
![CMA-21 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_21_Channel/actiCHamp_Plus/CMA-21_NO_REF.png)
:::

::: details CMA-32 (32 channels) - BrainAmp (with REF)
![CMA-32 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_32_Channel/BrainAmp/CMA-32_REF.png)
:::

::: details CMA-32 (32 channels) - actiCHamp (no REF)
![CMA-32 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_32_Channel/actiCHamp/CMA-32_NO_REF.png)
:::

::: details CMA-32 (32 channels) - actiCHamp Plus (no REF)
![CMA-32 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_32_Channel/actiCHamp_Plus/CMA-32_NO_REF.png)
:::

::: details CMA-64 (64 channels) - BrainAmp (with REF)
![CMA-64 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_64_Channel/BrainAmp/CMA-64_REF.png)
:::

::: details CMA-64 (64 channels) - actiCHamp (no REF)
![CMA-64 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_64_Channel/actiCHamp/CMA-64_NO_REF.png)
:::

::: details CMA-64 (64 channels) - actiCHamp Plus (no REF)
![CMA-64 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_64_Channel/actiCHamp_Plus/CMA-64_NO_REF.png)
:::

::: details CMA-96 (96 channels) - BrainAmp (with REF)
![CMA-96 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_96_Channel/BrainAmp/CMA-96_REF.png)
:::

::: details CMA-96 (96 channels) - actiCHamp (no REF)
![CMA-96 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_96_Channel/actiCHamp/CMA-96_NO_REF.png)
:::

::: details CMA-96 (96 channels) - actiCHamp Plus (no REF)
![CMA-96 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_96_Channel/actiCHamp_Plus/CMA-96_NO_REF.png)
:::

::: details CMA-128 (128 channels) - BrainAmp (with REF)
![CMA-128 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_128_Channel/BrainAmp/CMA-128_REF.png)
:::

::: details CMA-128 (128 channels) - actiCHamp (no REF)
![CMA-128 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_128_Channel/actiCHamp/CMA-128_NO_REF.png)
:::

::: details CMA-128 (128 channels) - actiCHamp Plus (no REF)
![CMA-128 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_128_Channel/actiCHamp_Plus/CMA-128_NO_REF.png)
:::

::: details CMA-160 (160 channels) - BrainAmp (with REF)
![CMA-160 REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_160_Channel/BrainAmp/CMA-160_REF.png)
:::

::: details CMA-160 (160 channels) - actiCHamp (no REF)
![CMA-160 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_160_Channel/actiCHamp/CMA-160_NO_REF.png)
:::

::: details CMA-160 (160 channels) - actiCHamp Plus (no REF)
![CMA-160 NO REF layout](/layouts/acticap/actiCAP_electrode_cap_CMA/actiCAP_UOL_160_Channel/actiCHamp_Plus/CMA-160_NO_REF.png)
:::

#### actiCAP Slim (AS)

::: details AS-32 (32 channels) - with REF
![AS-32 REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-32/AS-32_REF.png)
:::

::: details AS-32 (32 channels) - no REF
![AS-32 NO REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-32/AS-32_NO_REF.png)
:::

::: details AS-64 (64 channels) - with REF
![AS-64 REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-64/AS-64_REF.png)
:::

::: details AS-64 (64 channels) - no REF
![AS-64 NO REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-64/AS-64_NO_REF.png)
:::

::: details AS-96 (96 channels) - with REF
![AS-96 REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-96/AS-96_REF.png)
:::

::: details AS-96 (96 channels) - no REF
![AS-96 NO REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-96/AS-96_NO_REF.png)
:::

::: details AS-128 (128 channels) - with REF
![AS-128 REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-128/AS-128_REF.png)
:::

::: details AS-128 (128 channels) - no REF
![AS-128 NO REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-128/AS-128_NO_REF.png)
:::

::: details AS-160 (160 channels) - with REF
![AS-160 REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-160/AS-160_REF.png)
:::

::: details AS-160 (160 channels) - no REF
![AS-160 NO REF layout](/layouts/acticap/actiCAP_slim_for_actiCHamp_AS/AS-160/AS-160_NO_REF.png)
:::

#### actiCap for LiveAmp

::: details CLA-32 (32 channels) - actiCap for LiveAmp
![CLA-32 layout](/layouts/acticap/actiCap_for_LiveAmp_32/CLA-32.png)
:::

#### actiCap Slim (ASP) — for BrainAmpDC

::: details CACS-32 (32 channels) - for BrainAmpDC (with REF)
![CACS-32 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_BrainAmpDC_ASP/CACS-32/CACS-32_REF.png)
:::

::: details CACS-64 (64 channels) - for BrainAmpDC (with REF)
![CACS-64 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_BrainAmpDC_ASP/CACS-64/CACS-64_REF.png)
:::

::: details CACS-96 (96 channels) - for BrainAmpDC (with REF)
![CACS-96 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_BrainAmpDC_ASP/CACS-96/CACS-96_REF.png)
:::

::: details CACS-128 (128 channels) - for BrainAmpDC (with REF)
![CACS-128 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_BrainAmpDC_ASP/CACS-128/CACS-128_REF.png)
:::

#### actiCap Slim (ASP) — for LiveAmp

::: details CACS-32 (32 channels) - for LiveAmp (with REF)
![CACS-32 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_LiveAmp_ASP/CACS-32/CACS-32_REF.png)
:::

::: details CACS-64 (64 channels) - for LiveAmp (with REF)
![CACS-64 REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_LiveAmp_ASP/CACS-64/CACS-64_REF.png)
:::

#### actiCap Slim (ASP) — for actiCHamp Plus

::: details CACS-32 (32 channels) - for actiCHamp Plus (no REF)
![CACS-32 NO REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_actiChamp_Plus_ASP/CACS-32/CACS-32_NO_REF.png)
:::

::: details CACS-64 (64 channels) - for actiCHamp Plus (no REF)
![CACS-64 NO REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_actiChamp_Plus_ASP/CACS-64/CACS-64_NO_REF.png)
:::

::: details CACS-96 (96 channels) - for actiCHamp Plus (no REF)
![CACS-96 NO REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_actiChamp_Plus_ASP/CACS-96/CACS-96_NO_REF.png)
:::

::: details CACS-128 (128 channels) - for actiCHamp Plus (no REF)
![CACS-128 NO REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_actiChamp_Plus_ASP/CACS-128/CACS-128_NO_REF.png)
:::

::: details CACS-160 (160 channels) - for actiCHamp Plus (no REF)
![CACS-160 NO REF layout](/layouts/acticap/actiCap_slim_ASP/actiCap_slim_for_actiChamp_Plus_ASP/CACS-160/CACS-160_NO_REF.png)
:::

#### actiCap Snap — for BrainAmpDC

::: details CACS-21 (21 channels) - for BrainAmpDC (with REF)
![CACS-21 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_BrainAmpDC/CACS_21/CACS-21_REF.png)
:::

::: details CACS-32 (32 channels) - snap for BrainAmpDC (with REF)
![CACS-32 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_BrainAmpDC/CACS-32/CACS-32_REF.png)
:::

::: details CACS-64 (64 channels) - snap for BrainAmpDC (with REF)
![CACS-64 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_BrainAmpDC/CACS-64/CACS-64_REF.png)
:::

::: details CACS-96 (96 channels) - snap for BrainAmpDC (with REF)
![CACS-96 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_BrainAmpDC/CACS-96/CACS-96_REF.png)
:::

::: details CACS-128 (128 channels) - snap for BrainAmpDC (with REF)
![CACS-128 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_BrainAmpDC/CACS-128/CACS-128_REF.png)
:::

#### actiCap Snap — for LiveAmp

::: details CACS-21 (21 channels) - snap for LiveAmp (with REF)
![CACS-21 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_LiveAmp/CACS-21/CACS-21_REF.png)
:::

::: details CACS-32 (32 channels) - snap for LiveAmp (with REF)
![CACS-32 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_LiveAmp/CACS-32/CACS-32_REF.png)
:::

::: details CACS-64 (64 channels) - snap for LiveAmp (with REF)
![CACS-64 REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_LiveAmp/CACS-64/CACS-64_REF.png)
:::

#### actiCap Snap — for actiCHamp Plus

::: details CACS-21 (21 channels) - snap for actiCHamp Plus (no REF)
![CACS-21 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-21/CACS-21_NO_REF.png)
:::

::: details CACS-32 (32 channels) - snap for actiCHamp Plus (no REF)
![CACS-32 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-32/CACS-32_NO_REF.png)
:::

::: details CACS-64 (64 channels) - snap for actiCHamp Plus (no REF)
![CACS-64 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-64/CACS-64_NO_REF.png)
:::

::: details CACS-96 (96 channels) - snap for actiCHamp Plus (no REF)
![CACS-96 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-96/CACS-96_NO_REF.png)
:::

::: details CACS-128 (128 channels) - snap for actiCHamp Plus (no REF)
![CACS-128 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-128/CACS-128_NO_REF.png)
:::

::: details CACS-160 (160 channels) - snap for actiCHamp Plus (no REF)
![CACS-160 NO REF layout](/layouts/acticap/actiCap_snap_CACS_CAS_GACS/actiCap_slim_for_actiChamp_Plus/CACS-160/CACS-160_NO_REF.png)
:::

#### actiCap Snap for LiveAmp (CLACS)

::: details CLACS-32 (32 channels) - actiCap snap for LiveAmp
![CLACS-32 layout](/layouts/acticap/actiCap_snap_for_LiveAmp_CLACS/CLACS-32.png)
:::

### EasyCap

EasyCap montage configurations (M-series):

::: details easycapM1
![easycapM1 layout](/layouts/easycap/easycapM1.png)
:::

::: details easycapM3
![easycapM3 layout](/layouts/easycap/easycapM3.png)
:::

::: details easycapM7
![easycapM7 layout](/layouts/easycap/easycapM7.png)
:::

::: details easycapM10
![easycapM10 layout](/layouts/easycap/easycapM10.png)
:::

::: details easycapM11
![easycapM11 layout](/layouts/easycap/easycapM11.png)
:::

::: details easycapM14
![easycapM14 layout](/layouts/easycap/easycapM14.png)
:::

::: details easycapM15
![easycapM15 layout](/layouts/easycap/easycapM15.png)
:::

::: details easycapM16
![easycapM16 layout](/layouts/easycap/easycapM16.png)
:::

::: details easycapM17
![easycapM17 layout](/layouts/easycap/easycapM17.png)
:::

::: details easycapM20
![easycapM20 layout](/layouts/easycap/easycapM20.png)
:::

::: details easycapM22
![easycapM22 layout](/layouts/easycap/easycapM22.png)
:::

::: details easycapM23
![easycapM23 layout](/layouts/easycap/easycapM23.png)
:::

::: details easycapM24
![easycapM24 layout](/layouts/easycap/easycapM24.png)
:::

::: details easycapM25
![easycapM25 layout](/layouts/easycap/easycapM25.png)
:::

::: details easycapM64
![easycapM64 layout](/layouts/easycap/easycapM64.png)
:::

## Coordinate Systems

### Polar Coordinates (Source of Truth)

- **Incidence (inc)**: Vertical angle in degrees from the vertex (Cz = 0°)
  - 0° = vertex (top of head)
  - 90° = equator (ear level)
  - Values > 90° indicate positions below the equator (e.g., neck electrodes)

- **Azimuth (azi)**: Horizontal angle in degrees
  - 0° = right hemisphere
  - ±90° = midline (front/back)
  - ±180° = left hemisphere

### Cartesian Coordinates

Cartesian coordinates are computed on-demand from polar coordinates:

- **2D (x2, y2)**: Topographic projection for 2D plotting
  - By default, `preserve_radial_distance=true` maintains anatomical accuracy
  - Electrodes with inc > 90° appear outside the unit circle

- **3D (x3, y3, z3)**: Spherical projection for 3D visualization
  - Normalized to unit sphere

## Related Functions

Layout management and manipulation:

```@docs
EegFun.read_layout
EegFun.validate_layout
EegFun.polar_to_cartesian_xy!
EegFun.polar_to_cartesian_xyz!
EegFun.get_neighbours_xy!
EegFun.get_neighbours_xyz!
EegFun.rename_channel!
```

Layout visualization:

```@docs
EegFun.plot_layout_2d
EegFun.plot_layout_3d
```

## See Also

- [Layout plotting demo](../demos/plotting/plot_layout.md)
- [Layout types](../reference/types.md#Layout-Types)
- [Data structures](data-structures.md)
