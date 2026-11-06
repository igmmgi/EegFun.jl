# ICA Inspection GUIs

This demo walks through EegFun's interactive GUIs for inspecting Independent Component Analysis (ICA) decompositions, component spectra, and live artifact removal.

## Included Interactive Tools

| Function | Purpose |
| --- | --- |
| `plot_topography(ica)` | Paginated component scalp-map browser with component selection textbox, scale checkboxes, and double-click detail popup |
| `plot_ica_component_activation` | Synchronized scalp-map + scrollable time-course viewer with artifact category checkboxes (`EOG`, `ECG`, line/channel noise) |
| `plot_ica_component_spectrum` | Interactive component PSD viewer with live Linear / `Log10` / `dB` axis dropdown menus |
| `plot_databrowser(dat, ica)` | Full continuous/epoched browser with live per-component toggling and cleaned-trace (`Data - ICA`) overlay |

## Related GUIs

For interactive educational simulations of ICA math and geometry (`signal_example_ica_*`), see `tutorials/gui/ica_demos.jl` and the Teaching Demos section.



## Code Examples

::: details Show Code

```julia
# Demo: Interactive ICA Inspection GUIs
# Short walkthrough for testing EegFun's interactive ICA viewers:
# component topography browser, component activation viewer, component spectrum GUI,
# and live ICA component toggling inside the databrowser.

using EegFun
using GLMakie

#######################################################################
# PREPARE DATA & RUN ICA
#######################################################################

raw = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(raw, layout)

EegFun.highpass_filter!(dat, 1.0)
EegFun.rereference!(dat, :avg)
EegFun.is_extreme_value!(dat, 200)

ica = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_200), percentage_of_data = 50)
component_artifacts, _ = EegFun.identify_components(dat, ica; sample_selection = EegFun.samples_not(:is_extreme_value_200))

#######################################################################
# 1. INTERACTIVE ICA TOPOGRAPHY BROWSER (plot_topography)
#######################################################################

EegFun.plot_topography(ica, label_plot = false, point_plot = false)

#######################################################################
# 2. ICA COMPONENT ACTIVATION VIEWER (plot_ica_component_activation)
#######################################################################

# Controls to test:
# - Component paging buttons, selection textbox + "Apply"
# - Time position slider & amplitude scale controls (or arrow keys)
# - "Additional Channels" popup menu (overlay EOG / trigger / EEG channels)
# - Artifact category checkboxes (EOG, ECG, line noise, channel noise) when `artifacts` is passed
EegFun.plot_ica_component_activation(dat, ica, artifacts = component_artifacts)

#######################################################################
# 3. INTERACTIVE ICA COMPONENT SPECTRUM (plot_ica_component_spectrum)
#######################################################################

# Controls to test:
# - X-Axis menu: Linear (Hz) vs Log10 (Hz)
# - Y-Axis menu: Linear (μV²/Hz), Log10 (μV²/Hz), Decibels (dB)
EegFun.plot_ica_component_spectrum(dat, ica, component_selection = EegFun.components(1:10))

#######################################################################
# 4. DATABROWSER WITH LIVE ICA REMOVAL PREVIEW (plot_databrowser)
#######################################################################

# Toggle individual ICA components on/off and enable cleaned trace overlay (Data - ICA)
EegFun.plot_databrowser(dat, ica)
```

:::

## See Also

- [API Reference](../../reference/index.md)
