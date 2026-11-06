This demo walks through EegFun's interactive GUIs for inspecting Independent Component Analysis (ICA) decompositions, component spectra, and live artifact removal.

### Included Interactive Tools

| Function | Purpose |
| --- | --- |
| `plot_topography(ica)` | Paginated component scalp-map browser with component selection textbox, scale checkboxes, and double-click detail popup |
| `plot_ica_component_activation` | Synchronized scalp-map + scrollable time-course viewer with artifact category checkboxes (`EOG`, `ECG`, line/channel noise) |
| `plot_ica_component_spectrum` | Interactive component PSD viewer with live Linear / `Log10` / `dB` axis dropdown menus |
| `plot_databrowser(dat, ica)` | Full continuous/epoched browser with live per-component toggling and cleaned-trace (`Data - ICA`) overlay |

### Related GUIs

For interactive educational simulations of ICA math and geometry (`signal_example_ica_*`), see `tutorials/gui/ica_demos.jl` and the Teaching Demos section.

