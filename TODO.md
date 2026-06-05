# EegFun.jl Plotting Improvements TODO

- [x] **1. Standardize UI and Theme Font Sizes (`ui_fontsize` / `theme_fontsize`)**
  - Unify `:ui_fontsize` (or `:theme_fontsize`) across all interactive plots (`plot_erp`, `plot_epochs`, `plot_topography`, `plot_tf`, `plot_databrowser`).

- [x] **2. Expose Keyboard Navigation Step Sizes (`scroll_step`, `zoom_step`)**
  - Add `:zoom_step` and `:scroll_step` to `PLOT_ERP_KWARGS`, `PLOT_EPOCHS_KWARGS`, and `PLOT_TOPOGRAPHY_KWARGS` to configure arrow key scaling/scrolling.

- [x] **3. Unified Interactive Selection Styling (`selection_color`)**
  - Add `:selection_color` and `:selection_alpha` to `PLOT_ERP_KWARGS` and `PLOT_EPOCHS_KWARGS` for interactive time-window selections.



- [x] **5. Propagate `figure_padding` Everywhere**
  - Add `:figure_padding` to `PLOT_TF_KWARGS` and `PLOT_TOPOGRAPHY_KWARGS`.

- [x] **6. Scale Indicators for ERPs and Epochs**
  - Add `:show_scale_indicator` (with value and position) to `plot_erp` and `plot_epochs` as an alternative to y-axis ticks.
