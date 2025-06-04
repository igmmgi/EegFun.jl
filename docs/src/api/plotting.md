# Plotting

This section documents visualization functions for EEG data including ERPs, topographic maps, and spectral plots.

## ERP Plotting

### `plot_erp!(fig, ax, dat::ErpData, channels; kwargs...)`
Add ERP plot to existing figure/axis.

**Arguments:**
- `fig`: Figure object
- `ax`: Axis object  
- `dat::ErpData`: ERP data to plot
- `channels::Vector{Symbol}`: Channels to plot

### `plot_epochs(dat::EpochData, channels; kwargs...)`
Plot epoched EEG data showing individual trials and average.

**Returns:** Figure and Axis objects

### `plot_erp_image(dat::EpochData, channels; kwargs...)`
Create an ERP image plot showing trial-by-trial data.

**Returns:** Figure and Axis objects

## Topographic Plots

### `plot_topoplot!(fig, ax, dat, layout; kwargs...)`
Add topographic plot to existing figure/axis.

### `plot_topoplot(dat, layout; kwargs...)`
Create topographic plot from EEG data.

**Returns:** Figure and Axis objects

### `plot_layout_2d!(fig, ax, layout; kwargs...)`
Plot 2D EEG electrode layout with customizable head shape, electrode points, and labels.

### `create_convex_hull_graham(xpos, ypos, border_size)`
Create a convex hull around 2D points using Graham's Scan algorithm.

## Spectral Analysis Plots

### `plot_selected_spectrum(selected_data, channel; kwargs...)`
Plot power spectrum for selected data and channels.

**Returns:** Figure and Axis objects

### `plot_component_spectrum(ica_result, dat, comp_idx; kwargs...)`
Plot power spectrum of a specific ICA component.

**Returns:** Figure and Axis objects

### `plot_components_spectra(ica_result, dat, comp_indices; kwargs...)`
Plot power spectra for multiple ICA components.

**Returns:** Figure and Axis objects

## Utility Plotting

### `test_plot_eog_detection(dat, xlim, channel, detected)`
Test function for plotting EOG detection results.

**Returns:** Figure and Axis objects 