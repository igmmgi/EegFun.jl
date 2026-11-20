```@meta
CollapsedDocStrings = true
```

# Plotting, GUIs & Interactive Demos

All public visualization functions, plot layouts/themes, interactive inspection GUIs, and teaching demos in `EegFun.jl`, organized by subsystem.

## Index

```@index
Pages = ["plotting.md"]
```

## Themes, Plot Layouts & GUI Launcher

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/plots/utils/themes.jl",
    "src/plots/utils/layout_system.jl",
    "src/plots/plot_gui.jl",
]
```

## Preprocessing, Artifact, Layout & ICA Plots

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/plots/processing/plot_databrowser.jl",
    "src/plots/processing/plot_layout.jl",
    "src/plots/processing/plot_triggers.jl",
    "src/plots/processing/plot_filter.jl",
    "src/plots/processing/plot_artifacts.jl",
    "src/plots/processing/plot_channel_summary.jl",
    "src/plots/processing/plot_correlation_heatmap.jl",
    "src/plots/processing/plot_joint_probability.jl",
    "src/plots/processing/plot_ica.jl",
    "src/plots/processing/plot_ica_artifacts.jl",
]
```

## ERP, Epoch & Topography Plots

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/plots/erps/plot_epochs.jl",
    "src/plots/erps/plot_epochs_rejection.jl",
    "src/plots/erps/plot_erp.jl",
    "src/plots/erps/plot_erp_image.jl",
    "src/plots/erps/plot_erp_measurements.jl",
    "src/plots/erps/plot_erp_filter_gui.jl",
    "src/plots/erps/plot_erp_measurement_gui.jl",
    "src/plots/erps/plot_erp_stats.jl",
    "src/plots/erps/plot_gfp.jl",
    "src/plots/erps/plot_topography.jl",
    "src/plots/erps/plot_topography_3d.jl",
    "src/plots/erps/plot_topography_stats.jl",
]
```

## Time-Frequency & Power Spectrum Plots

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/plots/time_frequency/plot_power_spectrum.jl",
    "src/plots/time_frequency/plot_tf.jl",
    "src/plots/time_frequency/plot_tf_stats.jl",
    "src/plots/time_frequency/plot_topography_tf.jl",
]
```

## Decoding & RSA Plots

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/plots/decoding/plot_decoding.jl",
    "src/plots/rsa/plot_rsa.jl",
]
```

## Interactive Teaching Demos

```@autodocs
Modules = [EegFun]
Private = false
Order = [:function, :macro]
Pages = [
    "src/demos/signal_processing/signal_example_sampling.jl",
    "src/demos/signal_processing/signal_example_composition.jl",
    "src/demos/signal_processing/signal_example_dotproduct.jl",
    "src/demos/signal_processing/signal_example_convolution.jl",
    "src/demos/signal_processing/signal_example_tf.jl",
    "src/demos/ica/signal_example_ica_clt.jl",
    "src/demos/ica/signal_example_ica_mixture.jl",
    "src/demos/ica/signal_example_ica_separation.jl",
    "src/demos/ica/signal_example_ica_geometry.jl",
    "src/demos/ica/signal_example_ica_sphering.jl",
    "src/demos/ica/signal_example_ica_optimization.jl",
    "src/demos/erp/simulate_erp.jl",
    "src/demos/machine_learning/signal_example_decoding.jl",
    "src/demos/eog/signal_example_eog.jl",
    "src/demos/emg/signal_example_emg.jl",
    "src/demos/ecg/signal_example_ecg_einthoven.jl",
    "src/demos/ecg/signal_example_ecg_hrv.jl",
    "src/demos/eda/signal_example_eda_scr.jl",
    "src/demos/eda/signal_example_eda_habituation.jl",
]
```
