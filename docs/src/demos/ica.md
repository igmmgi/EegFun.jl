# ICA Analysis

Complete ICA workflow including component identification, artifact removal, and visualization for continuous and epoched data.

## Overview

Demonstrates Complete ICA workflow including component identification, artifact removal, and visualization for continuous and epoched data.

## Source Code

```julia
using EegFun

# read raw data
dat = EegFun.read_raw_data("./resources/data/example1.bdf");

# read and preprate layout file
layout_file = EegFun.read_layout("./resources/layouts/biosemi/biosemi72.csv");
EegFun.polar_to_cartesian_xy!(layout_file)

# create EegFun data structure (EegFun.ContinuousData)
dat = EegFun.create_eeg_dataframe(dat, layout_file);

# Some minimal preprocessing (average reference, highpass filter, and detect extreme values)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1)
EegFun.is_extreme_value!(dat, 200);

# Calculate EOG signals
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:Fp1, :Fp2]),
    channel_selection2 = EegFun.channels([:IO1, :IO2]),
    channel_out = :vEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.channel_difference!(
    dat,
    channel_selection1 = EegFun.channels([:F9]),
    channel_selection2 = EegFun.channels([:F10]),
    channel_out = :hEOG,
); # vertical EOG = mean(Fp1, Fp2) - mean(IO1, I02)
EegFun.detect_eog_onsets!(dat, 50, :vEOG, :is_vEOG)
EegFun.detect_eog_onsets!(dat, 30, :hEOG, :is_hEOG)

# ICA on continuous data
ica_result_infomax = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_200), percentage_of_data = 100)
ica_result_infomax_extended = EegFun.run_ica(
    dat;
    sample_selection = EegFun.samples_not(:is_extreme_value_200),
    percentage_of_data = 20,
    algorithm = :infomax_extended,
)

# ICA type plots
EegFun.plot_topography(ica_result_infomax, component_selection = EegFun.components(1:4));
EegFun.plot_topography(ica_result_infomax_extended, component_selection = EegFun.components(1:4));
EegFun.plot_ica_component_activation(dat, ica_result_infomax)
EegFun.plot_ica_component_spectrum(dat, ica_result_infomax, component_selection = EegFun.components(1:70))
EegFun.plot_databrowser(dat, ica_result_infomax)

# identify components
component_artifacts, component_metrics =
    EegFun.identify_components(dat, ica_result_infomax, sample_selection = EegFun.samples_not(:is_extreme_value_200))

# or individually
eog_comps, eog_comps_metrics_df =
    EegFun.identify_eog_components(dat, ica_result_infomax_extended, sample_selection = EegFun.samples_not(:is_extreme_value_200))
ecg_comps, ecg_comps_metrics_df =
    EegFun.identify_ecg_components(dat, ica_result_infomax_extended, sample_selection = EegFun.samples_not(:is_extreme_value_200))
line_noise_comps, line_noise_comps_metrics_df = EegFun.identify_line_noise_components(dat, ica_result_infomax_extended)
channel_noise_comps, channel_noise_comps_metrics_df = EegFun.identify_spatial_kurtosis_components(ica_result_infomax_extended)


# Get all identified component artifacts
all_comps = EegFun.get_all_ica_components(component_artifacts)
dat_ica_removed, ica_result_updated =
    EegFun.remove_ica_components(dat, ica_result_infomax_extended, component_selection = EegFun.components(all_comps))

# Reconstruct for sanity check
dat_ica_reconstructed, ica_result_restored =
    EegFun.restore_ica_components(dat_ica_removed, ica_result_updated, component_selection = EegFun.components(all_comps))

# Original should = reconstructed
EegFun.channel_data(dat) ≈ EegFun.channel_data(dat_ica_reconstructed)

# Plot component features
fig, ax = EegFun.plot_eog_component_features(eog_comps, eog_comps_metrics_df) # TODO: points sizes
fig, ax = EegFun.plot_ecg_component_features_(ecg_comps, ecg_comps_metrics_df)
fig, ax = EegFun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
fig, ax = EegFun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)


#################################
# Epoched DataFrameEeg
#################################
# Create some epoched data
epoch_cfg = [
    EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[2]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.2, 1.0)  # -200 to 1000 ms

# ICA on epoched data
ica_result_infomax = EegFun.run_ica(epochs[1]; sample_selection = EegFun.samples_not(:is_extreme_value_200))
ica_result_infomax_extended = EegFun.run_ica(epochs; sample_selection = EegFun.samples_not(:is_extreme_value_200))


# ICA type plots
EegFun.plot_topography(ica_result_infomax, component_selection = EegFun.components(1:4));
EegFun.plot_topography(ica_result_infomax_extended, component_selection = EegFun.components(1:4));
EegFun.plot_ica_component_activation(dat, ica_result_infomax_extended)
EegFun.plot_ica_component_spectrum(dat, ica_result_infomax, component_selection = EegFun.components(1:70))
EegFun.plot_databrowser(dat, ica_result_infomax)

# identify components
component_artifacts, component_metrics =
    EegFun.identify_components(dat, ica_result_infomax_extended, sample_selection = EegFun.samples_not(:is_extreme_value_200))

# Get all identified component artifacts
all_comps = EegFun.get_all_ica_components(component_artifacts)
dat_ica_removed, ica_result_updated = EegFun.remove_ica_components(dat, ica_result, component_selection = EegFun.components(all_comps))





# Reconstruct for sanity check
dat_ica_reconstructed, ica_result_restored =
    EegFun.restore_ica_components(dat_ica_removed, ica_result_updated, component_selection = EegFun.components(all_comps))

# Original should = reconstructed
EegFun.channel_data(dat) ≈ EegFun.channel_data(dat_ica_reconstructed)

# Plot component features
fig, ax = EegFun.plot_eog_component_features(eog_comps, eog_comps_metrics_df)
fig, ax = EegFun.plot_ecg_component_features_(ecg_comps, ecg_comps_metrics_df)
fig, ax = EegFun.plot_line_noise_components(line_noise_comps, line_noise_comps_metrics_df)
fig, ax = EegFun.plot_spatial_kurtosis_components(channel_noise_comps, channel_noise_comps_metrics_df)


GLMakie.closeall()
```

## See Also

- [API Reference](../reference/index.md)
