# Code Examples & Demos

Browse interactive code examples demonstrating EegFun.jl's capabilities.

## Data Import

- [BioSemi Import](import/biosemi_import.md) - Loading BioSemi BDF files
- [BrainVision Import](import/brainvision_import.md) - Loading BrainVision format files
- [EEGLAB Import](import/eeglab_import.md) - Loading EEGLAB .set files
- [FieldTrip Import](import/fieldtrip_import.md) - Loading FieldTrip .mat files
- [Data Creation](import/data.md) - Creating synthetic test data
- [Data Persistence](import/jld2.md) - Saving and loading with JLD2

## Data Processing

- [Filter](preprocessing/filter.md) - High-pass and low-pass filtering
- [Baseline Correction](preprocessing/baseline.md) - Removing baseline activity
- [Rereferencing](preprocessing/rereference.md) - Changing reference scheme
- [Resampling](preprocessing/resample.md) - Changing sampling rate
- [Mirror Padding](preprocessing/mirror.md) - Edge artifact reduction
- [Triggers](preprocessing/triggers.md) - Working with event markers
- [Epoch Extraction](preprocessing/epochs.md) - Extracting epochs from continuous data
- [Channel Operations](preprocessing/channel_operations.md) - Channel arithmetic
- [Channel Metrics](preprocessing/channel_metrics.md) - Computing channel statistics
- [Channel Repair](preprocessing/channel_repair.md) - Interpolating bad channels
- [Channel Summary](preprocessing/channel_summary.md) - Summarizing channel data
- [ICA](preprocessing/ica.md) - Independent Component Analysis

## Artifact Detection & Correction

- [Artifacts](artifacts/artifacts.md) - Detecting and managing artifacts
- [Artifact Detection](artifacts/artifact_detection.md) - Automatic epoch rejection

## ERP Analysis

- [ERP Measurements](erp/erp_measurements.md) - Quantifying ERP components
- [Global Field Power](erp/gfp.md) - GFP calculation
- [Grand Average](erp/grand_average.md) - Averaging across participants
- [Jackknife Average](erp/jackknife_average.md) - Leave-one-out averaging
- [Lateralised Readiness Potential](erp/lrp.md) - LRP calculation
- [Realignment](erp/realign.md) - Response-locked realignment
- [Condition Operations](erp/condition_operations.md) - Combining and differencing conditions

## Time-Frequency Analysis

- [TF Morlet](time_frequency/tf_morlet.md) - Morlet wavelet analysis
- [TF Multitaper](time_frequency/tf_multitaper.md) - Multitaper spectral estimation
- [TF STFT](time_frequency/tf_stft.md) - Short-Time Fourier Transform

## Statistics

- [Statistics](statistics/statistics.md) - Statistical testing
- [Decoding](statistics/decoding.md) - Multivariate pattern analysis (MVPA)
- [RSA](statistics/rsa.md) - Representational Similarity Analysis

## Visualization

- [Plot Databrowser](plotting/plot_databrowser.md) - Interactive data browsing
- [Plot ERP](plotting/plot_erp.md) - ERP waveforms
- [Plot Epochs](plotting/plot_epochs.md) - Individual epoch visualization
- [Plot Topography](plotting/plot_topography.md) - Topographic maps
- [Plot ERP Image](plotting/plot_erp_image.md) - ERP images
- [Plot ERP Measurements](plotting/plot_erp_measurements.md) - Measurement visualization
- [Plot ERP Measurement GUI](plotting/plot_erp_measurement_gui.md) - Interactive measurement GUI
- [Plot Channel Spectrum](plotting/plot_channel_spectrum.md) - Frequency spectra
- [Plot Channel Summary](plotting/plot_channel_summary.md) - Channel statistics
- [Plot Correlation Heatmap](plotting/plot_correlation_heatmap.md) - Channel correlations
- [Plot Filter](plotting/plot_filter.md) - Filter responses
- [Plot Joint Probability](plotting/plot_joint_probability.md) - Artifact detection
- [Plot Layout](plotting/plot_layout.md) - Electrode layouts
- [Plot Triggers](plotting/plot_triggers.md) - Event markers
- [Plot Artifacts](plotting/plot_artifacts.md) - Artifact visualization

## Workflows

- [Preprocessing Workflow](workflows/preprocessing_workflow.md) - End-to-end preprocessing pipeline
