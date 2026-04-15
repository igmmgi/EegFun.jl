# Code Examples & Demos

Browse interactive code examples demonstrating EegFun.jl's capabilities.

## Data Import

- [BioSemi Import](import/biosemi_import.md) - Loading BioSemi BDF files
- [BrainVision Import](import/brainvision_import.md) - Loading BrainVision format files
- [EEGLAB Import](import/eeglab_import.md) - Loading EEGLAB .set files
- [FieldTrip Import](import/fieldtrip_import.md) - Loading FieldTrip .mat files
- [Data Persistence](import/jld2.md) - Saving and loading with JLD2
- [BIDS Export](import/bids_export.md) - Exporting to BIDS-compliant directory structure
- [Data Structures](import/data.md) - Core data types and access functions

## Working with Data

- [Data Manipulation](data/data.md) - Data access, subsetting, and manipulation
- [Data Access](data/data_access.md) - Inspecting data with head, tail, viewer
- [Selection Helpers](data/selection_helpers.md) - Filtering with channels(), times(), epochs()

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
- [Layouts & Neighbours](preprocessing/layouts.md) - Electrode layouts and spatial neighbours
- [Analysis Settings](preprocessing/analysis_settings.md) - Saving and replaying preprocessing recipes

## Artifact Detection & Correction

- [Artifacts Overview](artifacts/artifacts.md) - Artifact types, detection, repair, and rejection strategies
- [Artifact Detection](artifacts/artifact_detection.md) - Automatic epoch rejection

## ERP Analysis

- [ERP Measurements](erp/erp_measurements.md) - Quantifying ERP components
- [Global Field Power](erp/gfp.md) - GFP calculation
- [Grand Average](erp/grand_average.md) - Averaging across participants
- [Jackknife Average](erp/jackknife_average.md) - Leave-one-out averaging
- [Lateralised Readiness Potential](erp/lrp.md) - LRP calculation
- [Realignment](erp/realign.md) - Response-locked realignment
- [Condition Operations](erp/condition_operations.md) - Combining and differencing conditions

## Visualization

- [Plot Artifacts](plotting/plot_artifacts.md) - Artifact visualization
- [Plot Channel Spectrum](plotting/plot_channel_spectrum.md) - Frequency spectra
- [Plot Channel Summary](plotting/plot_channel_summary.md) - Channel statistics
- [Plot Correlation Heatmap](plotting/plot_correlation_heatmap.md) - Channel correlations
- [Plot Epochs](plotting/plot_epochs.md) - Individual epoch visualization
- [Plot ERP](plotting/plot_erp.md) - ERP waveforms
- [Plot ERP Image](plotting/plot_erp_image.md) - ERP images
- [Plot ERP Measurements](plotting/plot_erp_measurements.md) - Measurement visualization
- [Plot Filter](plotting/plot_filter.md) - Filter responses
- [Plot Frequency Spectrum](plotting/plot_frequency_spectrum.md) - SpectrumData power spectra
- [Plot GFP](plotting/plot_gfp.md) - Global Field Power
- [Plot ICA](plotting/plot_ica.md) - ICA component topographies
- [Plot Joint Probability](plotting/plot_joint_probability.md) - Artifact detection
- [Plot Layout](plotting/plot_layout.md) - Electrode layouts
- [Plot Topography](plotting/plot_topography.md) - Topographic maps
- [Plot Triggers](plotting/plot_triggers.md) - Event markers
- [Saving Figures](plotting/save_figures.md) - Exporting figures (PNG, SVG, PDF)

## Interactive / GUI

- [Plot Databrowser](plotting/plot_databrowser.md) - Interactive data browsing
- [Plot ERP Filter GUI](plotting/plot_erp_filter_gui.md) - Interactive filter exploration
- [Plot ERP Measurement GUI](plotting/plot_erp_measurement_gui.md) - Interactive measurement GUI

## Teaching Demos

Interactive demonstrations for some signal processing and EEG concepts.

- [Nyquist Theorem](teaching/signal_processing/signal_example_sampling.md) - Sampling, aliasing, and signal reconstruction (linear vs. sinc)
- [Signal Composition](teaching/signal_processing/signal_example_composition.md) - Building complex waveforms from sine waves, noise, and filtering
- [Dot Product](teaching/signal_processing/signal_example_dotproduct.md) - How the DFT detects frequencies via dot products
- [Convolution](teaching/signal_processing/signal_example_convolution.md) - Sliding kernels including Morlet wavelets
- [Power Spectrum](teaching/signal_processing/signal_example_spectrum.md) - FFT power spectrum, frequency resolution, and spectral leakage
- [Time-Frequency](teaching/signal_processing/signal_example_tf.md) - Comparing Morlet, STFT, and multitaper TF methods
- [Simulate ERP](teaching/erp/simulate_erp.md) - Trial averaging, 1/f noise, and ERP component decomposition
- [Simulate Decoding](teaching/machine_learning/signal_example_decoding.md) - MVPA/SVM classification on synthetic EEG
- [ICA 0: Matrix Math (U = WX)](teaching/ica/signal_example_ica_0.md) - The fundamental equation underlying component analysis
- [ICA 1: What Is a Mixture?](teaching/ica/signal_example_ica_1.md) - Why EEG needs ICA: electrodes record weighted sums
- [ICA 2: Mixing & Unmixing](teaching/ica/signal_example_ica_2.md) - The cocktail party problem: mixing matrices and ICA separation
- [ICA 3: Blind Source Separation](teaching/ica/signal_example_ica_3.md) - Advanced: 3-source rotation geometry, scatter plots, and kurtosis
- [ICA: The Central Limit Theorem](teaching/ica/signal_example_ica_clt.md) - Deep Dive: Why mixing creates Gaussian blobs, and how to unmix them
- [ICA 4: Sphering (Whitening)](teaching/ica/signal_example_ica_4.md) - Preprocessing: morphing correlated blobs into perfect spheres
- [ICA 5: Inside the Black Box](teaching/ica/signal_example_ica_5.md) - The optimization landscape: how gradient ascent recovers the sources
- [ICA 6: Infomax (Information Theory)](teaching/ica/signal_example_ica_infomax.md) - The math of EEGLAB: how a neural network unpacks Gaussians into Maximum Entropy


## Specialized Visualization

- [Plot Decoding](plotting/plot_decoding.md) - MVPA decoding results
- [Plot RSA](plotting/plot_rsa.md) - RSA results (RDM, timecourse, models)
- [Plot Statistics](plotting/plot_statistics.md) - Statistical test results
- [Plot Time-Frequency](plotting/plot_tf.md) - Time-frequency heatmaps

## Workflows

- [Preprocessing Workflow](workflows/preprocessing_workflow.md) - End-to-end preprocessing pipeline
- [Batch Processing](workflows/batch_processing.md) - Loading and grouping multiple files
- [Pipeline Templates](workflows/pipeline_templates.md) - Generating scaffold files

## Time-Frequency Analysis

- [TF Morlet](time_frequency/tf_morlet.md) - Morlet wavelet analysis
- [TF Multitaper](time_frequency/tf_multitaper.md) - Multitaper spectral estimation
- [TF STFT](time_frequency/tf_stft.md) - Short-Time Fourier Transform
- [TF Operations](time_frequency/tf_operations.md) - Condition & channel operations on TF data

## Statistics

- [Statistics](statistics/statistics.md) - Statistical testing
- [TF Statistics](statistics/tf_stats_test.md) - Time-frequency statistical testing
- [Decoding](statistics/decoding.md) - Multivariate pattern analysis (MVPA)
- [RSA](statistics/rsa.md) - Representational Similarity Analysis

## Example Experiments

- [N170 (Face/Body)](experiments/n170.md) - Complete walkthrough: face-sensitive N170 component
- [Visual Attention (Posner Cueing)](experiments/visual-attention.md) - Complete walkthrough: endogenous attention with cluster statistics
