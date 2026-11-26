# Top-level module definition and file includes for EegFun.jl.

module EegFun

# Core dependencies
using BiosemiDataFormat
using BrainVisionDataFormat
using EuropeanDataFormat
using ExtensibleDataFormat
using FunctionalImageFormat
using CSV
using DataFrames
using DataDeps
using JLD2
using MAT

# Signal processing & GPU
using DSP
using LinearAlgebra
using SpecialFunctions
using FFTW
using KernelAbstractions

# Plotting and visualization
using Makie
using GeometryBasics
using ScatteredInterpolation

# Utilities
using AnovaFun
using Base64
using Dates
using Distributions
using DocStringExtensions
using Logging
using LoggingExtras
using NativeFileDialog
using OrderedCollections
using PrettyTables
using Printf
using ProgressMeter
using Random
using SparseArrays
using Statistics
using StatsAPI
using StatsBase
using TOML

# Machine learning (for decoding/MVPA)
using LIBSVM


# Public API declarations
include("public_api.jl")

# Core types
include("types/core.jl")
include("types/pipeline.jl")
include("types/statistics.jl")
include("types/statistics_tf.jl")
include("types/decoding.jl")
include("types/rsa.jl")

# Utilities and EegData operations
include("utils/error.jl")
include("utils/gpu.jl")
include("utils/create_data.jl")
include("utils/batch.jl")
include("utils/data/columns.jl")
include("utils/data/access.jl")
include("utils/data/frames.jl")
include("utils/data/math.jl")
include("utils/data/subset.jl")
include("utils/data/subset_files.jl")
include("utils/data/rename.jl")
include("utils/data/from_raw.jl")
include("utils/data/base_interface.jl")
include("utils/files.jl")
include("utils/logging.jl")
include("utils/misc.jl")
include("utils/viewer.jl")
include("utils/print.jl")
include("utils/tables.jl")

# Layouts, I/O, and configuration
include("layouts/layout.jl")
include("io/read_data.jl")
include("io/eeglab.jl")
include("io/fieldtrip.jl")
include("io/xdf.jl")
include("io/bids_export.jl")
include("io/bids_import.jl")
include("config/config.jl")

# Analysis — artifact detection
include("analysis/processing/artifact_detection.jl")
include("analysis/processing/artifact_detection/eog_onsets.jl")
include("analysis/processing/artifact_detection/extreme.jl")
include("analysis/processing/artifact_detection/step.jl")
include("analysis/processing/artifact_detection/flatline.jl")
include("analysis/processing/artifact_detection/peak_to_peak.jl")
include("analysis/processing/artifact_detection/bridged.jl")
include("analysis/processing/artifact_detection/epoch_rejection_types.jl")
include("analysis/processing/artifact_detection/epoch_detection.jl")
include("analysis/processing/artifact_detection/repair.jl")
include("analysis/processing/artifact_detection/subset_bad_data.jl")
include("analysis/processing/artifact_detection/autoreject.jl")

# Analysis — signal, channel, epoch, condition, and ICA processing
include("analysis/processing/signal/baseline.jl")
include("analysis/processing/signal/detrend.jl")
include("analysis/processing/signal/filter.jl")
include("analysis/processing/signal/cleanline.jl")
include("analysis/processing/signal/mirror.jl")
include("analysis/processing/signal/resample.jl")
include("analysis/processing/channels/channel_average.jl")
include("analysis/processing/channels/channel_delete.jl")
include("analysis/processing/channels/channel_difference.jl")
include("analysis/processing/channels/channel_repair.jl")
include("analysis/processing/channels/channel_summary.jl")
include("analysis/processing/channels/channel_metrics.jl")
include("analysis/processing/channels/rereference.jl")
include("analysis/processing/channels/csd.jl")
include("analysis/processing/epochs/triggers.jl")
include("analysis/processing/epochs/epochs.jl")
include("analysis/processing/epochs/apply.jl")
include("analysis/processing/conditions/condition_combine.jl")
include("analysis/processing/conditions/condition_difference.jl")
include("analysis/processing/conditions/condition_average.jl")
include("analysis/processing/conditions/grand_average.jl")
include("analysis/processing/ica/api.jl")
include("analysis/processing/ica/infomax.jl")
include("analysis/processing/ica/picard.jl")
include("analysis/processing/ica/amica.jl")
include("analysis/processing/ica/components.jl")
include("analysis/processing/ica/artifact_id.jl")

# Analysis — ERPs, RIDE, statistics, time-frequency, decoding, and RSA
include("analysis/erps/jackknife_average.jl")
include("analysis/erps/erp_measurements.jl")
include("analysis/erps/gfp.jl")
include("analysis/erps/lrp.jl")
include("analysis/erps/realign.jl")
include("analysis/ride/utils.jl")
include("analysis/ride/woody.jl")
include("analysis/ride/iterator.jl")
include("analysis/ride/decompose.jl")
include("analysis/ride/ride.jl")
include("analysis/statistics/core.jl")
include("analysis/statistics/statistics_tf_helpers.jl")
include("analysis/statistics/thresholding.jl")
include("analysis/statistics/clustering.jl")
include("analysis/statistics/permutations.jl")
include("analysis/statistics/inference.jl")
include("analysis/statistics/statistics.jl")
include("analysis/statistics/lmm_stats.jl")
include("analysis/statistics/statsapi.jl")
include("analysis/time_frequency/tf_utils.jl")
include("analysis/time_frequency/tf_morlet.jl")
include("analysis/time_frequency/tf_stft.jl")
include("analysis/time_frequency/tf_multitaper.jl")
include("analysis/time_frequency/tf_itc.jl")
include("analysis/time_frequency/freq_spectrum.jl")
include("analysis/time_frequency/baseline.jl")
include("analysis/decoding/liblinear.jl")
include("analysis/decoding/decoding.jl")
include("analysis/decoding/preparation.jl")
include("analysis/decoding/decoding_statistics.jl")
include("analysis/rsa/rsa_models.jl")
include("analysis/rsa/rsa.jl")
include("analysis/rsa/rsa_noise_ceiling.jl")
include("analysis/rsa/rsa_crossvalidation.jl")

# Automated pipelines and reporting
include("pipelines/utils/utils.jl")
include("pipelines/utils/repair_info.jl")
include("pipelines/utils/repair.jl")
include("pipelines/utils/summaries.jl")
include("pipelines/utils/template.jl")
include("pipelines/pipeline.jl")
include("pipelines/report.jl")

# Plotting and interactive GUIs
include("plots/utils/plot_misc.jl")
include("plots/utils/themes.jl")
include("plots/utils/plot_ui_utils.jl")
include("plots/utils/layout_system.jl")
include("plots/utils/shared_interactivity.jl")
include("plots/processing/plot_artifacts.jl")
include("plots/processing/plot_channel_summary.jl")
include("plots/processing/plot_correlation_heatmap.jl")
include("plots/processing/plot_databrowser.jl")
include("plots/processing/plot_ica.jl")
include("plots/processing/plot_ica_artifacts.jl")
include("plots/processing/plot_joint_probability.jl")
include("plots/processing/plot_layout.jl")
include("plots/processing/plot_triggers.jl")
include("plots/processing/plot_filter.jl")
include("plots/erps/plot_epochs.jl")
include("plots/erps/plot_epochs_rejection.jl")
include("plots/erps/plot_erp.jl")
include("plots/erps/plot_erp_image.jl")
include("plots/erps/plot_erp_measurements.jl")
include("plots/erps/plot_erp_measurement_gui.jl")
include("plots/erps/plot_erp_filter_gui.jl")
include("plots/erps/plot_gfp.jl")
include("plots/erps/plot_erp_stats.jl")
include("plots/erps/plot_topography_stats.jl")
include("plots/erps/plot_topography.jl")
include("plots/erps/plot_topography_3d.jl")
include("plots/time_frequency/plot_tf.jl")
include("plots/time_frequency/plot_tf_stats.jl")
include("plots/time_frequency/plot_power_spectrum.jl")
include("plots/time_frequency/plot_topography_tf.jl")
include("plots/decoding/plot_decoding.jl")
include("plots/rsa/plot_rsa.jl")
include("plots/plot_gui.jl")
include("plots/plot_gui_helpers.jl")
include("plots/utils/help_system.jl")
include("plots/utils/makie_recipes.jl")

# Interactive educational demos
include("demos/demo_utils.jl")
include("demos/signal_processing/signal_example_sampling.jl")
include("demos/signal_processing/signal_example_composition.jl")
include("demos/signal_processing/signal_example_dotproduct.jl")
include("demos/signal_processing/signal_example_convolution.jl")
include("demos/signal_processing/signal_example_tf.jl")
include("demos/ica/signal_example_ica_mixture.jl")
include("demos/ica/signal_example_ica_separation.jl")
include("demos/ica/signal_example_ica_geometry.jl")
include("demos/ica/signal_example_ica_clt.jl")
include("demos/ica/signal_example_ica_sphering.jl")
include("demos/ica/signal_example_ica_optimization.jl")
include("demos/erp/simulate_erp.jl")
include("demos/machine_learning/signal_example_decoding.jl")
include("demos/eog/signal_example_eog.jl")
include("demos/emg/signal_example_emg.jl")
include("demos/ecg/signal_example_ecg_einthoven.jl")
include("demos/ecg/signal_example_ecg_hrv.jl")
include("demos/eda/signal_example_eda_scr.jl")
include("demos/eda/signal_example_eda_habituation.jl")

# Precompilation
include("precompile.jl")

function __init__()
    # Note: the EegFun plot theme is applied per-figure (see `_eegfun_figure`), not globally,
    # so loading EegFun never changes the user's Makie theme.

    register(
        DataDep(
            "TutorialDataSets",
            "EegFun Tutorial Datasets",
            "https://zenodo.org/records/21770157/files/TutorialDataSets_EegFun.zip?download=1",
            "9c49c03b7185ba4a5f4aae6d1091ddbdafbda9f37ec10c03338e3dd1665d79ef",
            post_fetch_method = unpack,
        ),
    )
end

end # module EegFun
