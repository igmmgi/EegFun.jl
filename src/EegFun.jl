__precompile__(true)

module EegFun

# Core dependencies
using AnovaFun
using BiosemiDataFormat
using BrainVisionDataFormat
using CategoricalArrays
using CSV
using DataFrames
using JLD2
using MAT

# Signal processing
using DSP
using LinearAlgebra
using FFTW

# Plotting and visualization
using Makie
using CairoMakie
using GLMakie
using NativeFileDialog
using ScatteredInterpolation

# Utilities
using Dates
using Distributions
using DocStringExtensions
using Logging
using LoggingExtras
using OrderedCollections
using PrettyTables
using Printf
using ProgressMeter
using Random
using SparseArrays
using Statistics
using StatsBase
using TOML

# Machine learning (for decoding/MVPA)
using LIBSVM

# TODO: consider using Threads.@threads for parallel processing?
# using Base.Threads

# Core types
include("types/core.jl")
include("types/pipeline.jl")
include("types/statistics.jl")
include("types/decoding.jl")
include("types/rsa.jl")

# Utility functions
include("utils/create_data.jl")
include("utils/error.jl")
include("utils/batch.jl")
include("utils/data.jl")
include("utils/files.jl")
include("utils/logging.jl")
include("utils/misc.jl")
include("utils/print.jl")
include("utils/extern/load_csv.jl")

# Layout handling
include("layouts/layout.jl")

# Data import/export
include("io/eeglab.jl")

# Configuration system
include("config/config.jl")

# Analysis functions
include("analysis/processing/apply.jl")
include("analysis/processing/artifact_detection.jl")
include("analysis/processing/baseline.jl")
include("analysis/processing/channel_average.jl")
include("analysis/processing/channel_delete.jl")
include("analysis/processing/channel_difference.jl")
include("analysis/processing/channel_repair.jl")
include("analysis/processing/channel_summary.jl")
include("analysis/processing/channel_metrics.jl")
include("analysis/processing/condition_combine.jl")
include("analysis/processing/condition_difference.jl")
include("analysis/processing/epochs.jl")
include("analysis/processing/filter.jl")
include("analysis/processing/ica.jl")
include("analysis/processing/mirror.jl")
include("analysis/processing/rereference.jl")
include("analysis/processing/resample.jl")
include("analysis/processing/triggers.jl")
include("analysis/processing/read_data.jl")

include("analysis/erps/jackknife_average.jl")
include("analysis/erps/erp_measurements.jl")
include("analysis/erps/gfp.jl")
include("analysis/erps/grand_average.jl")
include("analysis/erps/lrp.jl")
include("analysis/erps/realign.jl")

# Statistics submodules (loaded in dependency order)
include("analysis/statistics/core.jl")
include("analysis/statistics/thresholding.jl")
include("analysis/statistics/clustering.jl")
include("analysis/statistics/permutations.jl")
include("analysis/statistics/inference.jl")
include("analysis/statistics/statistics.jl")

# time-frequency analysis
include("analysis/time_frequency/tf_morlet.jl")
include("analysis/time_frequency/tf_stft.jl")
include("analysis/time_frequency/tf_multitaper.jl")
include("analysis/time_frequency/baseline.jl")
include("analysis/time_frequency/utils/utils.jl")

# decoding analysis via libsvm
include("analysis/decoding/decoding.jl")
include("analysis/decoding/preparation.jl")
include("analysis/decoding/decoding_statistics.jl")

# RSA (Representational Similarity Analysis)
include("analysis/rsa/rsa_models.jl")
include("analysis/rsa/rsa.jl")
include("analysis/rsa/rsa_noise_ceiling.jl")
include("analysis/rsa/rsa_crossvalidation.jl")

# analysis pipelines
include("pipelines/utils/utils.jl")
include("pipelines/utils/template.jl")
include("pipelines/pipeline_v1.jl")
include("pipelines/pipeline_v2.jl")

# Plotting functions
include("plots/utils/plot_misc.jl")

# Layout system (must be included before other plotting functions)
include("plots/utils/layout_system.jl")

# Shared interactivity (must be included before plot_erp and plot_epochs)
include("plots/utils/shared_interactivity.jl")

include("plots/processing/plot_artifacts.jl")
include("plots/processing/plot_channel_summary.jl")
include("plots/processing/plot_correlation_heatmap.jl")
include("plots/processing/plot_databrowser.jl")
include("plots/processing/plot_ica.jl")
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
include("plots/erps/plot_statistics.jl")
include("plots/erps/plot_topography.jl")

include("plots/time_frequency/plot_time_frequency.jl")
include("plots/time_frequency/plot_spectrum.jl")
include("plots/time_frequency/plot_power_spectrum.jl")

include("plots/decoding/plot_decoding.jl")
include("plots/rsa/plot_rsa.jl")

include("plots/plot_gui.jl")
include("plots/utils/help_system.jl")

# Demos
include("demos/signal_example_1.jl")
include("demos/signal_example_2.jl")
include("demos/simulate_erp.jl")

end # module EegFun 
