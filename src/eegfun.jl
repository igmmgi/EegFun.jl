__precompile__(true)


module eegfun

# Core dependencies
using BiosemiDataFormat
using BrainVisionDataFormat
using CategoricalArrays
using CSV
using DataFrames
using JLD2

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
using MLJ
using MLJLinearModels
using LIBSVM

# Julia standard library
# using Base.Threads

# Core types
include("types/core.jl")
include("types/pipeline.jl")
include("types/statistics.jl")
include("types/decoding.jl")

# Utility functions
include("utils/error.jl")
include("utils/batch.jl")
include("utils/data.jl")
include("utils/files.jl")
include("utils/logging.jl")
include("utils/misc.jl")
include("utils/print.jl")

# Layout handling
include("layouts/layout.jl")

# Configuration system
include("config/config.jl")

# Analysis functions
include("analysis/apply.jl")
include("analysis/artifact_detection.jl")
include("analysis/baseline.jl")
include("analysis/channel_average.jl")
include("analysis/channel_difference.jl")
include("analysis/channel_repair.jl")
include("analysis/channel_summary.jl")
include("analysis/channel_metrics.jl")
include("analysis/condition_combine.jl")
include("analysis/condition_difference.jl")
include("analysis/decoding.jl")
include("analysis/erp_measurements.jl")
include("analysis/epochs.jl")
include("analysis/filter.jl")
include("analysis/statistics_ttest.jl")
include("analysis/statistics.jl")
include("analysis/gfp.jl")
include("analysis/grand_average.jl")
include("analysis/ica.jl")
include("analysis/jackknife_average.jl")
include("analysis/lrp.jl")
include("analysis/mirror.jl")
include("analysis/realign.jl")
include("analysis/rereference.jl")
include("analysis/resample.jl")
include("analysis/time_frequency.jl")
include("analysis/triggers.jl")

# analysis pipelines
include("pipelines/utils/utils.jl")
include("pipelines/utils/template.jl")
include("pipelines/pipeline.jl")

# Plotting functions
include("plots/utils/plot_misc.jl")

# Layout system (must be included before other plotting functions)
include("plots/utils/layout_system.jl")

# Shared interactivity (must be included before plot_erp and plot_epochs)
include("plots/utils/shared_interactivity.jl")

include("plots/plot_artifacts.jl")
include("plots/plot_statistics.jl")
include("plots/plot_channel_summary.jl")
include("plots/plot_correlation_heatmap.jl")
include("plots/plot_databrowser.jl")
include("plots/plot_decoding.jl")
include("plots/plot_epochs_rejection.jl")
include("plots/plot_epochs.jl")
include("plots/plot_erp.jl")
include("plots/plot_erp_measurements.jl")
include("plots/plot_erp_image.jl")
include("plots/plot_gfp.jl")
include("plots/plot_ica.jl")
include("plots/plot_joint_probability.jl")
include("plots/plot_layout.jl")
include("plots/plot_power_spectrum.jl")
include("plots/plot_time_frequency.jl")
include("plots/plot_topography.jl")
include("plots/plot_triggers.jl")
include("plots/plot_filter.jl")

# Plot utilities
include("plots/utils/help_system.jl")

# Examples
include("examples/signal_example_1.jl")
include("examples/signal_example_2.jl")
include("examples/plot_gui.jl")
include("examples/simulate_erp.jl")

end # module eegfun 
