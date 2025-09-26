__precompile__(true)

module eegfun

# Core dependencies
using BiosemiDataFormat
using BrainVisionDataFormat
using CSV
using DataFrames
using JLD2

# Signal processing
using DSP
using LinearAlgebra
using FFTW

# Plotting and visualization
using Makie
using ScatteredInterpolation

# Utilities
using DocStringExtensions
using Logging
using LoggingExtras
using OrderedCollections
using Printf
using Statistics
using StatsBase
using Random
using Dates
using TOML
using PrettyTables

# Julia standard library
# using Base.Threads

# Core types
include("types/types.jl")

# Utility functions
include("utils/error.jl")
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
include("analysis/analyse.jl")
include("analysis/artifact_detection.jl")
include("analysis/autoreject.jl")
include("analysis/baseline.jl")
include("analysis/channel_average.jl")
include("analysis/channel_difference.jl")
include("analysis/channel_summary.jl")
include("analysis/channel_metrics.jl")
include("analysis/epochs.jl")
include("analysis/filter.jl")
include("analysis/ica.jl")
include("analysis/preprocess.jl")
include("analysis/rereference.jl")
include("analysis/tf.jl")
include("analysis/triggers.jl")

# Plotting functions
include("plots/utils/plot_misc.jl")

# Layout system (must be included before other plotting functions)
include("plots/utils/layout_system.jl")

# Shared interactivity (must be included before plot_erp and plot_epochs)
include("plots/utils/shared_interactivity.jl")

include("plots/plot_channel_summary.jl")
include("plots/plot_correlation_heatmap.jl")
include("plots/plot_databrowser.jl")
include("plots/plot_epochs.jl")
include("plots/plot_erp.jl")
include("plots/plot_erp_image.jl")
include("plots/plot_ica.jl")
include("plots/plot_joint_probability.jl")
include("plots/plot_layout.jl")
include("plots/plot_power_spectrum.jl")
include("plots/plot_topography.jl")
include("plots/plot_triggers.jl")
include("plots/plot_filter.jl")

# Examples
include("examples/signal_example_1.jl")
include("examples/signal_example_2.jl")
include("examples/plot_my_data_gui.jl")

end # module eegfun 
