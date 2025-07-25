__precompile__(true)  

module eegfun

# Core dependencies
using BioSemiBDF
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
using Logging
using LoggingExtras
using OrderedCollections
using Printf
using StatsBase
using Random
using Dates
using TOML
using PrettyTables

# Julia standard library
using Base.Threads

# Core types
include("types/types.jl")

# Utility functions
include("utils/data.jl")
include("utils/error.jl")
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
include("analysis/autoreject.jl")
include("analysis/baseline.jl")
include("analysis/channel_difference.jl")
include("analysis/epochs.jl")
include("analysis/filter.jl")
include("analysis/ica.jl")
include("analysis/preprocess.jl")
include("analysis/rereference.jl")
include("analysis/tf.jl")

# Plotting functions
include("plots/constants.jl")
include("plots/plot_misc.jl")

include("plots/plot_channel_summary.jl")
include("plots/plot_correlation_heatmap.jl")
include("plots/plot_databrowser.jl")
include("plots/plot_epochs.jl")
include("plots/plot_erp.jl")
include("plots/plot_erp_grid.jl")
include("plots/plot_erp_image.jl")
include("plots/plot_grid_topo.jl")
include("plots/plot_ica.jl")
include("plots/plot_joint_probability.jl")
include("plots/plot_layout.jl")
include("plots/plot_power_spectrum.jl")
include("plots/plot_topography.jl")
include("plots/plot_trigger_overview.jl")
include("plots/plot_trigger_timing.jl")

end # module eegfun 