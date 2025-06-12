__precompile__(true)  # Temporarily disable precompilation due to Julia environment issues

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
# using GLMakie
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
include("utils/misc.jl")
include("utils/data.jl")
include("utils/error.jl")
include("utils/files.jl")
include("utils/print.jl")
include("utils/logging.jl")

# Layout handling
include("layouts/layout.jl")

# Configuration system
include("config/config.jl")

# Analysis functions
include("analysis/tf.jl")
include("analysis/channel_difference.jl")
include("analysis/epochs.jl")
include("analysis/filter.jl")
include("analysis/ica.jl")
include("analysis/rereference.jl")
include("analysis/analyse.jl")
include("analysis/autoreject.jl")
include("analysis/baseline.jl")
include("analysis/preprocess.jl")

# Plotting functions
include("plots/plot_ica.jl")
include("plots/plot_layout.jl")
include("plots/plot_misc.jl")
include("plots/plot_topo.jl")
include("plots/plot_databrowser.jl")
include("plots/plot_epochs.jl")
include("plots/plot_erp.jl")
include("plots/plot_erp_grid.jl")
include("plots/plot_erp_image.jl")
include("plots/plot_events.jl")
include("plots/plot_events_timing.jl")
include("plots/plot_grid_topo.jl")
include("plots/plot_power_spectrum.jl")

end # module eegfun 