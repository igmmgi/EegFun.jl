# Packages
# Data handling
using BioSemiBDF
using CSV
using DataFrames
using JLD2
using TOML

# Signal processing
using DSP
using LinearAlgebra

# Plotting and visualization
using GLMakie
using LibGEOS
using ScatteredInterpolation

# Utilities
using Logging
using OrderedCollections
using Printf
using StatsBase
using Random

# eegfun.jl
# Core types
include("types/types.jl")

# Utility functions
include("utils/utils.jl")

# Layout handling
include("layouts/layout.jl")

# Analysis functions
include("config/config.jl")
include("analysis/tf.jl")
include("analysis/channel_difference.jl")
include("analysis/epochs.jl")
include("analysis/filter.jl")
# include("analysis/hanning.jl")
include("analysis/ica.jl")
include("analysis/rereference.jl")
include("analysis/analyse.jl")
include("analysis/autoreject.jl")
include("analysis/baseline.jl")

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
include("plots/plot_grid_topo.jl")

# # Miscellaneous functions
# include("misc/keep.jl")
# include("misc/test.jl")
# include("misc/test_scripts.jl")
# include("misc/test_tf.jl")
