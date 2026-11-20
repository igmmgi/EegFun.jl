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

# Core types, utilities, layouts, I/O, and configuration
include("types.jl")
include("utils.jl")
include("layouts/layout.jl")
include("io.jl")
include("config/config.jl")

# Analysis subsystems
include("analysis/processing.jl")
include("analysis/erps.jl")
include("analysis/ride.jl")
include("analysis/statistics.jl")
include("analysis/time_frequency.jl")
include("analysis/decoding.jl")
include("analysis/rsa.jl")

# Automated pipelines, plotting, and interactive demos
include("pipelines.jl")
include("plots.jl")
include("demos.jl")

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
