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
using GLMakie
using CairoMakie
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
include("plots/plot_grid_topo.jl")
include("plots/plot_power_spectrum.jl")

# ===============================
# CORE TYPES AND DATA STRUCTURES
# ===============================
export ContinuousData
export EpochData
export ErpData
export AnalysisInfo
export EpochCondition
export ValidationResult
export ConfigParameter
export IcaPrms
export InfoIca
export IntervalIdx
export IntervalTime
export CoordXY
export Coord

# ===============================
# DATA I/O AND CONFIGURATION
# ===============================
export read_bdf
export read_layout
export create_eeg_dataframe
export load_config
export show_parameter_info
export generate_config_template
export print_config

# ===============================
# DATA PROPERTIES AND ACCESSORS
# ===============================
export channels
export times
export sample_rate
export n_samples
export n_channels
export n_epochs
export duration
export has_channels
export common_channels
export extra_channels
export reference
export filter_info
export data
export n_average
export datarange
export data_limits_x
export data_limits_y

# ===============================
# PREPROCESSING FUNCTIONS
# ===============================
export preprocess_eeg_data
export parse_epoch_conditions
export filter_data
export filter_data!
export rereference
export rereference!
export baseline
export baseline!
export diff_channel
export diff_channel!
export detect_eog_onsets!
export is_extreme_value!

# ===============================
# EPOCHING AND ERP ANALYSIS
# ===============================
export extract_epochs
export remove_bad_epochs
export average_epochs
export search_sequence
export mark_epoch_windows!
export plot_epochs_table
export epochs_to_dataframe
export get_mean_amplitude
export get_peak_latency

# ===============================
# ICA ANALYSIS
# ===============================
export run_ica
export infomax_ica
export remove_ica_components
export restore_original_data
export create_ica_data_matrix

# ===============================
# STATISTICAL ANALYSIS
# ===============================
export find_bad_channels
export repair_epochs!
export AutoRejectParams
export n_extreme_value
export channel_joint_probability
export compute_probability!
export trim_extremes

# ===============================
# TIME-FREQUENCY ANALYSIS
# ===============================
export generate_signal
export average_over_trials
export apply_tf_baseline_db
export apply_tf_baseline_perchange
export plot_tf

# ===============================
# PLOTTING FUNCTIONS
# ===============================
# Main plotting functions
export plot_databrowser
export plot_epochs
export plot_erp
export plot_erp_image
export plot_topoplot
export plot_grid_topo
export plot_events
export plot_joint_probability
export plot_channel_summary

# Spectral plotting
export plot_channel_spectrum
export plot_component_spectrum
export plot_ica_component_spectrum
export plot_power_spectrum

# ICA plotting
export plot_ica_components
export plot_component_topo

# Layout plotting
export plot_layout_2d
export plot_layout_3d
export plot_layout_2d!
export plot_layout_3d!

# ===============================
# LAYOUT AND SPATIAL UTILITIES
# ===============================
export polar_to_cartesian_xy!
export polar_to_cartesian_xyz!
export get_electrode_neighbours_xy
export get_electrode_neighbours_xyz
export print_neighbours_dict
export average_number_of_neighbours
export distance_xy
export squared_distance_xy
export distance_xyz
export squared_distance_xyz

# ===============================
# ARRAY AND DATA UTILITIES
# ===============================
export channel_number_to_channel_label
export consecutive
export splitgroups
export find_idx_range
export find_idx_start_end
export detrend
export extract_int
export validate_baseline_interval
export colmeans

# ===============================
# GEOMETRIC UTILITIES
# ===============================
export create_convex_hull
export orientation
export best_rect

# ===============================
# FILE AND SYSTEM UTILITIES
# ===============================
export get_files
export check_files_exist
export make_output_filename
export basename_without_ext
export get_git_commit
export get_eegfun_version

# ===============================
# PRINTING AND DISPLAY UTILITIES
# ===============================
export print_vector
export print_vector_

end # module eegfun 