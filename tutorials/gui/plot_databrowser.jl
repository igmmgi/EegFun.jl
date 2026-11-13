# Demo: Interactive Databrowser
# Shows the interactive databrowser for continuous and epoched data with ICA,
# file-picker launcher mode, analysis settings, custom channels, and JLD2 loading.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie
using JLD2

const DEMO_OUTPUT = "./tutorials/output/"
mkpath(DEMO_OUTPUT)

#######################################################################
# QUICK LAUNCHERS (NO DATA SETUP REQUIRED)
#######################################################################

# Quick demo with bundled raw dataset and full 72-channel layout (press "i" for shortcut help):
EegFun.plot_databrowser()

# Open with interactive file-picker launcher dialog:
# EegFun.plot_databrowser(launcher = true)

#######################################################################
# CONTINUOUS DATA (WITHOUT & WITH LAYOUT)
#######################################################################

# Without layout (waveform browsing works; topographic popups disabled)
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
dat = EegFun.create_eegfun_data(dat)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)
EegFun.plot_databrowser(dat)

# With 72-channel BioSemi layout (enables 2D/3D topographic popups & channel repair)
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)
EegFun.highpass_filter!(dat, 0.5)
EegFun.plot_databrowser(dat)

#######################################################################
# INTERACTIVE ANALYSIS SETTINGS ROUNDTRIP
#######################################################################

# Toggle filters/rereference and drag-select regions in the GUI; settings are returned on close:
_, _, analysis_settings = EegFun.plot_databrowser(dat)

# Apply captured interactive settings programmatically:
dat_new = EegFun.apply_analysis_settings(dat, analysis_settings)
EegFun.plot_databrowser(dat_new)

#######################################################################
# EXTRA CHANNELS & BOOLEAN ARTIFACT OVERLAYS
#######################################################################

dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
dat = EegFun.create_eegfun_data(dat, layout)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

# Flag extreme values, EOG channels/onsets, and epoch intervals (Boolean columns) + add a custom continuous trace
EegFun.calculate_eog_channels!(dat, EegFun.EogConfig())
EegFun.detect_eog_onsets!(dat, EegFun.EogConfig())
EegFun.is_extreme_value!(dat, 100)
EegFun.mark_epoch_intervals!(dat, [1, 2], [-0.2, 1.0])
dat.data[:, :sinewave] = 50.0 .* sin.(2π .* 10.0 .* dat.data[:, :time])

# Inspect under the "Extra Channels" menu in the databrowser
EegFun.plot_databrowser(dat)

# Custom line & selection styling
EegFun.plot_databrowser(dat; channel_line_width = 2, selection_color = (:green, 0.1))

#######################################################################
# ICA OVERLAY & LIVE COMPONENT REMOVAL PREVIEW
#######################################################################

ica_result = EegFun.run_ica(dat; sample_selection = EegFun.samples_not(:is_extreme_value_100), percentage_of_data = 25)
EegFun.plot_databrowser(dat, ica_result)

#######################################################################
# EPOCHED DATA & DIRECT JLD2 FILE LOADING
#######################################################################

epoch_cfg = [EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.5, 1.0))

EegFun.plot_databrowser(epochs[1])
EegFun.plot_databrowser(epochs[1], ica_result)

# Save and browse directly from JLD2 files on disk
jldsave(joinpath(DEMO_OUTPUT, "dat.jld2"), data = dat)
jldsave(joinpath(DEMO_OUTPUT, "epochs.jld2"), data = epochs)
jldsave(joinpath(DEMO_OUTPUT, "ica.jld2"), data = ica_result)

EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "dat.jld2"))
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "epochs.jld2"))
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "dat.jld2"), joinpath(DEMO_OUTPUT, "ica.jld2"))
EegFun.plot_databrowser(joinpath(DEMO_OUTPUT, "epochs.jld2"), joinpath(DEMO_OUTPUT, "ica.jld2"))
