# Demo: Interactive Viewers & Scrubbers
# Short walkthrough for testing widget-driven modes across EegFun's plotting functions:
# 2D & 3D topography time scrubbers, interactive power spectrum menus,
# trigger timing sliders, and interactive ERP / Epoch control panels.

using EegFun
using GLMakie

#######################################################################
# LOAD BUNDLED EXAMPLE DATA
#######################################################################

raw = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(raw, layout)
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 0.5)

erps = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
erps_lp = EegFun.lowpass_filter(erps, 30.0)

#######################################################################
# 1. 2D & 3D TOPOGRAPHY TIME SCRUBBER (scrubber = true)
#######################################################################

# Controls to test:
# - Time slider, "<" / "Play" ("Pause") / ">" buttons, and "Mean" button
# - In 3D mode: click-drag to rotate the 3D scalp mesh while scrubbing/playing!
EegFun.plot_topography(erps_lp[1])
EegFun.plot_topography_3d(erps_lp[1]; scrubber = true)

#######################################################################
# 2. INTERACTIVE CHANNEL POWER SPECTRUM (plot_channel_spectrum)
#######################################################################

# Controls to test:
# - X-Axis dropdown: Linear (Hz) / Log10 (Hz)
# - Y-Axis dropdown: Linear (μV²/Hz) / Log10 (μV²/Hz) / Decibels (dB)
EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fz, :Cz, :Pz, :Oz]))

#######################################################################
# 3. INTERACTIVE TRIGGER TIMING VIEWER (plot_trigger_timing)
#######################################################################

# Controls to test:
# - Bottom sliders for scrolling time position and adjusting visible window width
EegFun.plot_trigger_timing(dat; ignore_triggers = [3])

#######################################################################
# 4. INTERACTIVE ERP & EPOCHS CONTROL PANELS (plot_erp / plot_epochs)
#######################################################################

# Keyboard & mouse interactions to test in the figure window:
# - Press "c": Open the live Control Panel (Baseline Start/End + Apply, Condition checkboxes, Toggle Legend)
# - Press "l": Toggle legend visibility
# - Press "s": Open the Save Figure Settings dialog (PNG/PDF/SVG/JPG, DPI, cm dimensions)
# - Press "i": Show shortcut help popup
# - Drag-select a time interval + right-click inside it: Launch 2D/3D topography popup!
EegFun.plot_erp(erps_lp, channel_selection = EegFun.channels(:Cz))

epoch_cfg =
    [EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]), EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]])]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-0.2, 1.0))
EegFun.plot_epochs(epochs[1], channel_selection = EegFun.channels(:Cz))
