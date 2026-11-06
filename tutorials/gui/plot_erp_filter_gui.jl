# Demo: ERP Filter GUI
# Shows how to launch the interactive filter exploration GUI to compare
# IIR (Butterworth) and FIR filter settings on ERP and continuous waveforms.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# gui = EegFun.plot_erp_filter_gui("/path/to/participant_1_erps.jld2")

using EegFun
using GLMakie

#######################################################################
# 1. QUICK DEMO (BUNDLED ERP DATA)
#######################################################################

# Controls to test:
# - Channel dropdown, Filter Method (Butterworth IIR vs FIR Hamming)
# - Filter Function: filtfilt (2-pass zero-phase), filt_zerophase (1-pass FIR, MNE/EEGLAB style), filt (1-pass causal)
# - Mirror padding toggle, Low-pass & High-pass enable toggles, Cutoff & Order sliders (plus live FIR tap/delay readout)
gui = EegFun.plot_erp_filter_gui()

#######################################################################
# 2. MULTI-CONDITION ERP & SPECIFIC CHANNEL
#######################################################################

erps = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
# gui = EegFun.plot_erp_filter_gui(erps, channel = :Pz)

#######################################################################
# 3. CONTINUOUS DATA FILTER EXPLORATION
#######################################################################

# raw = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
# dat = EegFun.create_eegfun_data(raw)
# EegFun.rereference!(dat, :avg)
# gui = EegFun.plot_erp_filter_gui(dat, channel = :Cz)

# Keep window open when run as a script from terminal
if !isinteractive() && !isnothing(gui.screen)
    wait(gui.screen)
end
