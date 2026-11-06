# Demo: Plot GUI (Master Hub)
# Shows how to launch the central point-and-click Plot GUI for exploring
# single files or batch file patterns across all major EegFun plot types.

# Note: EegFun.example_path() resolves bundled example data paths.
# You can print these paths below to quickly paste or browse to them in the GUI.

using EegFun
using GLMakie

#######################################################################
# HELPFUL BUNDLED PATHS FOR TESTING IN THE GUI
#######################################################################

raw_bdf_path = EegFun.example_path("data/bdf/example1.bdf")
erp_jld2_path = EegFun.example_path("data/julia/erps/example1_erps.jld2")
layout_csv_path = EegFun.example_path("layouts/biosemi/biosemi72.csv")

@info "Example paths for testing in Plot GUI:" raw_bdf_path erp_jld2_path layout_csv_path

#######################################################################
# LAUNCH THE PLOT GUI
#######################################################################

# Controls to test:
# 1. Column 1 (File & Plot Selection):
#    - Select Directory / Select File (or paste `raw_bdf_path` / `erp_jld2_path`)
#    - Select Layout (`layout_csv_path` populates the Channel(s) dropdown)
#    - Plot Type: Data Browser, Epochs, ERP, Topography, ERP Image, Time-Frequency,
#      plus submenus for ICA >, ERP Analysis >, and Diagnostic Plots >
# 2. Column 2 (Data Filters):
#    - Participant, Condition, Epoch filters, Layout (Single / Single Avg / Grid / Topo), Channel(s)
# 3. Column 3 (Axis Settings):
#    - X / Y / Z limits, Baseline interval, TF baseline type, Invert Y, and PLOT button

gui = EegFun.plot_gui()
