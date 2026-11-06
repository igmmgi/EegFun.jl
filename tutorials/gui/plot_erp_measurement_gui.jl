# Demo: Interactive ERP Measurement GUI
# Explore ERP amplitude and latency measurements interactively with live visual feedback
# before running batch extraction with erp_measurements().

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# gui = EegFun.plot_erp_measurement_gui("/path/to/participant_1_erps.jld2")

using EegFun
using GLMakie

#######################################################################
# 1. QUICK DEMO (DEFAULT BUNDLED ERPS, 30 HZ LOW-PASS, :Cz, 300-500 MS)
#######################################################################

# Controls to test:
# - Channel menu & multi-channel ROI averaging
# - Measurement Type menu (mean amplitude, max/min peak amplitude & latency,
#   peak-to-peak, rectified/positive/negative area, fractional area/peak latency)
# - Analysis Window & Baseline Interval dual-handle sliders
# - Local peak interval & fractional threshold sliders + marker toggles
gui = EegFun.plot_erp_measurement_gui()

#######################################################################
# 2. PRE-CONFIGURED LATENCY OR MULTI-CHANNEL ROI EXPLORATION
#######################################################################

erps = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
erps_lp = EegFun.lowpass_filter(erps, 30.0)

# Explore fractional area latency across a centro-parietal ROI:
# gui = EegFun.plot_erp_measurement_gui(
#     erps_lp,
#     channel = [:Cz, :CPz, :Pz],
#     analysis_type = "fractional_area_latency",
#     analysis_interval = (0.25, 0.55),
#     baseline_interval = (-0.2, 0.0),
# )

# Keep window open when run as a script from terminal
if !isinteractive() && !isnothing(gui.screen)
    wait(gui.screen)
end
