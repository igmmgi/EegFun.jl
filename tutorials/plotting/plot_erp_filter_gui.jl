# Demo: ERP Filter GUI
# Shows how to launch the interactive filter exploration tool
# to compare the effect of different filter settings on ERP waveforms.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun
using GLMakie

# Launch with default demo data (bundled example ERPs, channel :Cz):
gui = EegFun.plot_erp_filter_gui()

# Alternatively, specify a channel or custom data / file:
# gui = EegFun.plot_erp_filter_gui(channel = :Pz)
# dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps.jld2"))
# gui = EegFun.plot_erp_filter_gui(dat, channel = :Cz)

# Keep window open when run from terminal
if !isinteractive() && !isnothing(gui.screen)
    wait(gui.screen)
end
