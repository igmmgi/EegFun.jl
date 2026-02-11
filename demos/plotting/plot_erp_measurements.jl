# Demo: Plotting ERPs with Measurement Overlays
# Computes and visualizes ERP measurements in a single call.

using EegFun

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")

# Mean amplitude
EegFun.plot_erp_measurements(dat, "mean_amplitude", analysis_interval = (0.3, 0.5), baseline_interval = (-0.2, 0.0))

# Peak latency with specific channels
EegFun.plot_erp_measurements(
    dat,
    "max_peak_latency",
    analysis_interval = (0.0, 1.0),
    baseline_interval = (-0.2, 0.0),
    channel_selection = EegFun.channels([:Cz, :Pz]),
)

# Grid layout
EegFun.plot_erp_measurements(dat, "max_peak_amplitude", analysis_interval = (0.0, 1.0), baseline_interval = (-0.2, 0.0), layout = :grid)

# Load from file path
EegFun.plot_erp_measurements(
    "./resources/data/julia/erps/example1_erps_good.jld2",
    "max_peak_amplitude",
    analysis_interval = (0.0, 1.0),
    baseline_interval = (-0.2, 0.0),
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Cz, :Pz]),
)
