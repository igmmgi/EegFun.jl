# Demo: Plotting ERPs with Measurement Overlays
# Computes and visualizes ERP measurements in a single call.

# Note: EegFun.example_path() resolves bundled example data paths.
# When using your own data, simply pass the file path directly, e.g.:
# dat = EegFun.read_raw_data("/path/to/your/data.bdf")

using EegFun

dat = EegFun.read_data(EegFun.example_path("data/julia/erps/example1_erps_good.jld2"))

# Mean amplitude
# EegFun.plot_erp_measurements(dat, "mean_amplitude", analysis_interval = (0.3, 0.5), baseline_interval = (-0.2, 0.0)) # not really useful as too crowded!

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
    EegFun.example_path("data/julia/erps/example1_erps_good.jld2"),
    "max_peak_amplitude",
    analysis_interval = (0.0, 1.0),
    baseline_interval = (-0.2, 0.0),
    condition_selection = EegFun.conditions([1]),
    channel_selection = EegFun.channels([:Cz, :Pz]),
    average_channels = true,
)
