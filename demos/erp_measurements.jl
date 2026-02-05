"""
Tutorial: ERP Measurement Options

This script provides an introduction to the ERP measurement capabilities 
in EegFun for extracting quantitative features from ERP data.

1. Amplitude measurements (mean, peak)
2. Latency measurements (peak, fractional)
3. Area/integral measurements
4. Peak-to-peak measurements
"""

using EegFun
dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")

# We can use the plot_erp_measurements_gui to explore the data and select the measurement parameters
EegFun.plot_erp_measurement_gui(dat)    # all conditions
EegFun.plot_erp_measurement_gui(dat[1]) # first condition

# ----------------------------------------------------------------------------
# Amplitude Measurements
# ----------------------------------------------------------------------------

# batch type analyses
input_dir = "./resources/data/julia/erps"
file_pattern = "erps_good"

# Mean amplitude in a time window
mean_amp = EegFun.erp_measurements(
    file_pattern,
    "max_peak_latency",
    input_dir = input_dir,
    condition_selection = EegFun.conditions([1, 2]),
    channel_selection = EegFun.channels(),  # all channels
    # channel_selection = EegFun.channels([:Pz, :Cz, :Fz]),
    analysis_interval = (0.6, 0.8),
    baseline_interval = (-0.2, 0.0),  # 200 ms pre-stimulus baseline
)

# the above results data AND saves the results to a csv file
