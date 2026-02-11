# Demo: Interactive ERP Measurement GUI
# Explore ERP measurements interactively with live visual feedback.
# Useful for teaching, parameter exploration, and visual validation
# before batch processing with erp_measurements().

using EegFun

dat = EegFun.read_data("./resources/data/julia/erps/example1_erps_good.jld2")

# All conditions overlaid
EegFun.plot_erp_measurement_gui(dat)

# Single condition
EegFun.plot_erp_measurement_gui(dat[1])

# Select some specific initial settings
EegFun.plot_erp_measurement_gui(
    dat[1],
    channel = :Cz,
    analysis_type = "mean_amplitude",
    analysis_interval = (0.3, 0.5),
    baseline_interval = (-0.2, 0.0),
)
