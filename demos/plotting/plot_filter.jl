# Demo: Filter Visualization
# Shows filter frequency response and impulse response visualization.

using EegFun

# Create lowpass IIR filter using create_filter
filter_info = EegFun.create_lowpass_filter(30.0, 256.0; filter_method = "iir")

# Plot filter response
EegFun.plot_filter_response(filter_info)

# Test with custom parameters
EegFun.plot_filter_response(
    filter_info,
    title = "Custom Lowpass Filter Plot",
    color = :blue,
    linewidth = 3,
    reference_lines = [-3, -12, -24],
    reference_color = :red,
    n_points = 1000,
)

# Create highpass IIR filter using create_filter
filter_info = EegFun.create_highpass_filter(1.0, 256.0; filter_method = "iir")
EegFun.plot_filter_response(filter_info, title = "High-pass Filter", color = :green)

# Create FIR filter using create_filter
filter_info = EegFun.create_lowpass_filter(40.0, 256.0; filter_method = "fir")
EegFun.plot_filter_response(filter_info, title = "FIR Lowpass Filter", color = :purple)

# Test additional filter with separate plotting
filter_info = EegFun.create_highpass_filter(0.5, 256.0; filter_method = "iir")
EegFun.plot_filter_response(filter_info, title = "High-pass Filter with Plot")





