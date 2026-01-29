# package
using EegFun
using GLMakie

###########################
# Filter Tests
###########################
function test_plot_filter()
    println("\n=== Testing Filter Plots ===")

    # Create lowpass IIR filter using create_filter
    filter_info = EegFun.create_filter("lp", "iir", 30.0, 256.0)

    # Plot filter response
    fig, ax = EegFun.plot_filter_response(filter_info)

    # Test with custom parameters
    fig, ax = EegFun.plot_filter_response(
        filter_info,
        title = "Custom Lowpass Filter Plot",
        actual_color = :blue,
        actual_linewidth = 3,
        reference_lines = [-3, -12, -24],
        reference_color = :red,
        n_points = 1000,
        display_plot = false,
    )

    # Create highpass IIR filter using create_filter
    filter_info2 = EegFun.create_filter("hp", "iir", 1.0, 256.0)

    fig, ax = EegFun.plot_filter_response(filter_info2, title = "High-pass Filter", actual_color = :green, display_plot = false)

    # Create FIR filter using create_filter
    filter_info3 = EegFun.create_filter("lp", "fir", 40.0, 256.0)

    fig, ax = EegFun.plot_filter_response(filter_info3, title = "FIR Lowpass Filter", actual_color = :purple, display_plot = false)

    # Test additional filter with separate plotting
    filter_info4 = EegFun.create_filter("hp", "iir", 0.5, 256.0)
    fig, ax = EegFun.plot_filter_response(filter_info4, title = "High-pass Filter with Plot", actual_color = :orange, display_plot = false)

    println("✓ Filter plots completed")
end
test_plot_filter()

###########################
# Power Spectrum Tests
###########################
function test_plot_power_spectrum()
    println("\n=== Testing Power Spectrum Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_channel_spectrum(dat, channel_selection = EegFun.channels([:Fp1]), title = "Fp1 Power Spectrum")
    fig, ax = EegFun.plot_channel_spectrum(
        dat,
        channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :F4]),
        title = "Frontal Channels Power Spectrum",
    )
    fig, ax = EegFun.plot_channel_spectrum(
        dat,
        channel_selection = EegFun.channels([:Fp1, :Fp2]),
        title = "Custom Power Spectrum",
        x_scale = :log10,
        y_scale = :log10,
        max_freq = 100.0,
        window_size = 512,
        line_width = 3,
        show_freq_bands = false,
    )
    fig, ax = EegFun.plot_channel_spectrum(
        dat,
        channel_selection = EegFun.channels([:Fp1]),
        title = "Power Spectrum with Hamming Window",
        window_function = EegFun.DSP.hamming,
        overlap = 0.75,
    )

    println("✓ Power spectrum plots completed")
end
test_plot_power_spectrum()


