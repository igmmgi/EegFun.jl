# package
using EegFun
using GLMakie
# using CairoMakie
# using BenchmarkTools

function get_data()

    # read **.bdf file 
    dat = EegFun.read_bdf("../Flank_C_3.bdf")
    layout = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv")
    dat = EegFun.create_eeg_dataframe(dat, layout)

    # neighbours
    EegFun.polar_to_cartesian_xy!(layout)
    EegFun.get_neighbours_xy!(layout, 40)
    EegFun.polar_to_cartesian_xyz!(layout)
    EegFun.get_neighbours_xyz!(layout, 40)

    # basic preprocessing
    EegFun.highpass_filter!(dat, 1)
    EegFun.rereference!(dat, :avg)

    # epoching
    epoch_cfg = [
        EegFun.EpochCondition(name = "ExampleEpoch1", trigger_sequences = [[1]]),
        EegFun.EpochCondition(name = "ExampleEpoch2", trigger_sequences = [[3]]),
    ]
    epochs = EegFun.EpochData[]
    for (idx, epoch) in enumerate(epoch_cfg)
        push!(epochs, EegFun.extract_epochs(dat, idx, epoch, -2, 4))
    end

    # ERP average
    erps = [EegFun.average_epochs(epoch) for epoch in epochs]

    return dat, epochs, erps, layout
end


###########################
# Channel Summary Tests
###########################
function test_channel_summary()

    # setup data and call function
    println("\n=== Testing Channel Summary Plots ===")
    dat, epochs, erps, layout = get_data()

    cs = EegFun.channel_summary(dat)

    # Plots
    fig, ax = EegFun.plot_channel_summary(cs, :max)
    fig, ax = EegFun.plot_channel_summary(cs, :max, ylabel = "Max Value", title = "Channel Max", bar_color = :blue)
    fig, ax = EegFun.plot_channel_summary(cs, :std, title = "Channel Standard Deviation")
    fig, ax = EegFun.plot_channel_summary(cs, :range, title = "Channel Range")

    println("✓ Channel summary plots completed")
end
test_channel_summary()

###########################
# Correlation Heatmap Tests
###########################
function test_correlation_heatmap()

    # setup data and call function
    println("\n=== Testing Correlation Heatmap Plots ===")
    dat, epochs, erps, layout = get_data()
    cm = EegFun.correlation_matrix(dat)

    # Plots
    fig, ax = EegFun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")

    cm = EegFun.correlation_matrix(dat, channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :F4, :C3, :C4]))
    fig, ax = EegFun.plot_correlation_heatmap(cm, title = "Frontal-Central Correlations", colorrange = (0, 1), colormap = :viridis)
    fig, ax = EegFun.plot_correlation_heatmap(cm, title = "Masked Correlations (0.3-0.7)", mask_range = (0.3, 0.7), colorrange = (-1, 1))

    println("✓ Correlation heatmap plots completed")
end
test_correlation_heatmap()

###########################
# Joint Probability Tests
###########################
function test_joint_probability()

    # setup data and call function
    println("\n=== Testing Joint Probability Plots ===")
    dat, epochs, erps, layout = get_data()

    jp = EegFun.channel_joint_probability(dat)

    # Plots 
    fig, ax = EegFun.plot_joint_probability(jp, title = "Channel Joint Probability")
    fig, ax = EegFun.plot_joint_probability(jp, title = "Joint Probability - Custom Style", bar_color = :red, sort_values = true)

    println("✓ Joint probability plots completed")
end
test_joint_probability()

###########################
# Layout Plot Tests
###########################
function test_layout_plots()
    println("\n=== Testing Layout Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_layout_2d(layout, title = "2D Electrode Layout")
    fig, ax = EegFun.plot_layout_3d(layout, title = "3D Electrode Layout")
    fig, ax = EegFun.plot_layout_2d(
        layout,
        title = "Custom 2D Layout",
        head_color = :red,
        point_color = :blue,
        point_markersize = 15,
        label_fontsize = 12,
    )

    println("✓ Layout plots completed")
end
test_layout_plots()

###########################
# Topography Tests
###########################
function test_topography_plots()
    println("\n=== Testing Topography Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_topography(dat)
    fig, ax = EegFun.plot_topography(dat, title = "Custom Topography", colormap = :jet, gridscale = 100, method = :spherical_spline)
    fig, ax = EegFun.plot_topography(dat, title = "Test Title", title_fontsize = 20, show_title = true)
    fig, ax = EegFun.plot_topography(dat, show_title = false)
    fig, ax = EegFun.plot_topography(
        dat,
        title = "Test with New Parameters",
        colorrange = (-2, 2),
        nan_color = :red,
        grid_visible = true,
        xlabel = "X Position",
        ylabel = "Y Position",
        label_fontsize = 16,
        tick_fontsize = 14,
        colorbar_fontsize = 14,
    )

    fig, ax = EegFun.plot_topography(dat, title = "Non-interactive Plot", interactive = false, display_plot = false)

    println("✓ Topography plots completed")
end
test_topography_plots()

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



###########################
# ERP Plot Tests
###########################
function test_erp_plots()
    println("\n=== Testing ERP Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_erp(erps[1], title = "ERP Waveform")
    fig, ax = EegFun.plot_erp(erps, title = "Multiple ERP Conditions")
    fig, ax = EegFun.plot_erp(erp, title = "Custom ERP Plot", color = :red, linewidth = 2, ylabel = "Amplitude (μV)")

    println("✓ ERP plots completed")
end
test_erp_plots()

###########################
# Epochs Plot Tests
###########################
function test_epoch_plots()
    println("\n=== Testing Epochs Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_epochs(epochs[1], title = "Epochs Plot")
    fig, ax = EegFun.plot_epochs(epochs[1], title = "Custom Epochs Plot", color = [:blue, :red], alpha = [0.7, 1.0])

    println("✓ Epochs plots completed")
end
test_epoch_plots()

###########################
# ERP Image Tests
###########################
function test_erp_image()
    println("\n=== Testing ERP Image Plots ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    fig, ax = EegFun.plot_erp_image(epochs[1], title = "ERP Image")
    fig, ax = EegFun.plot_erp_image(epochs[1], title = "ERP Image with Smoothing", boxcar_average = 5)
    fig, ax = EegFun.plot_erp_image(epochs[1], title = "Custom ERP Image", colormap = :viridis)

    println("✓ ERP image plots completed")
end
test_erp_image()

###########################
# Data Browser Tests
###########################
function test_databrowser()
    println("\n=== Testing Data Browser ===")
    dat, epochs, erps, layout = get_data()

    # Plots
    # fig, ax = EegFun.plot_databrowser(dat) 
    fig, ax = EegFun.plot_databrowser(epochs[1])

    println("✓ Data browser plots completed")
end
test_databrowser()
