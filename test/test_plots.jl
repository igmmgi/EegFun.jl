# package
using eegfun
using GLMakie
# using CairoMakie
# using BenchmarkTools

function get_data() 
    dat = eegfun.read_bdf("../Flank_C_3.bdf");
    layout = eegfun.read_layout("./data/layouts/biosemi72.csv");
    dat = eegfun.create_eeg_dataframe(dat, layout);
    return dat, layout
end

function get_epochs_data()
    dat, layout = get_data()
    # Create some test epochs using extract_epochs
    condition = eegfun.EpochCondition( name = "test", trigger_sequences = [[1]])
    epochs = eegfun.extract_epochs(dat, 1, condition, -0.1, 0.5)
    return epochs, layout
end

function get_erp_data()
    epochs, layout = get_epochs_data()
    erp = eegfun.average_epochs(epochs)
    return erp, layout
end


###########################
# Channel Summary Tests
###########################
function test_channel_summary()

    # setup data and call function
    println("\n=== Testing Channel Summary Plots ===")
    dat, layout = get_data()
    cs = eegfun.channel_summary(dat) 

    # Basic plot
    fig, ax = eegfun.plot_channel_summary(cs, :max)
    
    # With custom styling
    fig, ax = eegfun.plot_channel_summary(cs, :max, ylabel = "Max Value", title = "Channel Max", bar_color = :blue)
    
    # Test different metrics
    fig, ax = eegfun.plot_channel_summary(cs, :std, title = "Channel Standard Deviation")
    fig, ax = eegfun.plot_channel_summary(cs, :range, title = "Channel Range")
    
    println("✓ Channel summary plots completed")
end
test_channel_summary()

###########################
# Correlation Heatmap Tests
###########################
function test_correlation_heatmap()

    # setup data and call function
    println("\n=== Testing Correlation Heatmap Plots ===")
    dat, layout = get_data()
    cm = eegfun.correlation_matrix(dat)

    # Basic plot
    fig, ax = eegfun.plot_correlation_heatmap(cm, title = "Full Correlation Matrix")
    
    # Subset of channels
    cm_subset = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Fp1, :Fp2, :F3, :F4, :C3, :C4]))
    fig, ax = eegfun.plot_correlation_heatmap(cm_subset, 
        title = "Frontal-Central Correlations",
        colorrange = (0, 1),
        colormap = :viridis
    )
    
    # With masking
    fig, ax = eegfun.plot_correlation_heatmap(cm_subset, 
        title = "Masked Correlations (0.3-0.7)",
        mask_range = (0.3, 0.7),
        colorrange = (-1, 1)
    )
    
    println("✓ Correlation heatmap plots completed")
end
test_correlation_heatmap()

###########################
# Joint Probability Tests
###########################
function test_joint_probability()
    
    # setup data and call function
    println("\n=== Testing Joint Probability Plots ===")
    dat, layout = get_data()
    jp = eegfun.channel_joint_probability(dat)
    
    # Basic plot
    fig, ax = eegfun.plot_joint_probability(jp, title = "Channel Joint Probability")
    
    # With custom styling
    fig, ax = eegfun.plot_joint_probability(jp, title = "Joint Probability - Custom Style", bar_color = :red, sort_values = true)
    
    println("✓ Joint probability plots completed")
end
test_joint_probability()

###########################
# Layout Plot Tests
###########################
function test_layout_plots()
    println("\n=== Testing Layout Plots ===")
    dat, layout = get_data()
    
    # 2D layout
    fig, ax = eegfun.plot_layout_2d(layout, title = "2D Electrode Layout")
    
    # 3D layout
    fig, ax = eegfun.plot_layout_3d(layout, title = "3D Electrode Layout")
    
    # With custom styling
    fig, ax = eegfun.plot_layout_2d(layout, 
        title = "Custom 2D Layout",
        head_color = :red,
        point_color = :blue,
        point_markersize = 15,
        label_fontsize = 12
    )
    
    println("✓ Layout plots completed")
end
test_layout_plots()

###########################
# Topography Tests
###########################
function test_topography_plots()
    println("\n=== Testing Topography Plots ===")
    dat, layout = get_data()
    
    # Basic topography
    fig, ax = eegfun.plot_topography(dat)
    
    # With custom styling
    fig, ax = eegfun.plot_topography(dat,  
        title = "Custom Topography",
        colormap = :jet,
        gridscale = 100,
        method = :spherical_spline
    )
    
    # Test title parameters
    fig, ax = eegfun.plot_topography(dat, 
        title = "Test Title",
        title_fontsize = 20,
        show_title = true
    )
    
    # Test with title disabled
    fig, ax = eegfun.plot_topography(dat, 
        show_title = false
    )
    
    # Test new parameters
    fig, ax = eegfun.plot_topography(dat, 
        title = "Test with New Parameters",
        colorrange = (-2, 2),
        nan_color = :red,
        grid_visible = true,
        xlabel = "X Position",
        ylabel = "Y Position",
        label_fontsize = 16,
        tick_fontsize = 14,
        colorbar_fontsize = 14
    )
    
    # Test with interactive disabled
    fig, ax = eegfun.plot_topography(dat, 
        title = "Non-interactive Plot",
        interactive = false,
        display_plot = false
    )
    
    println("✓ Topography plots completed")
end
test_topography_plots()

###########################
# Filter Tests
###########################
function test_plot_filter()
    println("\n=== Testing Filter Plots ===")
    dat, layout = get_data()
    
    # Create FilterInfo struct for lowpass IIR filter
    filter_type = "lp"
    filter_method = "iir"
    cutoff_freq = 30.0
    sample_rate = dat.sample_rate
    order = 2
    transition_width = 0.1
    
    # Create filter prototype/object based on type
    filter_prototypes = Dict("hp" => eegfun.Highpass, "lp" => eegfun.Lowpass)
    filter_prototype = filter_prototypes[filter_type](cutoff_freq)
    transition_band = cutoff_freq * transition_width
    n_taps = nothing
    
    # Create filter with chosen method
    if filter_method == "iir"
        filter_object = eegfun.digitalfilter(filter_prototype, eegfun.Butterworth(order); fs = sample_rate)
    end
    
    filter_info = eegfun.FilterInfo(
        filter_type,
        filter_object,
        filter_method,
        Float64(cutoff_freq),
        Float64(sample_rate),
        order,
        n_taps,
        Float64(transition_band),
    )
    
    # Plot filter response
    fig, ax = eegfun.plot_filter_response(filter_info)
    
    # Test with custom parameters
    fig, ax = eegfun.plot_filter_response(filter_info, 
        title = "Custom Lowpass Filter Plot",
        actual_color = :blue,
        actual_linewidth = 3,
        reference_lines = [-3, -12, -24],
        reference_color = :red,
        n_points = 1000,
        display_plot = false
    )
    
    # Create FilterInfo struct for highpass IIR filter
    filter_type2 = "hp"
    cutoff_freq2 = 1.0
    transition_width2 = 0.25
    
    filter_prototype2 = filter_prototypes[filter_type2](cutoff_freq2)
    transition_band2 = cutoff_freq2 * transition_width2
    
    filter_object2 = eegfun.digitalfilter(filter_prototype2, eegfun.Butterworth(order); fs = sample_rate)
    
    filter_info2 = eegfun.FilterInfo(
        filter_type2,
        filter_object2,
        filter_method,
        Float64(cutoff_freq2),
        Float64(sample_rate),
        order,
        n_taps,
        Float64(transition_band2),
    )
    
    fig, ax = eegfun.plot_filter_response(filter_info2, 
        title = "High-pass Filter",
        actual_color = :green,
        display_plot = false
    )
    
    # Create FilterInfo struct for FIR filter
    filter_type3 = "lp"
    filter_method3 = "fir"
    cutoff_freq3 = 40.0
    transition_width3 = 0.1
    
    filter_prototype3 = filter_prototypes[filter_type3](cutoff_freq3)
    transition_band3 = cutoff_freq3 * transition_width3
    
    # Calculate number of taps (ensure not too small + next pow2 + odd) for FIR filter
    n_taps3 = Int(ceil(4.0 * sample_rate / transition_band3))
    n_taps3 = max(n_taps3, 101)
    n_taps3 = nextpow(2, n_taps3)
    n_taps3 += 1  # Always add 1 to make odd as nextpow2 is even
    
    filter_object3 = eegfun.digitalfilter(filter_prototype3, eegfun.FIRWindow(eegfun.hamming(n_taps3)); fs = sample_rate)
    
    filter_info3 = eegfun.FilterInfo(
        filter_type3,
        filter_object3,
        filter_method3,
        Float64(cutoff_freq3),
        Float64(sample_rate),
        0,  # order not applicable for FIR
        n_taps3,
        Float64(transition_band3),
    )
    
    fig, ax = eegfun.plot_filter_response(filter_info3, 
        title = "FIR Lowpass Filter",
        actual_color = :purple,
        display_plot = false
    )
    
    println("✓ Filter plots completed")
end
test_plot_filter()

###########################
# Power Spectrum Tests
###########################
function test_plot_power_spectrum()
    println("\n=== Testing Power Spectrum Plots ===")
    dat, layout = get_data()

    # Single channel spectrum
    fig, ax = eegfun.plot_channel_spectrum(dat, 
        channel_selection = eegfun.channels([:Fp1]),
        title = "Fp1 Power Spectrum",
        display_plot = false
    )
    
    # Multiple channels
    fig, ax = eegfun.plot_channel_spectrum(dat, 
        channel_selection = eegfun.channels([:Fp1, :Fp2, :F3, :F4]),
        title = "Frontal Channels Power Spectrum",
        display_plot = false
    )
    
    # With custom parameters
    fig, ax = eegfun.plot_channel_spectrum(dat, 
        channel_selection = eegfun.channels([:Fp1, :Fp2]),
        title = "Custom Power Spectrum",
        x_scale = :log10,
        y_scale = :log10,
        max_freq = 100.0,
        window_size = 512,
        line_width = 3,
        show_freq_bands = false,
        display_plot = true
    )
    
    # Test with different window function
    fig, ax = eegfun.plot_channel_spectrum(dat, 
        channel_selection = eegfun.channels([:Fp1]),
        title = "Power Spectrum with Hamming Window",
        window_function = eegfun.DSP.hamming,
        overlap = 0.75,
        display_plot = true
    )
    
    
    println("✓ Power spectrum plots completed")
end
test_plot_power_spectrum()



###########################
# ERP Plot Tests
###########################
function test_erp_plots()
    println("\n=== Testing ERP Plots ===")
        erp, layout = get_erp_data()
        
        # Single ERP
        fig, ax = eegfun.plot_erp(erp, title = "ERP Waveform")
        
        # Multiple ERPs
        erp2 = deepcopy(erp)
        erp2.condition = 2
        fig, ax = eegfun.plot_erp([erp, erp2], title = "Multiple ERP Conditions")
        
        # With custom styling
        fig, ax = eegfun.plot_erp(erp, 
            title = "Custom ERP Plot",
            color = :red,
            linewidth = 2,
            ylabel = "Amplitude (μV)"
        )
        
        println("✓ ERP plots completed")
end
test_erp_plots()

###########################
# Epochs Plot Tests
###########################
function test_epochs_plots()
    println("\n=== Testing Epochs Plots ===")
    try
        epochs, layout = get_epochs_data()
        
        # Basic epochs plot
        fig, ax = eegfun.plot_epochs(epochs, title = "Epochs Plot")
        
        # With custom styling
        fig, ax = eegfun.plot_epochs(epochs, 
            title = "Custom Epochs Plot",
            color = :blue,
            alpha = 0.7
        )
        
        println("✓ Epochs plots completed")
    catch e
        println("  (Skipping epochs plots - Error: $e)")
    end
end

###########################
# ERP Image Tests
###########################
function test_erp_image()
    println("\n=== Testing ERP Image Plots ===")
    try
        epochs, layout = get_epochs_data()
        
        # Basic ERP image
        fig, ax = eegfun.plot_erp_image(epochs, title = "ERP Image")
        
        # With smoothing
        fig, ax = eegfun.plot_erp_image(epochs, 
            title = "ERP Image with Smoothing",
            boxcar_average = 5,
            time_smoothing = 3
        )
        
        # With custom styling
        fig, ax = eegfun.plot_erp_image(epochs, 
            title = "Custom ERP Image",
            colormap = :viridis,
            ylim = (-50, 50)
        )
        
        println("✓ ERP image plots completed")
    catch e
        println("  (Skipping ERP image plots - Error: $e)")
    end
end

###########################
# Data Browser Tests
###########################
function test_databrowser()
    println("\n=== Testing Data Browser ===")
    try
        dat, layout = get_data()
        
        # Basic data browser
        fig, ax = eegfun.plot_databrowser(dat, title = "Data Browser")
        
        # With custom styling
        fig, ax = eegfun.plot_databrowser(dat, 
            title = "Custom Data Browser",
            channels = [:Fp1, :Fp2, :F3, :F4],
            color = :blue
        )
        
        println("✓ Data browser plots completed")
    catch e
        println("  (Skipping data browser plots - Error: $e)")
    end
end


