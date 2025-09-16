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
    # First, we need to create an EpochCondition
    condition = eegfun.EpochCondition(
        name = "test_condition",
        trigger_sequence = [1],  # Simple trigger
        reference_index = 1
    )
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

###########################
# Topography Tests
###########################
function test_topography()
    println("\n=== Testing Topography Plots ===")
    dat, layout = get_data()
    
    # Create some test data for topography
    test_values = randn(length(layout.channels))
    
    # Basic topography
    fig, ax = eegfun.plot_topography(layout, test_values, title = "Random Topography")
    
    # With custom styling
    fig, ax = eegfun.plot_topography(layout, test_values, 
        title = "Custom Topography",
        colormap = :jet,
        gridscale = 100,
        method = :spherical_spline
    )
    
    println("✓ Topography plots completed")
end

###########################
# Filter Tests
###########################
function test_filter()
    println("\n=== Testing Filter Plots ===")
    dat, layout = get_data()
    
    # Create a filter
    filter_design = eegfun.create_filter("lp", "iir", 30.0, dat.sample_rate)
    
    # Plot filter response
    fig, ax = eegfun.plot_filter(filter_design, title = "Filter Response")
    
    # With custom styling
    fig, ax = eegfun.plot_filter(filter_design, 
        title = "Custom Filter Plot",
        color = :red,
        linewidth = 3
    )
    
    println("✓ Filter plots completed")
end

###########################
# Power Spectrum Tests
###########################
function test_power_spectrum()
    println("\n=== Testing Power Spectrum Plots ===")
dat, layout = get_data()

    # Channel spectrum
    fig, ax = eegfun.plot_channel_spectrum(dat, :Fp1, title = "Fp1 Power Spectrum")
    
    # Multiple channels
    fig, ax = eegfun.plot_channel_spectrum(dat, [:Fp1, :Fp2, :F3, :F4], 
        title = "Frontal Channels Power Spectrum")
    
    # Component spectrum (if ICA is available)
    try
        epochs, _ = get_epochs_data()
        ica_result = eegfun.run_ica(epochs, n_components = 10)
        fig, ax = eegfun.plot_component_spectrum(ica_result, 1, title = "Component 1 Spectrum")
    catch
        println("  (Skipping component spectrum - ICA not available)")
    end
    
    println("✓ Power spectrum plots completed")
end

###########################
# ICA Plot Tests
###########################
function test_ica_plots()
    println("\n=== Testing ICA Plots ===")
    try
        epochs, layout = get_epochs_data()
        ica_result = eegfun.run_ica(epochs, n_components = 10)
        
        # Topoplot
        fig, ax = eegfun.plot_ica_topoplot(ica_result, 1, layout, title = "Component 1 Topoplot")
        
        # Component activation
        fig, ax = eegfun.plot_ica_component_activation(ica_result, 1, title = "Component 1 Activation")
        
        # Quality assessment plots
        fig, ax = eegfun.plot_eog_component_features(ica_result, title = "EOG Component Features")
        fig, ax = eegfun.plot_ecg_component_features(ica_result, title = "ECG Component Features")
        fig, ax = eegfun.plot_spatial_kurtosis_components(ica_result, title = "Spatial Kurtosis")
        fig, ax = eegfun.plot_line_noise_components(ica_result, title = "Line Noise Components")
        
        println("✓ ICA plots completed")
    catch e
        println("  (Skipping ICA plots - Error: $e)")
    end
end

###########################
# ERP Plot Tests
###########################
function test_erp_plots()
    println("\n=== Testing ERP Plots ===")
    try
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
    catch e
        println("  (Skipping ERP plots - Error: $e)")
    end
end

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

###########################
# Run All Tests
###########################
function test_all_plots()
    println("\n🚀 Running All Plot Tests...")
    
    test_channel_summary()
    test_correlation_heatmap()
    test_joint_probability()
    test_layout_plots()
    test_topography()
    test_filter()
    test_power_spectrum()
    test_ica_plots()
    test_erp_plots()
    test_epochs_plots()
    test_erp_image()
    test_databrowser()
    
    println("\n🎉 All plot tests completed!")
end

# Run a quick test by default
println("\nRunning quick test...")
test_channel_summary()
test_correlation_heatmap()
test_joint_probability()
