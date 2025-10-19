# Example demonstrating the interactive help system
# This example shows how to use the 'i' key to get help in interactive plots

using EEGFun
using Makie

# Create some sample data for demonstration
println("Creating sample EEG data...")

# Simulate some EEG data
n_channels = 32
n_samples = 1000
sample_rate = 250.0
time = collect(0:1/sample_rate:(n_samples-1)/sample_rate)

# Create some realistic-looking EEG data with artifacts
data = randn(n_channels, n_samples) * 10  # Base noise
data[1, :] .+= 50 * sin.(2π * 10 * time)  # Alpha rhythm in first channel
data[5, 200:300] .+= 100  # Artifact in channel 5
data[10, 400:500] .+= 80   # Another artifact

# Create channel names
channel_names = ["Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2",
                 "F7", "F8", "T3", "T4", "T5", "T6", "Fz", "Cz", "Pz", "Oz",
                 "Fpz", "Fcz", "Cpz", "Poz", "Fc1", "Fc2", "Cp1", "Cp2", "P1", "P2",
                 "Po1", "Po2"]

# Create ContinuousData
dat = ContinuousData(data, time, channel_names)

println("Sample data created with $(n_channels) channels and $(n_samples) samples")
println("\n" * "="^60)
println("INTERACTIVE HELP SYSTEM DEMONSTRATION")
println("="^60)
println("\nThis example demonstrates the new help system for interactive plots.")
println("In any interactive plot, press the 'i' key to see available controls.")
println("\nLet's create some plots and try the help system...")

# Example 1: ERP Plot
println("\n1. Creating ERP plot...")
println("   Press 'i' in the plot window to see ERP plot controls")
fig_erp, axes_erp = plot_erp(dat, layout = :grid)

# Example 2: Data Browser
println("\n2. Creating data browser...")
println("   Press 'i' in the plot window to see data browser controls")
fig_browser = plot_databrowser(dat)

# Example 3: Topography Plot
println("\n3. Creating topography plot...")
println("   Press 'i' in the plot window to see topography controls")
fig_topo = plot_topography(dat, time_point = 0.5)

# Example 4: Power Spectrum
println("\n4. Creating power spectrum plot...")
println("   Press 'i' in the plot window to see power spectrum controls")
fig_psd = plot_power_spectrum(dat, channels = 1:5)

println("\n" * "="^60)
println("HELP SYSTEM FEATURES:")
println("="^60)
println("• Press 'i' in any interactive plot to see available controls")
println("• Help information is context-aware for each plot type")
println("• Includes keyboard shortcuts, mouse interactions, and tips")
println("• Press 'i' again to hide the help")
println("\n" * "="^60)
println("Try pressing 'i' in any of the plot windows above!")
println("="^60)

# You can also show help programmatically
println("\nYou can also show help programmatically:")
println("show_plot_help(:erp)  # Show ERP plot help")
println("show_plot_help(:databrowser)  # Show data browser help")
println("show_plot_help(:topography)  # Show topography help")
println("show_plot_help(:power_spectrum)  # Show power spectrum help")
println("show_plot_help(:ica)  # Show ICA help")
println("show_plot_help(:epochs)  # Show epochs help")
println("show_plot_help(:erp_image)  # Show ERP image help")
println("show_plot_help(:triggers)  # Show triggers help")

# Demonstrate programmatic help
println("\nHere's the help for ERP plots:")
show_plot_help(:erp)
