#!/usr/bin/env julia

"""
Quick script to plot all EEG layout files and save images to their respective folders.

This script:
1. Finds all CSV layout files in the data/layouts directory
2. Reads each layout file using the eegfun package
3. Plots the layout in 2D using Makie
4. Saves the plot as a PNG image in the same folder as the layout file

Usage:
    julia plot_layouts.jl
"""

using Pkg
Pkg.activate("../../")

using eegfun
using CairoMakie

function main()
    println("Starting layout plotting script...")
    
    # Define the base directory for layouts (current directory)
    layouts_dir = "."
    
    if !isdir(layouts_dir)
        error("Layouts directory not found: $layouts_dir")
    end
    
    # Find all CSV files in subdirectories
    csv_files = String[]
    for (root, dirs, files) in walkdir(layouts_dir)
        for file in files
            if endswith(file, ".csv")
                push!(csv_files, joinpath(root, file))
            end
        end
    end
    
    if isempty(csv_files)
        println("No CSV layout files found in $layouts_dir")
        return
    end
    
    println("Found $(length(csv_files)) layout files:")
    for file in csv_files
        println("  - $file")
    end
    println()
    
    # Process each layout file
    for (i, csv_file) in enumerate(csv_files)
        try
            println("Processing $(i)/$(length(csv_files)): $csv_file")
            
            # Read the layout and convert to 2D, and plot
            layout = eegfun.read_layout(csv_file)
            eegfun._ensure_coordinates_2d!(layout)
            fig, ax = eegfun.plot_layout_2d(layout; display_plot=false)
            
            # Save the plot using CairoMakie with specific pixel dimensions
            output_file = replace(csv_file, ".csv" => ".png")
            CairoMakie.save(output_file, fig; size=(1000, 1000), px_per_unit=1)
            println("  ✓ Saved: $output_file")
            
        catch e
            println("  ✗ Error processing $csv_file: $e")
            continue
        end
    end
    
    println("\nLayout plotting completed!")
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end