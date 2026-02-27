# Script to generate demo screenshots for documentation
# Run locally with CairoMakie — NOT called by CI
#
# Usage:
#   julia --project=. docs/generate_demo_screenshots.jl
#
# Or via the doc_manager.jl interactive menu (option 8)

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using EegFun
using CairoMakie
CairoMakie.activate!()

# Output directory for demo screenshots
const IMAGES_DIR = joinpath(@__DIR__, "..", "images", "demos")
mkpath(IMAGES_DIR)

"""
    save_demo_figure(filename, fig; size=(800, 500))

Save a figure to the demo images directory.
"""
function save_demo_figure(filename::String, fig::Figure; size::Tuple{Int,Int} = (800, 500))
    filepath = joinpath(IMAGES_DIR, filename)
    save(filepath, fig; size = size)
    println("  ✓ Saved: $filepath")
    return filepath
end


# =====================================================================
# Demo screenshots — add new entries here
# =====================================================================

"""
    generate_all_screenshots()

Generate all demo screenshots. Add new screenshot functions to this list.
"""
function generate_all_screenshots()
    println("=== Generating Demo Screenshots ===\n")

    screenshots = [("plot_erp", generate_plot_erp)]

    generated = 0
    failed = 0

    for (name, func) in screenshots
        print("Generating: $name... ")
        try
            func()
            generated += 1
        catch e
            println("  ✗ FAILED: $e")
            failed += 1
        end
    end

    println("\n=== Done: $generated generated, $failed failed ===")
    println("Images saved to: $IMAGES_DIR")
end


# =====================================================================
# Individual screenshot generators
# =====================================================================

function generate_plot_erp()
    # Create synthetic ERP data
    dat = EegFun.create_test_erp_data()

    # Basic ERP plot
    result = EegFun.plot_erp([dat], channel_selection = EegFun.channels(:Ch1), display_plot = false)
    save_demo_figure("plot_erp_basic.png", result.fig)

    # Multi-channel grid layout
    result = EegFun.plot_erp([dat], channel_selection = EegFun.channels(:Ch1, :Ch2, :Ch3), layout = :grid, display_plot = false)
    save_demo_figure("plot_erp_grid.png", result.fig; size = (900, 700))
end


# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_all_screenshots()
end
