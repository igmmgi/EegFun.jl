using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))

# TODO: why is this needed 
# GitHub Actions error if not included
# Use CairoMakie for headless documentation builds (no OpenGL required)
# This prevents GLMakie from trying to initialize GLFW/OpenGL on GitHub Actions
using CairoMakie
CairoMakie.activate!()

using EegFun

# Define page structure following Diátaxis framework
pages = [
    "Home" => "index.md",
    "Tutorials" => [
        "Getting Started" => "tutorials/getting-started.md",
        "Basic Preprocessing" => "tutorials/basic-preprocessing.md",
        "ERP Analysis" => "tutorials/erp-analysis.md",
        "ICA Workflow" => "tutorials/ica-workflow.md",
    ],
    "How-To Guides" => [
        "Filter Data" => "how-to/filter-data.md",
        "Create Epochs" => "how-to/create-epochs.md",
        "Topographic Plots" => "how-to/topographic-plots.md",
    ],
    "Explanations" => [
        "Data Structures" => "explanations/data-structures.md",
        "ICA Concepts" => "explanations/ica.md",
        "Statistical Methods" => "explanations/statistics.md",
    ],
    "Reference" => [
        "Overview" => "reference/index.md",
        "Preprocessing" => "reference/preprocessing.md",
        "Types" => "reference/types.md",
        # TODO: Add more reference pages as needed
    ],
]

# Set up the documentation with DocumenterVitepress
makedocs(
    sitename = "EegFun.jl",
    modules = [EegFun],
    pages = pages,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/igmmgi/EegFun.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "igmmgi.github.io/EegFun.jl",
        md_output_path = ".",
        build_vitepress = false,  # TODO: ???
    ),
    warnonly = [:linkcheck, :cross_references, :missing_docs],  # Don't fail on warnings during development
    draft = false,
    source = "src",
    build = "build",
    checkdocs = :all,
)

# Post-Documenter: Fix theme imports for GitHub Actions build
run(`bash docs/fix_theme.sh`) # TODO: what is going on here?

println("\n Documentation build complete!")
println(" Markdown files generated in: docs/build/")

# Deploy configuration
deploydocs(repo = "github.com/igmmgi/EegFun.jl.git", devbranch = "main", push_preview = true)
