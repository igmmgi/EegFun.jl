using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))

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
    "Explanations" => ["Data Structures" => "explanations/data-structures.md", "Statistical Methods" => "explanations/statistics.md"],
    "Reference" => [
        "Overview" => "reference/index.md",
        "Types" => "reference/types.md",
        # TODO: Add more reference pages as needed
    ],
]

# Build and deploy documentation
makedocs(
    sitename = "EegFun.jl",
    modules = [EegFun],
    pages = pages,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/igmmgi/EegFun.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://igmmgi.github.io/EegFun.jl",
    ),
    warnonly = [:linkcheck, :cross_references, :missing_docs],
    draft = false,
    source = "src",
    build = "build",
)

DocumenterVitepress.deploydocs(repo = "github.com/igmmgi/EegFun.jl.git", push_preview = true, devbranch = "main", devurl = "dev")
