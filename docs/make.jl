using Documenter
using DocumenterVitepress

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))

using EegFun

# # Define page structure following Diátaxis framework
# pages = [
#     "Home" => "index.md",
#     "Tutorials" => [
#         "Getting Started" => "tutorials/getting-started.md",
#         "Basic Preprocessing" => "tutorials/basic-preprocessing.md",
#         "ERP Analysis" => "tutorials/erp-analysis.md",
#         "ICA Workflow" => "tutorials/ica-workflow.md",
#     ],
#     "How-To Guides" => [
#         "Filter Data" => "how-to/filter-data.md",
#         "Create Epochs" => "how-to/create-epochs.md",
#         "Topographic Plots" => "how-to/topographic-plots.md",
#     ],
#     "Explanations" => ["Data Structures" => "explanations/data-structures.md", "Statistical Methods" => "explanations/statistics.md"],
#     "Reference" => [
#         "Overview" => "reference/index.md",
#         "Types" => "reference/types.md",
#         # TODO: Add more reference pages as needed
#     ],
# ]

# Create version directory structure for DocumenterVitepress (workaround for siteinfo.js bug)
# mkpath(joinpath(@__DIR__, "build", "1"))

# # Build and deploy documentation (SIMPLE - like Makie!)
# makedocs(
#     sitename = "EegFun.jl",
#     modules = [EegFun],
#     pages = pages,
#     format = DocumenterVitepress.MarkdownVitepress(
#         repo = "github.com/igmmgi/EegFun.jl",
#         devbranch = "main",
#         devurl = "dev",
#         deploy_url = "https://igmmgi.github.io/EegFun.jl",
#     ),
#     warnonly = [:linkcheck, :cross_references, :missing_docs],
#     draft = false,
#     source = "src",
# )

# Pre-create versioned build directory (workaround for DocumenterVitepress siteinfo.js bug)
# mkpath(joinpath(@__DIR__, "build", "1"))

makedocs(;
    modules = [EegFun],
    authors = "igmmgi",
    repo = "https://github.com/igmmgi/EegFun.jl",
    sitename = "EegFun.jl",
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/igmmgi/EegFun.jl",
        devbranch = "main", # or master, trunk, ...
        devurl = "dev",
        deploy_url = "https://igmmgi.github.io/EegFun.jl",
    ),
    warnonly = [:linkcheck, :cross_references, :missing_docs],
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
    ],
)

# Deploy built VitePress site (use Documenter.deploydocs to avoid DocumenterVitepress versioning bugs)
Documenter.deploydocs(
    repo = "github.com/igmmgi/EegFun.jl",
    # target = joinpath(@__DIR__, "build", ".documenter", ".vitepress", "dist"),
    # branch = "gh-pages",
    # devbranch = "main",
    push_preview = true,
)
