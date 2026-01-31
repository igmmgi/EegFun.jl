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

# Build and deploy documentation (matching Makie's pattern)
makedocs(
    sitename = "EegFun.jl",
    modules = [EegFun],
    pages = pages,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/igmmgi/EegFun.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://igmmgi.github.io/EegFun.jl",
        build_vitepress = false,  # Patch theme imports before building
    ),
    warnonly = [:linkcheck, :cross_references, :missing_docs],
    draft = false,
    source = "src",
    build = "build",
)

# Fix @/ imports in generated theme (DocumenterVitepress limitation)
theme_file = joinpath(@__DIR__, "build", ".documenter", ".vitepress", "theme", "index.ts")
if isfile(theme_file)
    content = read(theme_file, String)
    content = replace(content, r"\"@/" => "\"./")
    content = replace(content, r"'@/" => "'./")
    write(theme_file, content)

    # Copy Vue components to build directory
    src_theme = joinpath(@__DIR__, "src", ".vitepress", "theme")
    dest_theme = dirname(theme_file)
    for vue_file in ["VersionPicker.vue", "AuthorBadge.vue", "Authors.vue"]
        src = joinpath(src_theme, vue_file)
        dest = joinpath(dest_theme, vue_file)
        if isfile(src)
            cp(src, dest, force = true)
        end
    end

    @info "Fixed theme imports and copied Vue components"

    # Build VitePress site
    @info "Building VitePress site..."
    cd(joinpath(@__DIR__, "build", ".documenter")) do
        run(`npm install`)
        run(`npx vitepress build .`)
    end
    @info "VitePress build complete"
end


DocumenterVitepress.deploydocs(repo = "github.com/igmmgi/EegFun.jl.git", devbranch = "main")
