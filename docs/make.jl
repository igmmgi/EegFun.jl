using Documenter

# Include the source files directly for now
push!(LOAD_PATH, "../src")

# Load the main module
include("../src/eegfun.jl")

makedocs(
    sitename = "eegfun",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://igmmgi.github.io/eegfun",
        edit_link = "main",
        size_threshold = nothing  # Disable size threshold checking
    ),
    modules = [eegfun],  # Enable automatic docstring extraction
    pages = [
        "Home" => "index.md",
        "Data Loading & Types" => "data_types.md",
        "Preprocessing" => "preprocessing.md",
        "Analysis" => "analysis.md",
        "ICA" => "ica.md",
        "Plotting" => "plotting.md",
        "Utilities" => "utilities.md"
    ],
    doctest = false,  # Disable doctests since we can't load the package
    clean = true,
    checkdocs = :none,  # Disable documentation checking entirely
    linkcheck = false,  # Disable link checking
    warnonly = [:cross_references]  # Only warn about cross-references, don't fail
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "github.com/igmmgi/eegfun.git",
    push_preview = true
)=# 