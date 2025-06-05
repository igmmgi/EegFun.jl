using Documenter

# Include the source files directly for now
# push!(LOAD_PATH, "../src")

makedocs(
    sitename = "eegfun",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://igmmgi.github.io/eegfun",
        edit_link = "main"
    ),
    # modules = [eegfun],  # Comment out since we can't load the full package
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Data Types" => "api/types.md",
            "Preprocessing" => "api/preprocessing.md",
            "Analysis" => "api/analysis.md",
            "Utilities" => "api/utilities.md",
            "Plotting" => "api/plotting.md"
        ],
        "Examples" => "examples.md"
    ],
    doctest = false,  # Disable doctests since we can't load the package
    clean = true,
    # checkdocs = :exports  # Comment out since we're not checking exports
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "github.com/igmmgi/eegfun.git",
    push_preview = true
)=# 