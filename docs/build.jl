using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using eegfun

# Set up the documentation
makedocs(
    sitename = "eegfun",
    modules = [eegfun],
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "Public API" => "api.md",
        ],
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
        size_threshold = nothing,  # Disable size threshold
    ),
    doctest = true,
    checkdocs = :all,
)

deploydocs(;
    repo = "github.com/igmmgi/eegfun.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "master"],
    push_preview = true,
)
