using Documenter

# Add the parent directory to the load path so we can load the local package
push!(LOAD_PATH, dirname(@__DIR__))
using eegfun

# Set up the documentation
makedocs(
    sitename = "eegfun",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
        size_threshold = nothing,  # Disable size threshold
    ),
    modules = [eegfun],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    doctest = true,
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/igmmgi/eegfun.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "master"],
    push_preview = true,
)
