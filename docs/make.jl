using Documenter

# Include the source files directly for now
# push!(LOAD_PATH, "../src")

makedocs(
    sitename = "EEGfun.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://your-username.github.io/EEGfun.jl",
        edit_link = "main"
    ),
    # modules = [EEGfun],  # Comment out since we can't load the full package
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

# deploydocs(
#     repo = "github.com/your-username/EEGfun.jl.git",
#     push_preview = true
# ) 