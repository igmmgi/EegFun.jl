using Documenter

# Add the parent directory to the load path
push!(LOAD_PATH, "..")

using eegfun

makedocs(
    sitename = "eegfun",
    format = Documenter.HTML(
        size_threshold = nothing,  # Disable size threshold
    ),
    modules = [eegfun],
    checkdocs = :none,  # Disable strict documentation checking
    warnonly = [:cross_references],  # Only warn about cross-references, don't fail
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
