# Visualization subsystem, layout engine, interactive GUIs, and Makie recipe includes.

# Core plot utilities, themes, layout engine, and shared interactivity (loaded first)
include("plots/utils/plot_misc.jl")
include("plots/utils/themes.jl")
include("plots/utils/plot_ui_utils.jl")
include("plots/utils/layout_system.jl")
include("plots/utils/shared_interactivity.jl")

# Feature plotting modules
include("plots/processing.jl")
include("plots/erps.jl")
include("plots/time_frequency.jl")
include("plots/decoding/plot_decoding.jl")
include("plots/rsa/plot_rsa.jl")

# Launcher GUI, shortcut help system, and Makie type recipes
include("plots/plot_gui.jl")
include("plots/plot_gui_helpers.jl")
include("plots/utils/help_system.jl")
include("plots/utils/makie_recipes.jl")
