# Demo: Plot Styling and Theming
# Shows how to use Makie themes and EegFun kwargs to customize plots.

using EegFun
using GLMakie

# ==============================================================================
# 1. Setup & Data Loading (same as other tutorials)
# ==============================================================================
# Load data
dat = EegFun.read_raw_data(EegFun.example_path("data/bdf/example1.bdf"))
layout = EegFun.read_layout(EegFun.example_path("layouts/biosemi/biosemi72.csv"))
EegFun.polar_to_cartesian_xy!(layout)
dat = EegFun.create_eegfun_data(dat, layout)

# Minimal preprocessing
EegFun.rereference!(dat, :avg)
EegFun.highpass_filter!(dat, 1.0)

# Create some epochs for ERP plots
epoch_cfg = [
    EegFun.EpochCondition(name = "Cond1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Cond2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "Cond3", trigger_sequences = [[3]]),
]
epochs = EegFun.extract_epochs(dat, epoch_cfg, (-2, 4))
erps = EegFun.average_epochs(epochs)

# ==============================================================================
# 2. Global vs. Local Themes
# ==============================================================================

# By default, when you run `using EegFun`, it automatically applies `EegFun.theme_eegfun()`!
EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]))

# However, you can still use `with_theme` to apply a different theme only to a specific block:
with_theme(theme_ggplot2()) do
    # This ERP plot will have the ggplot2 grey background and default colors
    EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]))
end

# ==============================================================================
# 3. Mixing Themes and Custom Palettes (Colormaps/Colors)
# ==============================================================================

# You can merge themes together. Here we use the base eegfun theme, but we override its default 
# color cycle (palette) for lines, and its colormap for topographies.
my_custom_colors = Theme(
    palette = (color = [:purple, :teal, :gold],), # Used by plot_erp
    colormap = :inferno                           # Used by plot_topography
)

with_theme(my_custom_colors) do
    # This ERP plot will now automatically use purple, teal, and gold for Cond 1, 2, and 3
    EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]))
    
    # This topoplot will automatically use the :inferno colormap from the theme
    EegFun.plot_topography(dat, interval_selection = EegFun.times(6), ylim = (-200, 200))
end

# ==============================================================================
# 4. Keyword Arguments (kwargs) Overriding Themes
# ==============================================================================

# Kwargs passed directly to EegFun functions ALWAYS take precedence over the theme.
with_theme(theme_dark()) do
    # Even though theme_dark has its own colors, passing `color` forces it to use these.
    EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]),
        color = [:red, :white, :blue], 
        linewidth = 3,
        plot_title = "My Override Plot")

    # Even though we are in a theme, we can force a different colormap directly via kwargs
    EegFun.plot_topography(dat, interval_selection = EegFun.times(6), 
        colormap = :inferno, ylim = (-200, 200))
end

# ==============================================================================
# 5. Tweaking Titles and Fonts
# ==============================================================================

# If you want to increase the base font size for an entire figure without using 
# a full theme block, you can use `theme_fontsize`:
EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]), 
    theme_fontsize = 24)

# You can also set specific overrides for elements like the main plot title:
EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]),
    plot_title = "Large Custom Title", 
    plot_title_fontsize = 36)

# ==============================================================================
# 6. Real-world example: A Publication-Ready Plot
# ==============================================================================

# Combining a clean, minimal theme with specific kwargs to make a nice figure
publication_theme = Theme(
    fontsize = 16,
    palette = (color = [:black, :dodgerblue, :firebrick],),
    linewidth = 2
)

with_theme(publication_theme) do
    EegFun.plot_erp(
        erps, 
        layout = :grid, 
        channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :Fz, :F4]),
        legend_channel = [:Fz],
        yreversed = true,                   # Negative up
        plot_title = "Publication ERPs",    
        plot_title_fontsize = 20
    )
end
