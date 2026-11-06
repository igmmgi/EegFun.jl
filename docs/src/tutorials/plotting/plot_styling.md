# Plot Styling

This demo shows how to customize plot appearance using Makie themes (`EegFun.theme_eegfun()`, `with_theme`, `Theme`) and function-level keyword arguments.

## Theming Hierarchy

EegFun applies `EegFun.theme_eegfun()` per figure rather than mutating your global Makie theme, with a clear precedence order:

1. **Function keyword arguments (`kwargs`)** — Highest priority; always override active themes (e.g., `color`, `colormap`, `linewidth`, `plot_title_fontsize`).
2. **User themes (`with_theme` / `set_theme!`)** — Override EegFun defaults for colors, palettes, fonts, and backgrounds (e.g., `theme_ggplot2()`, `theme_dark()`, or custom `Theme(...)`).
3. **`EegFun.theme_eegfun()`** — Built-in per-figure defaults when no custom theme or keyword override is provided.

## Common Customizations

- **Color palettes & colormaps**: Set `palette = (color = [...],)` for multi-condition line plots (`plot_erp`) and `colormap = :inferno` for scalp maps (`plot_topography`).
- **Font scaling**: Pass `theme_fontsize` for figure-wide font scaling or `plot_title_fontsize` for title-specific sizing.
- **Publication figures**: Combine a clean `Theme(fontsize = ..., palette = ..., linewidth = ...)` block with layout options (`layout = :grid`, `yreversed = true`, `legend_channel`).



## Code Examples

::: details Show Code

```julia
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

# Figures created by EegFun plotting functions use `EegFun.theme_eegfun()` by default.
# This is applied per figure, so `using EegFun` does not change your global Makie theme.
EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]))

# Anything you customise yourself (via `with_theme`, `set_theme!` or `update_theme!`) takes
# precedence over the EegFun defaults. For example, apply a different theme to a specific block:
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
    colormap = :inferno,                           # Used by plot_topography
)

with_theme(my_custom_colors) do
    # This ERP plot will now automatically use purple, teal, and gold for Cond 1, 2, and 3
    EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]))

    # This topoplot will automatically use the :inferno colormap from the theme
    EegFun.plot_topography(dat, interval_selection = EegFun.times(6), colorrange = (-200, 200))
end

# ==============================================================================
# 4. Keyword Arguments (kwargs) Overriding Themes
# ==============================================================================

# Kwargs passed directly to EegFun functions ALWAYS take precedence over the theme.
with_theme(theme_dark()) do
    # Even though theme_dark has its own colors, passing `color` forces it to use these.
    EegFun.plot_erp(
        erps,
        layout = :single,
        channel_selection = EegFun.channels([:Fp1]),
        color = [:red, :white, :blue],
        linewidth = 3,
        plot_title = "My Override Plot",
    )

    # Even though we are in a theme, we can force a different colormap directly via kwargs
    EegFun.plot_topography(dat, interval_selection = EegFun.times(6), colormap = :inferno, colorrange = (-200, 200))
end

# ==============================================================================
# 5. Tweaking Titles and Fonts
# ==============================================================================

# If you want to increase the base font size for an entire figure without using 
# a full theme block, you can use `theme_fontsize`:
EegFun.plot_erp(erps, layout = :single, channel_selection = EegFun.channels([:Fp1]), theme_fontsize = 24)

# You can also set specific overrides for elements like the main plot title:
EegFun.plot_erp(
    erps,
    layout = :single,
    channel_selection = EegFun.channels([:Fp1]),
    plot_title = "Large Custom Title",
    plot_title_fontsize = 36,
)

# ==============================================================================
# 6. Real-world example: A Publication-Ready Plot
# ==============================================================================

# Combining a clean, minimal theme with specific kwargs to make a nice figure
publication_theme = Theme(fontsize = 16, palette = (color = [:black, :dodgerblue, :firebrick],), linewidth = 2)

with_theme(publication_theme) do
    EegFun.plot_erp(
        erps,
        layout = :grid,
        channel_selection = EegFun.channels([:Fp1, :Fp2, :F3, :Fz, :F4]),
        legend_channel = [:Fz],
        yreversed = true,                   # Negative up
        plot_title = "Publication ERPs",
        plot_title_fontsize = 20,
    )
end
```

:::

## See Also

- [API Reference](../../reference/index.md)
