This demo shows how to customize plot appearance using Makie themes (`EegFun.theme_eegfun()`, `with_theme`, `Theme`) and function-level keyword arguments.

### Theming Hierarchy

EegFun applies `EegFun.theme_eegfun()` per figure rather than mutating your global Makie theme, with a clear precedence order:

1. **Function keyword arguments (`kwargs`)** — Highest priority; always override active themes (e.g., `color`, `colormap`, `linewidth`, `plot_title_fontsize`).
2. **User themes (`with_theme` / `set_theme!`)** — Override EegFun defaults for colors, palettes, fonts, and backgrounds (e.g., `theme_ggplot2()`, `theme_dark()`, or custom `Theme(...)`).
3. **`EegFun.theme_eegfun()`** — Built-in per-figure defaults when no custom theme or keyword override is provided.

### Common Customizations

- **Color palettes & colormaps**: Set `palette = (color = [...],)` for multi-condition line plots (`plot_erp`) and `colormap = :inferno` for scalp maps (`plot_topography`).
- **Font scaling**: Pass `theme_fontsize` for figure-wide font scaling or `plot_title_fontsize` for title-specific sizing.
- **Publication figures**: Combine a clean `Theme(fontsize = ..., palette = ..., linewidth = ...)` block with layout options (`layout = :grid`, `yreversed = true`, `legend_channel`).

