# ==============================================================================
# Theme Definitions for EegFun.jl
# ==============================================================================

"""
Default colormap for all EegFun plots. A diverging colormap is appropriate because
the vast majority of EEG visualizations display data centered around zero (topographies,
t-statistics, baseline-corrected TF power, ERP voltages, power differences). The few
truly sequential plots (e.g., RSA dissimilarity, confusion matrices) override this with
their own context-specific colormaps.
"""
const DEFAULT_COLORMAP = :coolwarm

"""
    theme_eegfun()

Returns a `Makie.Theme` optimized for EEG data visualization.

This theme incorporates the previous default keyword arguments of `EegFun.jl`,
providing a clean, publication-ready look out-of-the-box.

# Features
- Base `fontsize` = 16
- Default `colormap` = `:coolwarm`
- Categorical `palette` = Makie's colorblind-friendly `wong_colors`
- Clean `Axis` without gridlines
"""
function theme_eegfun()
    return Theme(
        fontsize = 16,
        linewidth = 2,
        colormap = DEFAULT_COLORMAP,
        palette = (color = Makie.wong_colors(),),
        Axis = (xgridvisible = false, ygridvisible = false),
    )
end
