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

"""
    _user_theme_overrides(current, default)

Return the subset of `current` (the active Makie theme) whose values differ from `default`
(Makie's pristine default theme), recursing into nested attributes. This is what the user
has customised via `set_theme!`, `update_theme!` or `with_theme`.
"""
function _user_theme_overrides(current::Makie.Attributes, default::Makie.Attributes)
    out = Makie.Attributes()
    for (key, val) in current
        if haskey(default, key)
            dval = default[key]
            if val isa Makie.Attributes && dval isa Makie.Attributes
                sub = _user_theme_overrides(val, dval)
                isempty(sub) || (out[key] = sub)
            elseif !isequal(Makie.to_value(val), Makie.to_value(dval))
                out[key] = Makie.to_value(val)
            end
        else
            out[key] = val isa Makie.Attributes ? val : Makie.to_value(val)
        end
    end
    return out
end

function _eegfun_effective_theme()
    # MAKIE_DEFAULT_THEME is internal to Makie; if it ever disappears, fall back to the plain EegFun theme.
    isdefined(Makie, :MAKIE_DEFAULT_THEME) || return theme_eegfun()
    overrides = _user_theme_overrides(Makie.current_default_theme(), Makie.MAKIE_DEFAULT_THEME)
    return merge(overrides, theme_eegfun())  # first argument takes priority
end

"""
    _eegfun_figure(args...; kwargs...)

Create a `Makie.Figure` styled with [`theme_eegfun`](@ref), without touching the global Makie theme.

Makie captures the active theme in the figure's scene at construction time, so axes, blocks and
plots subsequently added to the figure inherit it. Precedence is:
Makie defaults < `theme_eegfun()` < anything the user customised (`set_theme!`, `with_theme`, ...).

All EegFun plotting functions should create figures via this helper rather than `Figure(...)`.
"""
_eegfun_figure(args...; kwargs...) = Makie.with_theme(() -> Figure(args...; kwargs...), _eegfun_effective_theme())
