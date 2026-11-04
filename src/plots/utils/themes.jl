const DEFAULT_COLORMAP = :coolwarm

"""
    theme_eegfun()

Returns a `Makie.Theme` optimized for EEG data visualization.

# Features
- Base `fontsize` = 16, `linewidth` = 2
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
    isdefined(Makie, :MAKIE_DEFAULT_THEME) || return theme_eegfun()
    overrides = _user_theme_overrides(Makie.current_default_theme(), Makie.MAKIE_DEFAULT_THEME)
    return merge(overrides, theme_eegfun())
end

_eegfun_figure(args...; kwargs...) = Makie.with_theme(() -> Figure(args...; kwargs...), _eegfun_effective_theme())
