# Shared Makie layout, slider, animation, and biophysics helpers for interactive teaching demos.

######################################
######### Internal Functions #########
######################################

function _demo_adaptive_fonts(fig::Figure, base_width::Real, specs::Vararg{Tuple{Int,Int},N}) where {N}
    obs = ntuple(i -> Observable(specs[i][1]), Val(N))
    on(fig.scene.viewport) do area
        scale = area.widths[1] / base_width
        @inbounds for i = 1:N
            def_val, min_val = specs[i]
            obs[i][] = max(min_val, round(Int, def_val * scale))
        end
    end
    return obs
end

function _demo_adaptive_fonts(fig::Figure, base_width::Real, specs::Vararg{Tuple{Real,Real},N}) where {N}
    specs_f = ntuple(i -> (Float64(specs[i][1]), Float64(specs[i][2])), Val(N))
    obs = ntuple(i -> Observable(specs_f[i][1]), Val(N))
    bw = Float64(base_width)
    on(fig.scene.viewport) do area
        scale = area.widths[1] / bw
        @inbounds for i = 1:N
            def_val, min_val = specs_f[i]
            obs[i][] = Float64(max(round(Int, min_val), round(Int, def_val * scale)))
        end
    end
    return obs
end

function _demo_labelled_slider(
    parent,
    col::Int,
    header,
    range_vals,
    startval,
    fmt,
    fontsize;
    obs = nothing,
    transform = identity,
    on_change = nothing,
    kwargs...,
)
    Label(parent[1, col], header; fontsize = fontsize, halign = :center)
    sl = Slider(parent[2, col]; range = range_vals, startvalue = startval, kwargs...)
    lbl = Label(parent[3, col], fmt(startval); fontsize = fontsize, halign = :center)
    if !isnothing(obs) || !isnothing(on_change)
        on(sl.value) do v
            !isnothing(obs) && (obs[] = transform(v))
            lbl.text = fmt(v)
            !isnothing(on_change) && on_change(v)
        end
    end
    return sl, lbl
end

function _demo_hslider!(
    parent,
    row::Int,
    col_start::Int,
    label_str::String,
    range_vals,
    startval,
    fmt;
    obs = nothing,
    fontsize = 12,
    sl_width = 110,
    val_width = nothing,
    lbl_halign = :center,
    val_halign = :left,
    bold::Bool = true,
)
    Label(parent[row, col_start], text = label_str, font = bold ? :bold : :regular, fontsize = fontsize, halign = lbl_halign)
    sl =
        isnothing(sl_width) ? Slider(parent[row, col_start+1], range = range_vals, startvalue = startval) :
        Slider(parent[row, col_start+1], range = range_vals, startvalue = startval, width = sl_width)
    val_kw = isnothing(val_width) ? (;) : (; width = val_width)
    Label(parent[row, col_start+2], text = lift(fmt, sl.value), fontsize = fontsize, halign = val_halign; val_kw...)
    !isnothing(obs) && on(val -> (obs[] = val), sl.value)
    return sl
end

function _demo_run_async!(
    step_fn!::Function,
    fig::Figure,
    active_ref::Base.RefValue{Bool},
    items;
    dt::Real = 1.0 / 30.0,
    on_finish! = nothing,
)
    active_ref[] && return
    active_ref[] = true
    @async begin
        try
            for item in items
                (!active_ref[] || !events(fig).window_open[]) && break
                step_fn!(item)
                sleep(dt)
            end
            if events(fig).window_open[] && !isnothing(on_finish!)
                on_finish!()
            end
        finally
            active_ref[] = false
        end
    end
end

function _demo_batch_update!(body::Function, guard::Base.RefValue{Bool}, refresh_fn = nothing)
    guard[] = true
    try
        body()
    finally
        guard[] = false
    end
    !isnothing(refresh_fn) && refresh_fn()
    return nothing
end

function _demo_eval_scr_unit(t::Real, t_on::Real, t_r::Real, t_h::Real)
    u = t - t_on
    u <= 0.0 && return 0.0
    u <= t_r && return 0.5 * (1.0 - cos(π * u / t_r))
    a_half = 1.6783469900166605
    x = (u - t_r) / t_h
    return (1.0 + a_half * x) * exp(-a_half * x)
end

function _demo_draw_diff_amp!(ax::Axis, y_amp::Real; half_h::Real = 0.45, out_x::Real = 1.50, label_above::Bool = false)
    ya = Float32(y_amp)
    hh = Float32(half_h)
    poly!(
        ax,
        [Point2f(-0.70f0, ya + hh), Point2f(-0.70f0, ya - hh), Point2f(0.70f0, ya)],
        color = "#f8fafc",
        strokewidth = 2.0,
        strokecolor = "#1e293b",
    )
    text!(ax, -0.46, ya + 0.5f0 * hh, text = "+", align = (:center, :center), font = :bold, fontsize = 17, color = "#dc2626")
    text!(ax, -0.46, ya - 0.5f0 * hh, text = "–", align = (:center, :center), font = :bold, fontsize = 19, color = "#2563eb")
    if label_above
        text!(
            ax,
            0.0,
            ya + hh + 0.65f0,
            text = "Difference Amplifier",
            align = (:center, :center),
            font = :bold,
            fontsize = 16,
            color = "#1e293b",
        )
    else
        text!(
            ax,
            0.0,
            ya - hh - 0.13f0,
            text = "Difference Amplifier",
            align = (:center, :top),
            font = :bold,
            fontsize = 11.0,
            color = "#1e293b",
        )
        text!(
            ax,
            0.0,
            ya - hh - 0.33f0,
            text = "Vout = V(+) – V(-)",
            align = (:center, :top),
            font = :italic,
            fontsize = 8.5,
            color = "#64748b",
        )
    end
    ox = Float32(out_x)
    lines!(ax, [Point2f(0.70f0, ya), Point2f(ox, ya)], color = "#1e293b", linewidth = 2.1)
    scatter!(ax, [Point2f(ox + 0.05f0, ya)], marker = :circle, markersize = 9.5, color = :white, strokecolor = "#1e293b", strokewidth = 2.0)
    return nothing
end
