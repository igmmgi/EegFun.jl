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

function _demo_labelled_slider(parent, col::Int, header, range_vals, startval, fmt, fontsize; kwargs...)
    Label(parent[1, col], header; fontsize = fontsize, halign = :center)
    sl = Slider(parent[2, col]; range = range_vals, startvalue = startval, kwargs...)
    lbl = Label(parent[3, col], fmt(startval); fontsize = fontsize, halign = :center)
    return sl, lbl
end

"""
    signal_example_sampling()

Interactive Nyquist–Shannon sampling theorem demonstration.

Displays a continuous sine wave alongside optional linear (join-the-dots) and ideal
Whittaker–Shannon `sinc` reconstructions as duration, frequency, phase, and sample rate vary.

# Example
```julia
using GLMakie, EegFun
EegFun.signal_example_sampling()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax::Axis`: The axis object containing the plots
"""
function signal_example_sampling()
    _check_interactive_demo_backend("Signal Sampling demo") || return nothing

    fig = _eegfun_figure(size = (1350, 800), title = "Signal Example 1")
    title_font, label_font, tick_font, slider_font, slider_width, slider_height =
        _demo_adaptive_fonts(fig, 1350, (20, 14), (16, 12), (13, 10), (16, 12), (250, 180), (22, 15))

    ax = Axis(
        fig[1, 1:6],
        title = "Interactive Signal Generator",
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        xgridvisible = true,
        ygridvisible = true,
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
    )

    sig_dur = Observable(1.0)
    freq = Observable(10.0)
    phase = Observable(0.0)
    samp_rate = Observable(100.0)
    show_sampled = Observable(false)
    show_sinc = Observable(false)

    base_time = Observable(Float64[])
    samp_time = Observable(Float64[])
    recon_time = Observable(Float64[])
    my_sin1 = Observable(Float64[])
    my_sin2 = Observable(Float64[])
    my_sin3 = Observable(Float64[])

    function sinc_reconstruct(t_samp::Vector{Float64}, y_samp::Vector{Float64}, t_out::AbstractRange, inv_dt::Float64)
        out = Vector{Float64}(undef, length(t_out))
        @inbounds for i in eachindex(t_out)
            t = t_out[i]
            acc = 0.0
            for j in eachindex(t_samp)
                acc += y_samp[j] * sinc((t - t_samp[j]) * inv_dt)
            end
            out[i] = acc
        end
        return out
    end

    function update_signal()
        dur = sig_dur[]
        sr = samp_rate[]
        f = freq[]
        ph = phase[]

        base_t = 0:0.0001:dur
        samp_t = 0:(1/sr):dur
        bt = collect(base_t)
        st = collect(samp_t)

        base_time[] = bt
        samp_time[] = st
        my_sin1[] = @. sin(2π * f * bt + ph)
        s2 = @. sin(2π * f * st + ph)
        my_sin2[] = s2

        t_recon = range(0.0, dur, length = 2000)
        recon_time[] = collect(t_recon)
        my_sin3[] = sinc_reconstruct(st, s2, t_recon, sr)

        xlims!(ax, 0, dur)
        ylims!(ax, -1.2, 1.2)
    end

    for obs in (sig_dur, freq, phase, samp_rate)
        on(obs) do _
            update_signal()
        end
    end
    update_signal()

    lines!(ax, base_time, my_sin1, color = :blue, linewidth = 2, label = "Signal")
    sampled_line = lines!(ax, samp_time, my_sin2, color = :red, linewidth = 2, label = "Linear", visible = false)
    sample_points = scatter!(ax, samp_time, my_sin2, color = :red, markersize = 16, visible = false)
    sinc_line = lines!(ax, recon_time, my_sin3, color = :green, linewidth = 2, label = "Sinc", visible = false)
    sinc_points = scatter!(ax, samp_time, my_sin2, color = :green, markersize = 16, visible = false)

    on(sampled_line.visible) do val
        sample_points.visible = val
    end
    on(show_sampled) do val
        sampled_line.visible = val
    end
    on(sinc_line.visible) do val
        sinc_points.visible = val
    end
    on(show_sinc) do val
        sinc_line.visible = val
    end

    axislegend(ax, position = :rt)

    slider_layout = GridLayout(fig[2, 1:6], tellheight = false, valign = :top, padding = (1, 5, 1, 5), rowgap = 1, colgap = 5)
    slider_configs = (
        (range = 1.0:1.0:10.0, startvalue = 1.0, obs = sig_dur, fmt = x -> "Duration: $(round(x, digits=1)) secs"),
        (range = 1.0:1.0:100.0, startvalue = 10.0, obs = freq, fmt = x -> "Frequency: $(round(Int, x)) Hz"),
        (range = (-π):(π/4):π, startvalue = 0.0, obs = phase, fmt = x -> "Phase Angle: $(round(x, digits=2))"),
        (range = 1.0:1.0:300.0, startvalue = 100.0, obs = samp_rate, fmt = x -> "Sample Rate: $(round(Int, x)) Hz"),
    )

    for (i, cfg) in enumerate(slider_configs)
        sl = Slider(slider_layout[1, i], range = cfg.range, startvalue = cfg.startvalue, width = slider_width, height = slider_height)
        lbl = Label(slider_layout[2, i], text = cfg.fmt(cfg.startvalue), width = slider_width, fontsize = slider_font)
        on(sl.value) do val
            cfg.obs[] = val
            lbl.text = cfg.fmt(val)
        end
    end

    show_sampled_checkbox = Checkbox(slider_layout[1, 5], checked = false, halign = :center)
    Label(slider_layout[2, 5], text = "Linear", fontsize = slider_font, halign = :center)
    on(show_sampled_checkbox.checked) do val
        show_sampled[] = val
    end

    show_sinc_checkbox = Checkbox(slider_layout[1, 6], checked = false, halign = :center)
    Label(slider_layout[2, 6], text = "Sinc", fontsize = slider_font, halign = :center)
    on(show_sinc_checkbox.checked) do val
        show_sinc[] = val
    end

    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))

    _display_demo_figure(fig)
    return fig, ax
end
