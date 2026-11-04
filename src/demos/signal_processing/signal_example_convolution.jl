"""
    signal_example_convolution()

Interactive discrete convolution demo showing Gaussian, Boxcar, and Morlet wavelet kernels
sliding across a noisy sinusoid and accumulating the output or amplitude envelope.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_convolution()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::Tuple`: `(ax_sig, ax_out)` — the two axis objects
"""
function signal_example_convolution()
    _check_interactive_demo_backend("Signal Convolution demo") || return nothing

    FS = 512.0
    EPOCH = 2.0

    fig = _eegfun_figure(size = (1350, 650), figure_padding = (8, 8, 8, 8), title = "Convolution Demo")
    lbl_sz, tick_sz, title_sz, ctrl_sz = _demo_adaptive_fonts(fig, 1350, (15, 11), (12, 10), (18, 13), (15, 11))

    ax_common = (
        xlabel = "Time (s)",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible = true,
        ygridvisible = true,
    )
    ax_sig = Axis(fig[1, 1]; title = "Signal + Kernel at current position", ylabel = "Amplitude", ax_common...)
    ax_out = Axis(fig[2, 1]; title = "Convolution Output", ylabel = "Output", ax_common...)
    linkxaxes!(ax_sig, ax_out)

    kernel_pos = Observable(0.0)
    kernel_width = Observable(0.1)
    kernel_type = Observable(:gaussian)
    wavelet_freq = Observable(8.0)
    wavelet_cycles = Observable(5.0)
    sig_freq = Observable(4.0)
    noise_lvl = Observable(0.3)

    kernel_options = ["Gaussian", "Boxcar", "Wavelet (Morlet)"]

    t_arr = collect(range(0.0, EPOCH - 1.0 / FS, step = 1.0 / FS))
    n = length(t_arr)
    base_noise = randn(n)

    sig_data = Observable(Float64[])
    kern_t = Observable(Float64[])
    kern_vals = Observable(Float64[])
    shd_x = Observable(Float64[])
    out_t = Observable(Float64[])
    out_data = Observable(Float64[])

    function build_kernel(type::Symbol, width_s::Float64, wfreq::Float64, ncycles::Float64)
        span_s = type === :wavelet ? ncycles / wfreq : width_s
        half = round(Int, span_s / 2.0 * FS)
        σ = span_s / 6.0
        τ = ((-half):half) ./ FS

        if type === :gaussian
            k = @. exp(-(τ^2) / (2.0 * σ^2))
            k ./= sum(k)
            return k, nothing, τ
        elseif type === :boxcar
            m = length(τ)
            return fill(1.0 / m, m), nothing, τ
        else
            gauss = @. exp(-(τ^2) / (2.0 * σ^2))
            inv_norm = 2.0 / sum(gauss)
            k_r = @. gauss * cos(2π * wfreq * τ) * inv_norm
            k_i = @. gauss * sin(2π * wfreq * τ) * inv_norm
            return k_r, k_i, τ
        end
    end

    function full_conv(sig::Vector{Float64}, k::Vector{Float64})
        ns = length(sig)
        nk = length(k)
        out = Vector{Float64}(undef, ns)
        h = (nk - 1) ÷ 2
        @inbounds for i = 1:ns
            acc = 0.0
            for j = 1:nk
                idx = i - h + j - 1
                if 1 <= idx <= ns
                    acc = muladd(sig[idx], k[nk-j+1], acc)
                end
            end
            out[i] = acc
        end
        return out
    end

    function update()
        sf = sig_freq[]
        nl = noise_lvl[]
        sig = nl > 0.0 ? (@. sin(2π * sf * t_arr) + nl * base_noise) : (@. sin(2π * sf * t_arr))

        ktype = kernel_type[]
        k_r, k_i, τ = build_kernel(ktype, kernel_width[], wavelet_freq[], wavelet_cycles[])

        pos_idx = clamp(round(Int, kernel_pos[] * EPOCH * FS), 1, n)

        if k_i !== nothing
            out_r = full_conv(sig, k_r)
            out_i = full_conv(sig, k_i)
            full_out = @. hypot(out_r, out_i)
        else
            full_out = full_conv(sig, k_r)
        end

        kern_centre = t_arr[pos_idx]
        k_t = collect(kern_centre .+ τ)
        inv_peak = 1.0 / max(maximum(abs, k_r), 1e-8)
        k_v = k_r .* inv_peak

        half = (length(k_r) - 1) ÷ 2
        i_lo = max(1, pos_idx - half)
        i_hi = min(n, pos_idx + half)

        sig_data[] = sig
        kern_t[] = k_t
        kern_vals[] = k_v
        shd_x[] = t_arr[i_lo:i_hi]
        out_t[] = t_arr[1:pos_idx]
        out_data[] = full_out[1:pos_idx]

        ax_out.title =
            ktype === :wavelet ? "Amplitude Envelope — |Wavelet convolution at $(round(Int, wavelet_freq[])) Hz|" :
            "Convolution Output (built up as kernel slides)"

        xlims!(ax_sig, 0.0, EPOCH)
        xlims!(ax_out, 0.0, EPOCH)
        lim = 1.5 + nl * 1.5
        ylims!(ax_sig, -lim, lim)
        if ktype === :wavelet
            ylims!(ax_out, 0.0, 1.1)
        else
            max_o = max(maximum(abs, full_out), 0.01)
            ylims!(ax_out, -max_o * 1.2, max_o * 1.2)
        end
    end

    is_applying_preset = Ref(false)
    for obs in (kernel_pos, kernel_width, kernel_type, wavelet_freq, wavelet_cycles, sig_freq, noise_lvl)
        on(obs) do _
            is_applying_preset[] || update()
        end
    end
    update()

    lines!(ax_sig, t_arr, sig_data, color = :black, linewidth = 2)
    band!(ax_sig, shd_x, @lift(fill(-3.0, length($shd_x))), @lift(fill(3.0, length($shd_x))), color = (:black, 0.07))
    band!(ax_sig, kern_t, @lift(zeros(length($kern_t))), kern_vals, color = (:black, 0.20))
    lines!(ax_sig, kern_t, kern_vals, color = :black, linewidth = 2, label = "Kernel")
    vlines!(ax_sig, @lift([$kernel_pos * EPOCH]), color = :black, linestyle = :dash, linewidth = 1.5)
    axislegend(ax_sig, position = :rt)

    lines!(ax_out, out_t, out_data, color = :black, linewidth = 2)
    hlines!(ax_out, [0.0], color = :black, linewidth = 0.8, linestyle = :dash)
    vlines!(ax_out, @lift([$kernel_pos * EPOCH]), color = :black, linestyle = :dash, linewidth = 1.5)

    out_text = @lift begin
        pos_s = $kernel_pos * EPOCH
        "Output at t = $(round(pos_s, digits=2)) s = $(isempty($out_data) ? "–" : round($out_data[end], digits=3))"
    end
    text!(ax_out, 0.5, 0.97, text = out_text, space = :relative, fontsize = lbl_sz, color = :black, align = (:center, :top))

    ctrl = GridLayout(fig[3, 1], colgap = 16)
    Label(ctrl[1, 1:4], "Kernel Position", fontsize = ctrl_sz, halign = :center)
    sl_pos = Slider(ctrl[2, 1:4], range = 0.0:0.001:1.0, startvalue = 0.0)
    lbl_pos = Label(ctrl[3, 1:4], "pos: 0.00 s", fontsize = ctrl_sz, halign = :center)

    inner = GridLayout(ctrl[4, 1:4], colgap = 16)
    sl_kw, lbl_kw = _demo_labelled_slider(inner, 1, "Kernel Width", 0.01:0.01:2.0, 0.1, v -> "$(round(Int, v * 1000)) ms", ctrl_sz)
    sl_cy, lbl_cy = _demo_labelled_slider(inner, 2, "Cycles", 1.0:1.0:10.0, 5.0, v -> "$(round(Int, v))", ctrl_sz)
    sl_wf, lbl_wf = _demo_labelled_slider(inner, 3, "Wavelet Freq", 1.0:1.0:40.0, 8.0, v -> "$(round(Int, v)) Hz", ctrl_sz)
    sl_sf, lbl_sf = _demo_labelled_slider(inner, 4, "Signal Freq", 1.0:1.0:20.0, 4.0, v -> "$(round(Int, v)) Hz", ctrl_sz)
    sl_nl, lbl_nl = _demo_labelled_slider(inner, 5, "Noise", 0.0:0.05:1.0, 0.3, v -> "$(round(v, digits=2))", ctrl_sz)

    Label(inner[1, 6], "Kernel type", fontsize = ctrl_sz, halign = :center)
    menu_kt = Menu(inner[2, 6], options = kernel_options, default = "Gaussian", fontsize = ctrl_sz)
    Label(inner[3, 6], "", fontsize = ctrl_sz, halign = :center)

    function update_kw_label()
        lbl_kw.text =
            kernel_type[] === :wavelet ? "→ eff. $(round(Int, wavelet_cycles[] / wavelet_freq[] * 1000)) ms" :
            "$(round(Int, kernel_width[] * 1000)) ms"
    end

    on(sl_pos.value) do v
        kernel_pos[] = v
        lbl_pos.text = "pos: $(round(v * EPOCH, digits=2)) s"
    end
    on(sl_kw.value) do v
        kernel_width[] = v
        update_kw_label()
    end
    on(sl_cy.value) do v
        wavelet_cycles[] = v
        lbl_cy.text = "$(round(Int, v))"
        update_kw_label()
    end
    on(sl_wf.value) do v
        wavelet_freq[] = v
        lbl_wf.text = "$(round(Int, v)) Hz"
        update_kw_label()
    end
    on(sl_sf.value) do v
        sig_freq[] = v
        lbl_sf.text = "$(round(Int, v)) Hz"
    end
    on(sl_nl.value) do v
        noise_lvl[] = v
        lbl_nl.text = "$(round(v, digits=2))"
    end
    on(menu_kt.selection) do sel
        kernel_type[] = sel == "Gaussian" ? :gaussian : sel == "Boxcar" ? :boxcar : :wavelet
        update_kw_label()
    end

    for col = 1:4
        colsize!(ctrl, col, Relative(1 / 4))
    end
    for col = 1:6
        colsize!(inner, col, Relative(1 / 6))
    end

    preset_grid = GridLayout(fig[4, 1], tellwidth = false, colgap = 10, padding = (10, 10, 10, 10))
    Label(preset_grid[1, 1], "Presets", fontsize = ctrl_sz, font = :bold, halign = :right)
    btn_play = Button(preset_grid[1, 2], label = "▶ Play Animation", width = 180)
    btn_denoise = Button(preset_grid[1, 3], label = "Smooth Denoise (Gaussian)", width = 240)
    btn_wavelet = Button(preset_grid[1, 4], label = "Amplitude Envelope (Wavelet)", width = 260)

    is_playing = Ref(false)
    on(btn_play.clicks) do _
        is_playing[] && return
        is_playing[] = true
        @async begin
            try
                for p = 0.0:0.01:1.0
                    (!is_playing[] || !events(fig).window_open[]) && break
                    set_close_to!(sl_pos, p)
                    sleep(0.05)
                end
            finally
                is_playing[] = false
            end
        end
    end

    on(btn_denoise.clicks) do _
        is_playing[] = false
        is_applying_preset[] = true
        try
            menu_kt.i_selected[] = 1
            set_close_to!(sl_kw, 0.15)
            set_close_to!(sl_sf, 4.0)
            set_close_to!(sl_nl, 0.5)
            set_close_to!(sl_pos, 0.0)
        finally
            is_applying_preset[] = false
        end
        update()
    end

    on(btn_wavelet.clicks) do _
        is_playing[] = false
        is_applying_preset[] = true
        try
            menu_kt.i_selected[] = 3
            set_close_to!(sl_wf, 8.0)
            set_close_to!(sl_sf, 8.0)
            set_close_to!(sl_nl, 0.0)
            set_close_to!(sl_pos, 0.0)
        finally
            is_applying_preset[] = false
        end
        update()
    end

    rowsize!(fig.layout, 1, Relative(0.38))
    rowsize!(fig.layout, 2, Relative(0.38))
    rowsize!(fig.layout, 3, Relative(0.18))
    rowsize!(fig.layout, 4, Relative(0.06))

    _display_demo_figure(fig)
    return fig, (ax_sig, ax_out)
end
