"""
    signal_example_ica_separation()

Interactive 2-source / 2-electrode ICA separation demo illustrating forward mixing matrices,
manual subtraction weights, and automatic recovery via Infomax ICA.

## See Also
- [`signal_example_ica_mixture`](@ref), [`signal_example_ica_geometry`](@ref),
  [`signal_example_ica_clt`](@ref), [`signal_example_ica_sphering`](@ref),
  [`signal_example_ica_optimization`](@ref)

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_separation()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica_separation()
    _check_interactive_demo_backend("ICA Separation demo") || return nothing

    FS = 256.0
    T = 4.0
    t_arr = collect(range(0.0, T - 1.0 / FS, step = 1.0 / FS))
    n = length(t_arr)

    fig = _eegfun_figure(size = (1400, 850), figure_padding = (10, 10, 10, 10))
    lbl_sz_obs, tick_sz_obs, title_sz_obs, ctrl_sz_obs, mat_sz_obs =
        _demo_adaptive_fonts(fig, 1400, (14, 10), (12, 8), (15, 11), (13, 10), (18, 9))

    mix_amount_obs = Observable(0.5)
    brain_freq_obs = Observable(10.0)
    blink_size_obs = Observable(3.0)
    noise_lvl_obs = Observable(0.05)
    manual_w1_obs = Observable(0.0)
    manual_w2_obs = Observable(0.0)
    use_ica_obs = Observable(false)

    s1_buf, s2_buf = zeros(n), zeros(n)
    e1_buf, e2_buf = zeros(n), zeros(n)
    r1_buf, r2_buf = zeros(n), zeros(n)
    noise1_buf, noise2_buf = randn(n), randn(n)

    src1_obs, src2_obs = Observable(s1_buf), Observable(s2_buf)
    mix1_obs, mix2_obs = Observable(e1_buf), Observable(e2_buf)
    rec1_obs, rec2_obs = Observable(r1_buf), Observable(r2_buf)

    e1_title_obs = Observable("E1")
    e2_title_obs = Observable("E2")
    rec_title_obs = Observable("Recovered")

    mix_mat_obs = Observable([1.0 0.5; 0.5 1.0])
    inv_mat_obs = Observable([1.0 0.0; 0.0 1.0])
    unmix_mat_obs = Observable([1.0 0.0; 0.0 1.0])
    inv_title_obs = Observable("Exact Inverse M⁻¹")
    unmix_title_obs = Observable("Your Weights")

    ica_c1_buf, ica_c2_buf = zeros(n), zeros(n)
    ica_weights_buf = Matrix{Float64}(I, 2, 2)

    function _make_sources!(s1, s2, freq, blink_amp, noise)
        fill!(s1, 0.0)
        fill!(s2, 0.0)
        σ = 0.12
        inv_2σ2 = 1.0 / (2.0 * σ^2)
        for center in (0.5, 1.4, 2.3, 3.2)
            @. s1 += 1.5 * sin(2π * freq * (t_arr - center)) * exp(-((t_arr - center)^2) * inv_2σ2)
        end
        width = round(Int, 0.25 * FS)
        for onset_time in (0.9, 1.8, 2.8)
            onset = round(Int, onset_time * FS)
            @inbounds for k = 1:width
                idx = onset + k
                if 1 <= idx <= n
                    s2[idx] += blink_amp * (k / 20.0) * exp(-(k - 20.0) / 30.0)
                end
            end
        end
        if noise > 0.0
            @. s1 += noise * 0.3 * noise1_buf
            @. s2 += noise * 0.3 * noise2_buf
        end
    end

    function _make_matrix_equation!(parent, mat_obs, title_obs, left_vec, right_vec)
        Label(parent[0, 1:7], title_obs; fontsize = title_sz_obs, font = :bold, valign = :bottom, padding = (0, 0, 5, 0))

        ax = Axis(parent[1, 2]; aspect = DataAspect(), yreversed = true, width = 90, height = 90)
        hidedecorations!(ax)
        heatmap!(ax, @lift(abs.($mat_obs)); colormap = [:white, RGBf(0.2, 0.4, 0.7)], colorrange = (0.0, 2.0))

        for r = 1:2, c = 1:2
            t_obj = lift(mat_obs) do m
                val = m[r, c]
                isnan(val) && return "∅"
                abs(val) >= 10.0 && return string(round(Int, val))
                rv = round(val; digits = 2)
                return isinteger(rv) ? string(Int(rv)) : string(rv)
            end
            c_obj = lift(m -> abs(m[r, c]) > 1.0 ? :white : :black, mat_obs)
            text!(ax, c, r; text = t_obj, align = (:center, :center), fontsize = @lift($mat_sz_obs * 0.7), color = c_obj, font = :bold)
        end

        bracket_sz = @lift($mat_sz_obs * 3.5)
        op_sz = @lift($mat_sz_obs * 1.5)
        Label(parent[1, 1], "["; fontsize = bracket_sz, font = :regular, color = :gray80)
        Label(parent[1, 3], "]"; fontsize = bracket_sz, font = :regular, color = :gray80)
        Label(parent[1, 4], "×"; fontsize = op_sz, font = :bold)

        v_l = GridLayout(parent[1, 5], tellheight = false)
        Label(v_l[1, 1], "["; fontsize = bracket_sz, font = :regular, color = :gray80)
        Label(v_l[1, 2], left_vec; fontsize = mat_sz_obs, font = :bold)
        Label(v_l[1, 3], "]"; fontsize = bracket_sz, font = :regular, color = :gray80)

        Label(parent[1, 6], "="; fontsize = op_sz, font = :bold)

        v_r = GridLayout(parent[1, 7], tellheight = false)
        Label(v_r[1, 1], "["; fontsize = bracket_sz, font = :regular, color = :gray80)
        Label(v_r[1, 2], right_vec; fontsize = mat_sz_obs, font = :bold)
        Label(v_r[1, 3], "]"; fontsize = bracket_sz, font = :regular, color = :gray80)

        colgap!(parent, 3)
    end

    function update()
        _make_sources!(s1_buf, s2_buf, brain_freq_obs[], blink_size_obs[], noise_lvl_obs[])
        m = mix_amount_obs[]
        @. e1_buf = s1_buf + m * s2_buf
        @. e2_buf = m * s1_buf + s2_buf

        mstr = round(m; digits = 2)
        e1_title_obs[] = "E1 = 1.00 × Brain + $(mstr) × Blink"
        e2_title_obs[] = "E2 = $(mstr) × Brain + 1.00 × Blink"
        mix_mat_obs[] = [1.0 m; m 1.0]

        det_val = 1.0 - m^2
        is_singular = abs(det_val) < 1e-4
        if is_singular
            inv_mat_obs[] = fill(NaN, 2, 2)
            inv_title_obs[] = "Singular Matrix\n(M⁻¹ Does Not Exist)"
        else
            inv_mat_obs[] = [1.0/det_val -m/det_val; -m/det_val 1.0/det_val]
            inv_title_obs[] = "Exact Inverse M⁻¹"
        end

        if use_ica_obs[]
            copyto!(r1_buf, ica_c1_buf)
            copyto!(r2_buf, ica_c2_buf)
            unmix_mat_obs[] = ica_weights_buf
            unmix_title_obs[] = is_singular ? "Degenerate Case!\n(Sensors are Identical!)" : "ICA Found Matrix"
            rec_title_obs[] = "Recovered by ICA  ✓"
        else
            w1, w2 = manual_w1_obs[], manual_w2_obs[]
            @. r1_buf = e1_buf - w1 * e2_buf
            @. r2_buf = e2_buf - w2 * e1_buf
            unmix_mat_obs[] = [1.0 -round(w1; digits = 2); -round(w2; digits = 2) 1.0]
            unmix_title_obs[] = "Your Weights"
            rec_title_obs[] = (w1 == 0 && w2 == 0) ? "Recovered  ·  try the sliders or click Unmix!" : "Recovered (manual unmixing)"
        end

        for obs in (src1_obs, src2_obs, mix1_obs, mix2_obs, rec1_obs, rec2_obs)
            notify(obs)
        end
    end

    ax_args = (
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        titlesize = title_sz_obs,
        xlabelsize = lbl_sz_obs,
        ylabelsize = lbl_sz_obs,
        xticklabelsize = tick_sz_obs,
        yticklabelsize = tick_sz_obs,
    )
    y_lims = (-7, 7)

    ax_s1 = Axis(fig[1, 1]; title = "Source 1: Brain Oscillation (alpha)", ax_args...)
    ax_s2 = Axis(fig[2, 1]; title = "Source 2: Eye Blink Artifact", ax_args...)
    ax_e1 = Axis(fig[1, 2]; title = e1_title_obs, ax_args...)
    ax_e2 = Axis(fig[2, 2]; title = e2_title_obs, ax_args...)
    ax_r1 = Axis(fig[1, 3]; title = rec_title_obs, ax_args...)
    ax_r2 = Axis(fig[2, 3]; title = "", ax_args...)
    for ax in (ax_s1, ax_s2, ax_e1, ax_e2, ax_r1, ax_r2)
        ylims!(ax, y_lims)
    end
    linkxaxes!(ax_s1, ax_s2, ax_e1, ax_e2, ax_r1, ax_r2)

    lines!(ax_s1, t_arr, src1_obs; color = :black, linewidth = 2.0)
    lines!(ax_s2, t_arr, src2_obs; color = :black, linewidth = 2.0)
    lines!(ax_e1, t_arr, mix1_obs; color = :black, linewidth = 1.5)
    lines!(ax_e2, t_arr, mix2_obs; color = :black, linewidth = 1.5)
    lines!(ax_r1, t_arr, rec1_obs; color = :black, linewidth = 2.0)
    lines!(ax_r2, t_arr, rec2_obs; color = :black, linewidth = 2.0)

    hdr_sz = @lift($title_sz_obs + 2)
    Label(fig[0, 1], "What's really happening", fontsize = hdr_sz, font = :bold)
    Label(fig[0, 2], "What electrodes record", fontsize = hdr_sz, font = :bold)
    Label(fig[0, 3], "What ICA recovers", fontsize = hdr_sz, font = :bold)

    _make_matrix_equation!(GridLayout(fig[3, 1]), mix_mat_obs, Observable("Mixing M"), "S1\nS2", "E1\nE2")
    _make_matrix_equation!(GridLayout(fig[3, 2]), inv_mat_obs, inv_title_obs, "E1\nE2", "S1\nS2")
    _make_matrix_equation!(GridLayout(fig[3, 3]), unmix_mat_obs, unmix_title_obs, "E1\nE2", "R1\nR2")

    ctrl = GridLayout(fig[4, 1:3]; colgap = 120)
    rowgap!(ctrl, 0)
    sl_m, lbl_m = _demo_labelled_slider(ctrl, 1, "Mix Amount", 0.0:0.05:1.0, 0.5, v -> "$(round(v; digits=2))", ctrl_sz_obs; height = 10)
    sl_f, lbl_f = _demo_labelled_slider(ctrl, 2, "Brain Freq", 4.0:1.0:20.0, 10.0, v -> "$(Int(v)) Hz", ctrl_sz_obs; height = 10)
    sl_b, lbl_b = _demo_labelled_slider(ctrl, 3, "Blink Size", 0.5:0.5:5.0, 3.0, v -> "$(round(v; digits=1))", ctrl_sz_obs; height = 10)
    sl_n, lbl_n = _demo_labelled_slider(ctrl, 4, "Noise", 0.0:0.05:1.0, 0.05, v -> "$(round(v; digits=2))", ctrl_sz_obs; height = 10)
    sl_w1, lbl_w1 =
        _demo_labelled_slider(ctrl, 5, "Subtract E2→E1", 0.0:0.05:1.5, 0.0, v -> "$(round(v; digits=2))", ctrl_sz_obs; height = 10)
    sl_w2, lbl_w2 =
        _demo_labelled_slider(ctrl, 6, "Subtract E1→E2", 0.0:0.05:1.5, 0.0, v -> "$(round(v; digits=2))", ctrl_sz_obs; height = 10)

    rowsize!(ctrl, 1, Fixed(28))
    rowsize!(ctrl, 2, Fixed(12))
    rowsize!(ctrl, 3, Fixed(20))
    btn_unmix = Button(ctrl[2, 7]; label = "Unmix!", fontsize = @lift($ctrl_sz_obs + 2))
    btn_reset = Button(ctrl[2, 8]; label = "Reset", fontsize = ctrl_sz_obs)

    rowgap!(fig.layout, 0, 15)
    rowgap!(fig.layout, 2, 20)
    rowgap!(fig.layout, 3, 30)
    rowsize!(fig.layout, 0, Relative(0.06))
    rowsize!(fig.layout, 1, Relative(0.24))
    rowsize!(fig.layout, 2, Relative(0.24))
    rowsize!(fig.layout, 3, Relative(0.25))
    rowsize!(fig.layout, 4, Relative(0.18))
    for c = 1:3
        colsize!(fig.layout, c, Relative(0.33))
    end

    on(sl_m.value) do v
        mix_amount_obs[] = v
        lbl_m.text = "$(round(v; digits=2))"
    end
    on(sl_f.value) do v
        brain_freq_obs[] = v
        lbl_f.text = "$(Int(v)) Hz"
    end
    on(sl_b.value) do v
        blink_size_obs[] = v
        lbl_b.text = "$(round(v; digits=1))"
    end
    on(sl_n.value) do v
        noise_lvl_obs[] = v
        lbl_n.text = "$(round(v; digits=2))"
        randn!(noise1_buf)
        randn!(noise2_buf)
    end
    on(sl_w1.value) do v
        manual_w1_obs[] = v
        lbl_w1.text = "$(round(v; digits=2))"
    end
    on(sl_w2.value) do v
        manual_w2_obs[] = v
        lbl_w2.text = "$(round(v; digits=2))"
    end

    for b in (mix_amount_obs, brain_freq_obs, blink_size_obs, noise_lvl_obs, manual_w1_obs, manual_w2_obs)
        on(b) do _
            use_ica_obs[] = false
            update()
        end
    end

    dummy_layout = Layout(DataFrame(label = [:E1, :E2], inc = [0, 0], azi = [0, 0]), nothing, nothing, nothing)
    on(btn_unmix.clicks) do _
        dat = vcat(mix1_obs[]', mix2_obs[]')
        res = infomax_ica(dat, dummy_layout, "demo"; n_components = 2)
        rec = res.unmixing * ((dat .- res.mean) ./ res.scale)
        s1r, s2r = src1_obs[], src2_obs[]
        swap = abs(cor(vec(rec[2, :]), s1r)) > abs(cor(vec(rec[1, :]), s1r))
        c1, c2 = swap ? (rec[2, :], rec[1, :]) : (rec[1, :], rec[2, :])
        flip1 = cor(c1, s1r) < 0
        flip2 = cor(c2, s2r) < 0
        flip1 && (c1 .*= -1)
        flip2 && (c2 .*= -1)
        copyto!(ica_c1_buf, c1 .- mean(c1))
        copyto!(ica_c2_buf, c2 .- mean(c2))
        W = res.unmixing[(swap ? [2, 1] : [1, 2]), :]
        flip1 && (W[1, :] .*= -1)
        flip2 && (W[2, :] .*= -1)
        ica_weights_buf .= W
        w1e = clamp(-W[1, 2] / W[1, 1], 0.0, 1.5)
        w2e = clamp(-W[2, 1] / W[2, 2], 0.0, 1.5)
        use_ica_obs[] = true
        set_close_to!(sl_w1, w1e)
        set_close_to!(sl_w2, w2e)
        update()
    end

    on(btn_reset.clicks) do _
        use_ica_obs[] = false
        set_close_to!(sl_w1, 0.0)
        set_close_to!(sl_w2, 0.0)
        update()
    end

    update()
    _display_demo_figure(fig)
    return fig
end
