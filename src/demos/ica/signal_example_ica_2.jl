"""
    signal_example_ica_2()

Interactive Demo — Signal Mixing & ICA Unmixing.

Demonstrates the fundamental concept behind ICA: EEG electrodes record
*mixtures* of underlying sources, and ICA can separate them back out.

Uses a simple **2-source / 2-electrode** scenario:
- **Source 1**: A brain oscillation (alpha bursts at a controllable frequency)
- **Source 2**: Eye blink artifacts (sharp upward deflections)

A single "Mix Amount" slider controls how much the two sources bleed into
each electrode. Users can attempt manual unmixing before using the **Unmix!**
button, which runs the Infomax ICA algorithm to automatically recover the
original sources from the mixtures.

## Layout

| Column 1 (wide) | Column 2 (wide) | Column 3 (wide) |
|------------------|-------------------|-------------------|
| Source 1: Brain  | Electrode 1 (mix) | Recovered 1       |
| Source 2: Blink  | Electrode 2 (mix) | Recovered 2       |

Below the plots: mixing matrix, unmixing matrix, and match scores.

## Controls

| Control          | Range     | Description                                |
|:-----------------|:----------|:-------------------------------------------|
| Mix Amount       | 0 – 1     | How much sources bleed across electrodes   |
| Brain Freq       | 4 – 20 Hz | Frequency of the oscillatory source        |
| Blink Size       | 0.5 – 5   | Amplitude of blink artifacts               |
| Noise            | 0 – 1     | Additive sensor noise                      |
| Subtract E2→E1   | 0 – 1     | Manual: remove E2 contribution from E1     |
| Subtract E1→E2   | 0 – 1     | Manual: remove E1 contribution from E2     |
| Unmix!           | button    | Run Infomax ICA to find optimal weights    |
| Reset            | button    | Clear ICA result                           |

## Teaching Notes

Think of it like two microphones in a room with two people talking.
Each microphone picks up *both* voices, just at different volumes.
ICA figures out who said what — even though neither microphone
recorded a single voice cleanly.

In EEG terms:
- The "people talking" are your brain signal and your eye blinks
- The "microphones" are your electrodes
- ICA is the clever algorithm that untangles them

## See Also

- [`signal_example_ica_1`](@ref) — simpler: "What is a mixture?" with just signal + weight
- [`signal_example_ica_3`](@ref) — more advanced 3-source version with rotation geometry and scatter plots
- [`signal_example_ica_4`](@ref) — sphering: morphing data into perfect spheres
- [`signal_example_ica_5`](@ref) — inside the black box: optimization landscape and gradient ascent

# Examples
```julia
using EegFun
EegFun.signal_example_ica_2()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica_2()
    # ── 1. Setup & Constants ──────────────────────────────────────────────────
    FS    = 256.0
    T     = 4.0
    t_arr = collect(range(0.0, T - 1.0 / FS, step = 1.0 / FS))
    n     = length(t_arr)

    fig = Figure(size = (1800, 900), figure_padding = (10, 10, 10, 10))

    # ── 2. Observables (State & UI) ───────────────────────────────────────────
    lbl_sz_obs   = Observable(14)
    tick_sz_obs  = Observable(12)
    title_sz_obs = Observable(16)
    ctrl_sz_obs  = Observable(13)
    mat_sz_obs   = Observable(18)

    on(fig.scene.viewport) do area
        sf             = area.widths[1] / 1800
        lbl_sz_obs[]   = max(10, round(Int, 14 * sf))
        tick_sz_obs[]  = max(8, round(Int, 12 * sf))
        title_sz_obs[] = max(11, round(Int, 15 * sf))
        ctrl_sz_obs[]  = max(10, round(Int, 13 * sf))
        mat_sz_obs[]   = max(9, round(Int, 18 * sf))
    end

    # Parameters
    mix_amount_obs = Observable(0.5)
    brain_freq_obs = Observable(10.0)
    blink_size_obs = Observable(3.0)
    noise_lvl_obs  = Observable(0.05)
    manual_w1_obs  = Observable(0.0)
    manual_w2_obs  = Observable(0.0)
    use_ica_obs    = Observable(false)

    # Data Buffers
    s1_buf, s2_buf = zeros(n), zeros(n)
    e1_buf, e2_buf = zeros(n), zeros(n)
    r1_buf, r2_buf = zeros(n), zeros(n)
    noise1_buf, noise2_buf = randn(n), randn(n)  # stable noise — only re-drawn when Noise slider moves

    src1_obs, src2_obs = Observable(s1_buf), Observable(s2_buf)
    mix1_obs, mix2_obs = Observable(e1_buf), Observable(e2_buf)
    rec1_obs, rec2_obs = Observable(r1_buf), Observable(r2_buf)

    # UI Feedback
    e1_title_obs  = Observable("E1")
    e2_title_obs  = Observable("E2")
    rec_title_obs = Observable("Recovered")

    mix_mat_obs     = Observable([1.0 0.5; 0.5 1.0])
    inv_mat_obs     = Observable([1.0 0.0; 0.0 1.0])
    unmix_mat_obs   = Observable([1.0 0.0; 0.0 1.0])
    inv_title_obs   = Observable("Exact Inverse M⁻¹")
    unmix_title_obs = Observable("Your Weights")

    # ICA state
    ica_c1_buf, ica_c2_buf = zeros(n), zeros(n)
    ica_weights_buf        = Matrix{Float64}(I, 2, 2)

    # ── 3. Internal Helper Functions ──────────────────────────────────────────

    function _make_sources!(s1, s2, freq, blink_amp, noise)
        fill!(s1, 0.0)
        fill!(s2, 0.0)
        for center in [0.5, 1.4, 2.3, 3.2]
            σ = 0.12
            @. s1 += 1.5 * sin(2π * freq * (t_arr - center)) * exp(-((t_arr - center)^2) / (2 * σ^2))
        end
        for onset_time in [0.9, 1.8, 2.8]
            onset = round(Int, onset_time * FS)
            width = round(Int, 0.25 * FS)
            for k = 1:width
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

    function _apply_mixing!(e1, e2, s1, s2, m)
        @. e1 = s1 + m * s2
        @. e2 = m * s1 + s2
    end

    function _make_matrix_equation!(parent, mat_obs, title_obs, left_vec, right_vec)
        # Clear title for better spacing
        Label(parent[0, 1:7], title_obs; fontsize = title_sz_obs, font = :bold, valign = :bottom, padding = (0, 0, 5, 0))

        ax = Axis(parent[1, 2]; aspect = DataAspect(), yreversed = true, width = 90, height = 90)
        hidedecorations!(ax) # Hide ticks and labels

        heatmap!(ax, @lift(abs.($mat_obs)); colormap = [:white, RGBf(0.2, 0.4, 0.7)], colorrange = (0.0, 2.0))

        for r = 1:2, c = 1:2
            t_obj = lift(mat_obs) do m
                val = m[r, c]
                if isnan(val)
                    "∅"
                elseif abs(val) >= 10.0
                    string(round(Int, val))
                else
                    rv = round(val; digits = 2)
                    isinteger(rv) ? string(Int(rv)) : string(rv)
                end
            end
            c_obj = lift(mat_obs) do m
                abs(m[r, c]) > 1.0 ? :white : :black
            end
            text!(ax, c, r; text = t_obj, align = (:center, :center), fontsize = @lift($mat_sz_obs * 0.7), color = c_obj, font = :bold)
        end

        # Symbols & Brackets
        Label(parent[1, 1], "["; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)
        Label(parent[1, 3], "]"; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)
        Label(parent[1, 4], "×"; fontsize = @lift($mat_sz_obs * 1.5), font = :bold)

        # Left Vector
        v_l = GridLayout(parent[1, 5], tellheight = false)
        Label(v_l[1, 1], "["; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)
        Label(v_l[1, 2], left_vec; fontsize = mat_sz_obs, font = :bold)
        Label(v_l[1, 3], "]"; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)

        Label(parent[1, 6], "="; fontsize = @lift($mat_sz_obs * 1.5), font = :bold)

        # Right Vector
        v_r = GridLayout(parent[1, 7], tellheight = false)
        Label(v_r[1, 1], "["; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)
        Label(v_r[1, 2], right_vec; fontsize = mat_sz_obs, font = :bold)
        Label(v_r[1, 3], "]"; fontsize = @lift($mat_sz_obs * 3.5), font = :regular, color = :gray80)

        colgap!(parent, 3)
    end

    function _labelled_slider!(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header; fontsize = ctrl_sz_obs, halign = :center, valign = :bottom)
        sl  = Slider(parent[2, col]; range = range_vals, startvalue = startval, height = 10)
        lbl = Label(parent[3, col], fmt(startval); fontsize = ctrl_sz_obs, halign = :center, valign = :top)
        return sl, lbl
    end

    function update()
        _make_sources!(s1_buf, s2_buf, brain_freq_obs[], blink_size_obs[], noise_lvl_obs[])
        m = mix_amount_obs[]
        _apply_mixing!(e1_buf, e2_buf, s1_buf, s2_buf, m)

        mstr           = round(m; digits = 2)
        e1_title_obs[] = "E1 = 1.00 × Brain + $(mstr) × Blink"
        e2_title_obs[] = "E2 = $(mstr) × Brain + 1.00 × Blink"
        mix_mat_obs[]  = [1.0 m; m 1.0]

        det_val = 1.0 - m^2
        is_singular = abs(det_val) < 1e-4

        if is_singular
            inv_mat_obs[]   = fill(NaN, 2, 2)
            inv_title_obs[] = "Singular Matrix\n(M⁻¹ Does Not Exist)"
        else
            inv_mat_obs[]   = [1.0/det_val -m/det_val; -m/det_val 1.0/det_val]
            inv_title_obs[] = "Exact Inverse M⁻¹"
        end

        if use_ica_obs[]
            copyto!(r1_buf, ica_c1_buf)
            copyto!(r2_buf, ica_c2_buf)
            unmix_mat_obs[]   = ica_weights_buf
            unmix_title_obs[] = is_singular ? "Degenerate Case!\n(Sensors are Identical!)" : "ICA Found Matrix"
            rec_title_obs[]   = "Recovered by ICA  ✓"
        else
            w1, w2            = manual_w1_obs[], manual_w2_obs[]
            @. r1_buf         = e1_buf - w1 * e2_buf
            @. r2_buf         = e2_buf - w2 * e1_buf
            unmix_mat_obs[]   = [1.0 -round(w1; digits = 2); -round(w2; digits = 2) 1.0]
            unmix_title_obs[] = "Your Weights"
            rec_title_obs[]   = (w1 == 0 && w2 == 0) ? "Recovered  ·  try the sliders or click Unmix!" : "Recovered (manual unmixing)"
        end
        notify(src1_obs)
        notify(src2_obs)
        notify(mix1_obs)
        notify(mix2_obs)
        notify(rec1_obs)
        notify(rec2_obs)

    end

    # ── 4. UI Layout ──────────────────────────────────────────
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
    ylims!(ax_s1, y_lims)
    ax_s2 = Axis(fig[2, 1]; title = "Source 2: Eye Blink Artifact", ax_args...)
    ylims!(ax_s2, y_lims)
    ax_e1 = Axis(fig[1, 2]; title = e1_title_obs, ax_args...)
    ylims!(ax_e1, y_lims)
    ax_e2 = Axis(fig[2, 2]; title = e2_title_obs, ax_args...)
    ylims!(ax_e2, y_lims)
    ax_r1 = Axis(fig[1, 3]; title = rec_title_obs, ax_args...)
    ylims!(ax_r1, y_lims)
    ax_r2 = Axis(fig[2, 3]; title = "", ax_args...)
    ylims!(ax_r2, y_lims)
    linkxaxes!(ax_s1, ax_s2, ax_e1, ax_e2, ax_r1, ax_r2)

    lines!(ax_s1, t_arr, src1_obs; color = :black, linewidth = 2.0)
    lines!(ax_s2, t_arr, src2_obs; color = :black, linewidth = 2.0)
    lines!(ax_e1, t_arr, mix1_obs; color = :black, linewidth = 1.5)
    lines!(ax_e2, t_arr, mix2_obs; color = :black, linewidth = 1.5)
    lines!(ax_r1, t_arr, rec1_obs; color = :black, linewidth = 2.0)
    lines!(ax_r2, t_arr, rec2_obs; color = :black, linewidth = 2.0)

    Label(fig[0, 1], "What's really happening", fontsize = @lift($title_sz_obs + 2), font = :bold)
    Label(fig[0, 2], "What electrodes record", fontsize = @lift($title_sz_obs + 2), font = :bold)
    Label(fig[0, 3], "What ICA recovers", fontsize = @lift($title_sz_obs + 2), font = :bold)

    _make_matrix_equation!(GridLayout(fig[3, 1]), mix_mat_obs, Observable("Mixing M"), "S1\nS2", "E1\nE2")
    _make_matrix_equation!(GridLayout(fig[3, 2]), inv_mat_obs, inv_title_obs, "E1\nE2", "S1\nS2")
    _make_matrix_equation!(GridLayout(fig[3, 3]), unmix_mat_obs, unmix_title_obs, "E1\nE2", "R1\nR2")

    ctrl = GridLayout(fig[4, 1:3]; colgap = 120)
    rowgap!(ctrl, 0)
    sl_m, lbl_m = _labelled_slider!(ctrl, 1, "Mix Amount", 0.0:0.05:1.0, 0.5, v -> "$(round(v; digits=2))")
    sl_f, lbl_f = _labelled_slider!(ctrl, 2, "Brain Freq", 4.0:1.0:20.0, 10.0, v -> "$(Int(v)) Hz")
    sl_b, lbl_b = _labelled_slider!(ctrl, 3, "Blink Size", 0.5:0.5:5.0, 3.0, v -> "$(round(v; digits=1))")
    sl_n, lbl_n = _labelled_slider!(ctrl, 4, "Noise", 0.0:0.05:1.0, 0.05, v -> "$(round(v; digits=2))")
    sl_w1, lbl_w1 = _labelled_slider!(ctrl, 5, "Subtract E2→E1", 0.0:0.05:1.5, 0.0, v -> "$(round(v; digits=2))")
    sl_w2, lbl_w2 = _labelled_slider!(ctrl, 6, "Subtract E1→E2", 0.0:0.05:1.5, 0.0, v -> "$(round(v; digits=2))")

    rowsize!(ctrl, 1, Fixed(28))
    rowsize!(ctrl, 2, Fixed(12))
    rowsize!(ctrl, 3, Fixed(20))
    btn_unmix = Button(ctrl[2, 7]; label = "Unmix!", fontsize = @lift($ctrl_sz_obs + 2))
    btn_reset = Button(ctrl[2, 8]; label = "Reset", fontsize = ctrl_sz_obs)

    # Global Layout
    rowgap!(fig.layout, 1, 15)
    rowgap!(fig.layout, 3, 20)
    rowgap!(fig.layout, 4, 30)
    rowsize!(fig.layout, 0, Relative(0.06))
    rowsize!(fig.layout, 1, Relative(0.24))
    rowsize!(fig.layout, 2, Relative(0.24))
    rowsize!(fig.layout, 3, Relative(0.25))
    rowsize!(fig.layout, 4, Relative(0.18))
    for c = 1:3
        colsize!(fig.layout, c, Relative(0.33))
    end

    # ── 5. Wiring ─────────────────────────────────────────────────────────────
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
    on(btn_unmix.clicks) do _
        dat = vcat(mix1_obs[]', mix2_obs[]')
        res = infomax_ica(
            dat,
            Layout(DataFrame(label = [:E1, :E2], inc = [0, 0], azi = [0, 0]), nothing, nothing, nothing),
            "demo";
            n_components = 2,
        )
        rec = res.unmixing * ((dat .- res.mean) ./ res.scale)
        s1r, s2r = src1_obs[], src2_obs[]
        cor1, cor2 = abs(cor(vec(rec[1, :]), s1r)), abs(cor(vec(rec[2, :]), s1r))
        swap = cor2 > cor1
        c1, c2 = swap ? (rec[2, :], rec[1, :]) : (rec[1, :], rec[2, :])
        flip1 = cor(c1, s1r) < 0
        flip2 = cor(c2, s2r) < 0
        if flip1
            ;
            c1 .*= -1;
        end
        if flip2
            ;
            c2 .*= -1;
        end
        copyto!(ica_c1_buf, c1 .- mean(c1))
        copyto!(ica_c2_buf, c2 .- mean(c2))
        W = res.unmixing[(swap ? [2, 1] : [1, 2]), :]
        if flip1
            ;
            W[1, :] .*= -1;
        end
        if flip2
            ;
            W[2, :] .*= -1;
        end
        ica_weights_buf .= W
        w1e, w2e = clamp(-W[1, 2] / W[1, 1], 0, 1.5), clamp(-W[2, 1] / W[2, 2], 0, 1.5)
        use_ica_obs[] = true   # set before set_close_to! to avoid intermediate manual-mode update()
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
    display(fig)
    return fig
end
