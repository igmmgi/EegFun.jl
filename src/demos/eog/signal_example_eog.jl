"""
    signal_example_heog(; kwargs...)

Interactive Horizontal Electrooculography (hEOG) simulator illustrating corneo-retinal standing
dipoles, volume-conducted electric fields, common-mode rejection (CMR), and bipolar outer-canthi
derivation (\$V_{\\text{out}} = V_{(+)} - V_{(-)}\$).

# Example
```julia
using GLMakie, EegFun
signal_example_heog()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax_eyes::Axis`: The anatomical dipole & circuit schematic axis
- `ax_eog::Axis`: The potential vs. time signal recording axis
"""
function signal_example_heog(; kwargs...)
    _check_interactive_demo_backend("Interactive hEOG demo") || return nothing
    _set_window_title("Interactive Bipolar hEOG Recording Simulator")

    fig = _eegfun_figure(size = (1350, 850), title = "Interactive Bipolar hEOG Recording Simulator", backgroundcolor = :white)
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1350, (18, 14), (14, 11), (12, 10), (13, 11))

    angle_deg = Observable(30.0)
    current_display_angle = Observable(30.0)
    duration_ms = Observable(700.0)
    noise_amp = Observable(2.0)
    show_fields = Observable(true)
    playback_time = Observable(-1.0)

    CL = Point2f(-2.2, 0.5)
    CR = Point2f(2.2, 0.5)
    R = 1.15f0

    unit_circ = [Point2f(cos(t), sin(t)) for t in range(0, 2π, length = 64)]

    ax_eyes = Axis(fig[1, 1], title = "Corneo-Retinal Dipoles & Bipolar Amplifier Circuit", titlesize = title_font, aspect = DataAspect())
    hidedecorations!(ax_eyes)
    hidespines!(ax_eyes)
    xlims!(ax_eyes, -5.2, 5.2)
    ylims!(ax_eyes, -3.2, 2.5)

    eog_colormap = cgrad([RGBf(0.05, 0.22, 0.95), RGBf(0.55, 0.10, 0.75), RGBf(0.95, 0.08, 0.20)])
    n_slices = 64
    s_vals = collect(Float32, range(-R, R, length = n_slices + 1))
    w_vals = [sqrt(max(0.0f0, R * R - s * s)) for s in s_vals]
    slice_colors = [eog_colormap[clamp((s_vals[i] + s_vals[i+1] + 2 * R) / (4 * R), 0.0f0, 1.0f0)] for i = 1:n_slices]

    function get_gradient_polys(C::Point2f, deg::Real)
        sθ, cθ = sincos(Float32(deg2rad(deg)))
        u = Point2f(sθ, cθ)
        v = Point2f(cθ, -sθ)
        polys = Vector{Vector{Point2f}}(undef, n_slices)
        @inbounds for i = 1:n_slices
            s1, s2 = s_vals[i], s_vals[i+1]
            w1, w2 = w_vals[i], w_vals[i+1]
            c1 = C + s1 * u
            c2 = C + s2 * u
            polys[i] = Point2f[c1-w1*v, c1+w1*v, c2+w2*v, c2-w2*v]
        end
        return polys
    end

    poly!(ax_eyes, lift(deg -> get_gradient_polys(CL, deg), current_display_angle), color = slice_colors, strokewidth = 0)
    poly!(ax_eyes, lift(deg -> get_gradient_polys(CR, deg), current_display_angle), color = slice_colors, strokewidth = 0)

    for C in (CL, CR)
        poly!(ax_eyes, [C + R * p for p in unit_circ], color = :transparent, strokewidth = 2.0, strokecolor = "#1e293b")
    end

    cornea_angles = collect(range(-45.0, 45.0, length = 25))
    cornea_r_bulge = [Float32(R + 0.22 * cosd(ca)^2) for ca in cornea_angles]

    function get_eye_features(C::Point2f, deg::Real)
        sθ, cθ = sincos(Float32(deg2rad(deg)))
        u = Point2f(sθ, cθ)
        n_ca = length(cornea_angles)
        cornea_pts = Vector{Point2f}(undef, 2 * n_ca)
        @inbounds for i = 1:n_ca
            sa, ca_v = sincos(Float32(deg2rad(deg + cornea_angles[i])))
            cornea_pts[i] = C + cornea_r_bulge[i] * Point2f(sa, ca_v)
            cornea_pts[2*n_ca-i+1] = C + R * Point2f(sa, ca_v)
        end
        return (;
            cornea_pts,
            iris_center = C + 0.82f0 * R * u,
            pupil_center = C + 0.85f0 * R * u,
            dipole_seg = Point2f[C-0.95f0*R*u, C+1.40f0*R*u],
            plus_pos = C + 1.62f0 * R * u,
            minus_pos = C - 1.20f0 * R * u,
        )
    end

    feat_L = lift(deg -> get_eye_features(CL, deg), current_display_angle)
    feat_R = lift(deg -> get_eye_features(CR, deg), current_display_angle)

    for feat in (feat_L, feat_R)
        poly!(ax_eyes, lift(f -> f.cornea_pts, feat), color = RGBAf(0.2, 0.85, 0.98, 0.75), strokecolor = "#0284c7", strokewidth = 2.0)
        poly!(ax_eyes, lift(f -> [f.iris_center + 0.38f0 * p for p in unit_circ], feat), color = "#1e293b")
        poly!(ax_eyes, lift(f -> [f.pupil_center + 0.18f0 * p for p in unit_circ], feat), color = :black)
        linesegments!(ax_eyes, lift(f -> f.dipole_seg, feat), color = "#10b981", linewidth = 2.5)
        text!(
            ax_eyes,
            lift(f -> f.plus_pos[1], feat),
            lift(f -> f.plus_pos[2], feat),
            text = "+",
            align = (:center, :center),
            font = :bold,
            fontsize = 20,
            color = "#e11d48",
        )
        text!(
            ax_eyes,
            lift(f -> f.minus_pos[1], feat),
            lift(f -> f.minus_pos[2], feat),
            text = "–",
            align = (:center, :center),
            font = :bold,
            fontsize = 22,
            color = "#2563eb",
        )
    end

    lines!(ax_eyes, [CL[1], CL[1]], [CL[2], CL[2] + 1.8 * R], color = :gray60, linestyle = :dash, linewidth = 1.2)
    lines!(ax_eyes, [CR[1], CR[1]], [CR[2], CR[2] + 1.8 * R], color = :gray60, linestyle = :dash, linewidth = 1.2)

    arc_pts = lift(current_display_angle) do deg
        r_arc = 1.65f0
        angles = deg == 0 ? (0.0:1.0:0.0) : range(0.0, deg, step = deg >= 0 ? 1.0 : -1.0)
        [CR + Point2f(r_arc * sind(a), r_arc * cosd(a)) for a in angles]
    end
    lines!(ax_eyes, arc_pts, color = "#475569", linewidth = 1.8)

    text!(
        ax_eyes,
        lift(deg -> CR[1] + 1.95 * sin(deg * π / 360), current_display_angle),
        lift(deg -> CR[2] + 1.95 * cos(deg * π / 360), current_display_angle),
        text = lift(deg -> "$(deg >= 0 ? "+" : "")$(round(deg, digits=1))°", current_display_angle),
        align = (:center, :center),
        font = :bold,
        fontsize = 14,
        color = "#1e293b",
    )

    phi_range = collect(range(0.12π, 0.88π, length = 35))
    field_lines = lift(current_display_angle, show_fields) do deg, show
        show || return Point2f[]
        sθ, cθ = sincos(deg2rad(deg))
        pts = Vector{Point2f}(undef, 2 * 3 * 2 * (length(phi_range) + 1))
        idx = 0
        @inbounds for C in (CL, CR), r0 in (1.35, 1.85, 2.45), sgn in (1.0, -1.0)
            for phi in phi_range
                sp, cp = sincos(phi)
                r_val = r0 * sp * sp
                xp = sgn * r_val * sp
                yp = r_val * cp
                idx += 1
                pts[idx] = C + Point2f(xp * cθ + yp * sθ, -xp * sθ + yp * cθ)
            end
            idx += 1
            pts[idx] = Point2f(NaN32, NaN32)
        end
        return pts
    end
    lines!(ax_eyes, field_lines, color = RGBAf(0.2, 0.5, 0.35, 0.45), linestyle = :dash, linewidth = 1.2)

    poly!(ax_eyes, Rect2f(-4.3, 0.2, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, -4.3, 1.0, text = "[-] Left Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")
    poly!(ax_eyes, Rect2f(4.08, 0.2, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, 4.19, 1.0, text = "[+] Right Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")

    poly!(
        ax_eyes,
        [Point2f(-0.7, -1.2), Point2f(-0.7, -2.2), Point2f(0.7, -1.7)],
        color = "#f8fafc",
        strokewidth = 2.0,
        strokecolor = "#1e293b",
    )
    text!(ax_eyes, -0.45, -1.40, text = "+", align = (:center, :center), font = :bold, fontsize = 18, color = "#dc2626")
    text!(ax_eyes, -0.45, -2.00, text = "–", align = (:center, :center), font = :bold, fontsize = 20, color = "#2563eb")
    text!(ax_eyes, 0.0, -0.55, text = "Difference Amplifier", align = (:center, :center), font = :bold, fontsize = 16, color = "#1e293b")

    lines!(ax_eyes, [Point2f(0.7, -1.7), Point2f(1.2, -1.7)], color = "#1e293b", linewidth = 2.0)
    scatter!(ax_eyes, [Point2f(1.25, -1.7)], marker = :circle, markersize = 10, color = :white, strokecolor = "#1e293b", strokewidth = 2.0)

    lines!(
        ax_eyes,
        [Point2f(4.08, 0.5), Point2f(3.5, 0.5), Point2f(3.5, -0.85), Point2f(-1.0, -0.85), Point2f(-1.0, -1.40), Point2f(-0.70, -1.40)],
        color = "#dc2626",
        linewidth = 2.0,
    )
    lines!(
        ax_eyes,
        [Point2f(-4.08, 0.5), Point2f(-3.5, 0.5), Point2f(-3.5, -2.00), Point2f(-0.70, -2.00)],
        color = "#2563eb",
        linewidth = 2.0,
    )

    box_text = lift(current_display_angle) do deg
        vd = 5.0 * deg
        vr, vl = 0.5 * vd, -0.5 * vd
        fmt_s(v) = "$(v >= 0 ? "+" : "")$(round(v, digits=1)) µV"
        "V(+) = $(fmt_s(vr))   |   V(–) = $(fmt_s(vl))   |   Vout = $(fmt_s(vd))"
    end
    poly!(ax_eyes, Rect2f(-3.3, -2.88, 6.6, 0.56), color = "#f8fafc", strokewidth = 1.2, strokecolor = "#cbd5e1")
    text!(ax_eyes, 0.0, -2.60, text = box_text, align = (:center, :center), font = :bold, fontsize = 15.5, color = "#0f172a")

    ctrl_layout = fig[2, 1] = GridLayout(tellwidth = false, halign = :center)
    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(slider_sub[1, 1], text = "Gaze:", font = :bold, fontsize = label_font)
    sl_angle = Slider(slider_sub[1, 2], range = -45.0:1.0:45.0, startvalue = 30.0, width = 75)
    Label(slider_sub[1, 3], text = lift(v -> "$(round(v, digits=1))°", sl_angle.value), width = 58, fontsize = label_font)

    Label(slider_sub[1, 5], text = "Duration:", font = :bold, fontsize = label_font)
    sl_dur = Slider(slider_sub[1, 6], range = 100.0:25.0:1200.0, startvalue = 700.0, width = 70)
    Label(slider_sub[1, 7], text = lift(v -> "$(round(Int, v)) ms", sl_dur.value), width = 68, fontsize = label_font)

    Label(slider_sub[1, 9], text = "Noise:", font = :bold, fontsize = label_font)
    sl_noise = Slider(slider_sub[1, 10], range = 0.0:1.0:100.0, startvalue = 2.0, width = 70)
    Label(slider_sub[1, 11], text = lift(v -> "$(round(Int, v)) µV", sl_noise.value), width = 58, fontsize = label_font)

    colgap!(slider_sub, 6)
    colsize!(slider_sub, 4, Fixed(60))
    colsize!(slider_sub, 8, Fixed(60))
    rowgap!(ctrl_layout, 12)

    btn_sub = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    btn_play = Button(btn_sub[1, 1], label = "Play Saccade", buttoncolor = "#dcfce7", fontsize = btn_font, width = 150, height = 36)
    chk_fields = Checkbox(btn_sub[1, 3], checked = true)
    Label(btn_sub[1, 4], text = "Field Lines", fontsize = label_font)
    chk_mono = Checkbox(btn_sub[1, 6], checked = false)
    Label(btn_sub[1, 7], text = "Monopolar (CMR)", fontsize = label_font)

    colgap!(btn_sub, 6)
    colsize!(btn_sub, 2, Fixed(18))
    colsize!(btn_sub, 5, Fixed(18))

    ax_eog = Axis(
        fig[1:2, 2],
        title = "Bipolar hEOG Recording [µV]",
        xlabel = "Time [s]",
        ylabel = "Potential [µV]",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
    )
    xlims!(ax_eog, 0.0, 2.0)
    ylims!(ax_eog, -300, 320)
    hlines!(ax_eog, [0.0], color = :gray70, linestyle = :dash, linewidth = 1.2)

    t_vec = collect(range(0.0, 2.0, length = 600))
    base_noise = 0.8 .* randn(length(t_vec))
    common_eeg = @. 18.0 * sin(2π * 10.0 * t_vec) + 6.0 * cos(2π * 3.5 * t_vec)

    waves = lift(angle_deg, noise_amp, duration_ms) do deg, n_amp, dur_ms
        amp = 5.0 * deg
        t_on = 0.40
        t_off = t_on + dur_ms / 1000.0
        n_pts = length(t_vec)
        w_diff = Vector{Float64}(undef, n_pts)
        w_r = Vector{Float64}(undef, n_pts)
        w_l = Vector{Float64}(undef, n_pts)
        @inbounds for i = 1:n_pts
            tt = t_vec[i]
            pulse = (1.0 / (1.0 + exp(-(tt - t_on) / 0.015))) * (1.0 / (1.0 + exp((tt - t_off) / 0.015)))
            diff_i = muladd(amp, pulse, n_amp * base_noise[i])
            w_diff[i] = diff_i
            half_i = 0.5 * diff_i
            eeg_i = common_eeg[i]
            w_r[i] = half_i + eeg_i
            w_l[i] = -half_i + eeg_i
        end
        return (w_diff, w_r, w_l)
    end

    line_R = lines!(ax_eog, t_vec, lift(w -> w[2], waves), color = "#dc2626", linewidth = 1.8, label = "V(+) Right Temple", visible = false)
    line_L = lines!(ax_eog, t_vec, lift(w -> w[3], waves), color = "#2563eb", linewidth = 1.8, label = "V(–) Left Temple", visible = false)
    lines!(ax_eog, t_vec, lift(w -> w[1], waves), color = :black, linewidth = 2.4, label = "Vout Bipolar Diff")

    leg = axislegend(ax_eog, position = :lt, labelsize = 12, patchsize = (20, 10), backgroundcolor = RGBAf(1, 1, 1, 0.85))
    leg.blockscene.visible[] = false

    lines!(
        ax_eog,
        lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -300), Point2f(pt, 320)], playback_time),
        color = "#e11d48",
        linewidth = 2.0,
        linestyle = :dash,
    )

    hlines!(ax_eog, lift(deg -> [5.0 * deg], angle_deg), color = :gray60, linestyle = :dot, linewidth = 1.2)
    text!(
        ax_eog,
        lift(d -> 0.40 + (d / 2000.0), duration_ms),
        lift(deg -> 5.0 * deg + (deg >= 0 ? 25.0 : -35.0), angle_deg),
        text = lift(deg -> begin
            val = round(5.0 * deg, digits = 1)
            "ΔV = $(val >= 0 ? "+" : "")$(val) µV ($(round(deg, digits=1))°)"
        end, angle_deg),
        align = (:center, :center),
        font = :bold,
        fontsize = 14,
        color = :black,
    )

    is_animating = Ref(false)

    on(sl_angle.value) do val
        angle_deg[] = val
        is_animating[] || (current_display_angle[] = val)
    end
    on(val -> (duration_ms[] = val), sl_dur.value)
    on(val -> (noise_amp[] = val), sl_noise.value)
    on(val -> (show_fields[] = val), chk_fields.checked)
    on(chk_mono.checked) do val
        line_R.visible = val
        line_L.visible = val
        leg.blockscene.visible[] = val
    end

    on(btn_play.clicks) do _
        is_animating[] && return
        is_animating[] = true
        @async begin
            try
                target = angle_deg[]
                t_on = 0.40
                t_off = t_on + duration_ms[] / 1000.0
                n_frames = 60
                dt = 2.0 / n_frames
                for frame = 0:n_frames
                    (!is_animating[] || !events(fig).window_open[]) && break
                    t_curr = frame * dt
                    playback_time[] = t_curr
                    sig_on = 1.0 / (1.0 + exp(-(t_curr - t_on) / 0.015))
                    sig_off = 1.0 / (1.0 + exp((t_curr - t_off) / 0.015))
                    current_display_angle[] = target * sig_on * sig_off
                    sleep(1.0 / 30.0)
                end
                if events(fig).window_open[]
                    playback_time[] = -1.0
                    current_display_angle[] = target
                end
            finally
                is_animating[] = false
            end
        end
    end

    rowsize!(fig.layout, 1, Relative(0.80))
    rowsize!(fig.layout, 2, Relative(0.20))
    colsize!(fig.layout, 1, Relative(0.55))
    colsize!(fig.layout, 2, Relative(0.45))

    _display_demo_figure(fig)
    return fig, ax_eyes, ax_eog
end
