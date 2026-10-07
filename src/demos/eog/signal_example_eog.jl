"""
    signal_example_heog()

Interactive Horizontal Electrooculography (hEOG) Simulator — Bipolar Recording & Corneo-Retinal Dipole Demo.

Creates an interactive psychophysiology demonstration illustrating how eye rotations
are recorded via a bipolar difference amplifier at the outer canthi (temples).

## Educational Background

The human eye functions as a steady electrical dipole:
- **Corneo-Retinal Standing Potential**: The cornea is electrically positive (+) relative
  to the retina (-), generating an intrinsic dipole of approximately 0.4–1.0 mV driven by the
  metabolically active retinal pigment epithelium (RPE).
- **Volume Conduction**: As the eyeball rotates, this dipole re-orients relative to scalp/facial
  electrodes.
- **Bipolar Derivation**:
  - Electrode placed at the right outer canthus connects to the non-inverting (+) input: \$V_{(+)}\$.
  - Electrode placed at the left outer canthus connects to the inverting (–) input: \$V_{(-)}\$.
  - Difference amplifier output: \$V_{\\text{out}} = V_{(+)} - V_{(-)} = V_{\\text{right}} - V_{\\text{left}}\$.
- **Polarity Convention**:
  - Gaze **Right** (+ angle): Positive cornea approaches the right electrode (\$V_{(+)} > 0\$) and
    negative retina approaches the left electrode (\$V_{(-)} < 0\$) → **positive deflection**
    (e.g., +30° produces ≈ +150 µV).
  - Gaze **Left** (- angle): Positive cornea approaches the left electrode (\$V_{(-)} > 0\$) and
    negative retina approaches the right electrode (\$V_{(+)} < 0\$) → **negative deflection**
    (e.g., -15° produces ≈ -75 µV).
- **Linearity**: Over typical visual angles (±30°), the recorded potential is approximately proportional
  to the sine of the rotation angle (\$V \\approx 5.0\\,\\mu\\text{V}/\\text{degree}\$).

## Controls

| Control | Description |
|:---|:---|
| **Gaze Angle Slider** | Adjusts rotation angle from -45° (left) to +45° (right). Real-time update of dipoles & waveform. |
| **Duration Slider** | Adjusts fixation / saccade hold duration (100–1200 ms). |
| **Noise Slider** | Injects realistic physiological baseline noise (0–100 µV). |
| **▶ Play Saccade** | Plays an animated real-time saccade sequence (0° → target angle → 0°) with live cursor sweep. |

## Interactive Backend (GLMakie)

This interactive demo requires **GLMakie** for real-time reactivity, dragging sliders,
clicking preset buttons, and smooth saccadic playback animations:

```julia
using GLMakie
using EegFun

signal_example_heog()
```

If `CairoMakie` is active, a friendly warning will prompt you to activate `GLMakie`.

# Example
```julia
using GLMakie
using EegFun

signal_example_heog()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI.
- `ax_eyes::Axis`: The anatomical dipole & circuit schematic axis.
- `ax_eog::Axis`: The potential vs. time signal recording axis.
"""
function signal_example_heog()
    # Explicitly activate GLMakie before constructing the Figure to prevent CairoMakie from capturing it
    if isdefined(Main, :GLMakie)
        try
            Main.GLMakie.activate!(inline = false)
        catch
        end
    end

    if string(Makie.current_backend()) == "CairoMakie" || (!_is_glmakie_available() && !isdefined(Main, :GLMakie))
        @minimal_warning """
        Interactive demo requires GLMakie for live sliders, buttons, and animations.
        Currently active backend: $(Makie.current_backend())
        Please activate GLMakie in your session:
            using GLMakie
            GLMakie.activate!(inline = false)
        """
    end

    _set_window_title("Interactive Bipolar hEOG Recording Simulator")

    fig = Figure(
        size = (1350, 850),
        title = "Interactive Bipolar hEOG Recording Simulator",
        backgroundcolor = :white,
    )

    # ── Adaptive Font and UI Sizing ──────────────────────────────────────────
    function setup_adaptive_sizing(fig)
        title_font = Observable(18)
        label_font = Observable(14)
        tick_font = Observable(12)
        btn_font = Observable(13)

        on(fig.scene.viewport) do area
            scale = area.widths[1] / 1350
            title_font[] = max(14, round(Int, 18 * scale))
            label_font[] = max(11, round(Int, 14 * scale))
            tick_font[]  = max(10, round(Int, 12 * scale))
            btn_font[]   = max(11, round(Int, 13 * scale))
        end

        return title_font, label_font, tick_font, btn_font
    end

    title_font, label_font, tick_font, btn_font = setup_adaptive_sizing(fig)

    # ── Observables ──────────────────────────────────────────────────────────
    angle_deg = Observable(30.0)             # Target gaze angle in degrees
    current_display_angle = Observable(30.0) # Angle currently rendered (for animation)
    duration_ms = Observable(700.0)          # Saccade fixation / hold duration in milliseconds
    noise_amp = Observable(2.0)              # Additive noise standard deviation [µV]
    show_fields = Observable(true)           # Electric field lines toggle
    playback_time = Observable(-1.0)         # < 0: static, >= 0: time cursor position [s]

    # ── Geometry Constants ───────────────────────────────────────────────────
    CL = Point2f(-2.2, 0.5) # Left Eye center
    CR = Point2f(2.2, 0.5)  # Right Eye center
    R = 1.15                # Eyeball radius

    n_pts = 64
    theta_circ = range(0, 2π, length = n_pts)

    # ── Left Column [1, 1]: Anatomy & Circuit Schematic ─────────────────────
    ax_eyes = Axis(
        fig[1, 1],
        title = "Corneo-Retinal Dipoles & Bipolar Amplifier Circuit",
        titlesize = title_font,
        aspect = DataAspect(),
    )
    hidedecorations!(ax_eyes)
    hidespines!(ax_eyes)
    xlims!(ax_eyes, -5.2, 5.2)
    ylims!(ax_eyes, -3.2, 2.5)

    # ── Corneo-Retinal Potential Gradient (Deep Blue -> Purple -> Bright Red) ─
    eog_colormap = cgrad([RGBf(0.05, 0.22, 0.95), RGBf(0.55, 0.10, 0.75), RGBf(0.95, 0.08, 0.20)])

    function get_gradient_polys(C, deg)
        rad = deg * π / 180.0
        u = Point2f(sin(rad), cos(rad))    # dipole axis pointing from retina to cornea
        v = Point2f(cos(rad), -sin(rad))   # perpendicular axis

        n_slices = 64
        s_vals = range(-R, R, length = n_slices + 1)
        polys = Vector{Point2f}[]
        colors = RGBf[]

        for i in 1:n_slices
            s1 = s_vals[i]
            s2 = s_vals[i + 1]
            w1 = sqrt(max(0.0f0, Float32(R^2 - s1^2)))
            w2 = sqrt(max(0.0f0, Float32(R^2 - s2^2)))

            p1 = C .+ s1 .* u .- w1 .* v
            p2 = C .+ s1 .* u .+ w1 .* v
            p3 = C .+ s2 .* u .+ w2 .* v
            p4 = C .+ s2 .* u .- w2 .* v

            push!(polys, [p1, p2, p3, p4])
            t_color = (s1 + s2 + 2*R) / (4*R)
            push!(colors, eog_colormap[clamp(t_color, 0.0f0, 1.0f0)])
        end
        return polys, colors
    end

    slices_L = lift(current_display_angle) do deg
        get_gradient_polys(CL, deg)
    end
    slices_R = lift(current_display_angle) do deg
        get_gradient_polys(CR, deg)
    end

    # Render dynamic gradient slices for both eyeballs (64 thin polygons for smooth color transition)
    for i in 1:64
        poly!(ax_eyes, lift(s -> s[1][i], slices_L), color = lift(s -> s[2][i], slices_L), strokewidth = 0)
        poly!(ax_eyes, lift(s -> s[1][i], slices_R), color = lift(s -> s[2][i], slices_R), strokewidth = 0)
    end

    # Eyeball boundary circle rings
    poly!(ax_eyes, [CL .+ Point2f(R * cos(t), R * sin(t)) for t in theta_circ], color = :transparent, strokewidth = 2.0, strokecolor = "#1e293b")
    poly!(ax_eyes, [CR .+ Point2f(R * cos(t), R * sin(t)) for t in theta_circ], color = :transparent, strokewidth = 2.0, strokecolor = "#1e293b")

    # Feature generator function for each eye at rotation `deg`
    function get_eye_features(C, deg)
        rad = deg * π / 180.0
        u = Point2f(sin(rad), cos(rad))

        cornea_angles = range(-45.0, 45.0, length = 25)
        cornea_pts = Point2f[]
        for ca in cornea_angles
            ca_rad = (deg + ca) * π / 180.0
            r_bulge = R + 0.22 * cos(ca * π / 180.0)^2
            push!(cornea_pts, C .+ Point2f(r_bulge * sin(ca_rad), r_bulge * cos(ca_rad)))
        end
        for ca in reverse(cornea_angles)
            ca_rad = (deg + ca) * π / 180.0
            push!(cornea_pts, C .+ Point2f(R * sin(ca_rad), R * cos(ca_rad)))
        end

        iris_center  = C .+ 0.82 * R .* u
        pupil_center = C .+ 0.85 * R .* u
        dipole_start = C .- 0.95 * R .* u # negative retinal pole
        dipole_end   = C .+ 1.40 * R .* u # positive corneal pole
        plus_pos     = C .+ 1.62 * R .* u
        minus_pos    = C .- 1.20 * R .* u

        return (; cornea_pts, iris_center, pupil_center, dipole_start, dipole_end, plus_pos, minus_pos, u)
    end

    feat_L = lift(current_display_angle) do deg
        get_eye_features(CL, deg)
    end
    feat_R = lift(current_display_angle) do deg
        get_eye_features(CR, deg)
    end

    # Cornea bulge (electropositive translucent cyan dome)
    poly!(ax_eyes, lift(f -> f.cornea_pts, feat_L), color = RGBAf(0.2, 0.85, 0.98, 0.75), strokecolor = "#0284c7", strokewidth = 2.0)
    poly!(ax_eyes, lift(f -> f.cornea_pts, feat_R), color = RGBAf(0.2, 0.85, 0.98, 0.75), strokecolor = "#0284c7", strokewidth = 2.0)

    # Irises and pupils
    poly!(ax_eyes, lift(f -> [f.iris_center .+ Point2f(0.38 * cos(t), 0.38 * sin(t)) for t in theta_circ], feat_L), color = "#1e293b")
    poly!(ax_eyes, lift(f -> [f.iris_center .+ Point2f(0.38 * cos(t), 0.38 * sin(t)) for t in theta_circ], feat_R), color = "#1e293b")
    poly!(ax_eyes, lift(f -> [f.pupil_center .+ Point2f(0.18 * cos(t), 0.18 * sin(t)) for t in theta_circ], feat_L), color = :black)
    poly!(ax_eyes, lift(f -> [f.pupil_center .+ Point2f(0.18 * cos(t), 0.18 * sin(t)) for t in theta_circ], feat_R), color = :black)

    # Dipole arrows (green vector pointing along gaze direction)
    linesegments!(ax_eyes, lift(f -> [f.dipole_start, f.dipole_end], feat_L), color = "#10b981", linewidth = 2.5)
    linesegments!(ax_eyes, lift(f -> [f.dipole_start, f.dipole_end], feat_R), color = "#10b981", linewidth = 2.5)

    # Dipole charge signs: (+) at cornea, (–) at retina
    text!(ax_eyes, lift(f -> f.plus_pos[1], feat_L), lift(f -> f.plus_pos[2], feat_L), text = "+", align = (:center, :center), font = :bold, fontsize = 20, color = "#e11d48")
    text!(ax_eyes, lift(f -> f.plus_pos[1], feat_R), lift(f -> f.plus_pos[2], feat_R), text = "+", align = (:center, :center), font = :bold, fontsize = 20, color = "#e11d48")
    text!(ax_eyes, lift(f -> f.minus_pos[1], feat_L), lift(f -> f.minus_pos[2], feat_L), text = "–", align = (:center, :center), font = :bold, fontsize = 22, color = "#2563eb")
    text!(ax_eyes, lift(f -> f.minus_pos[1], feat_R), lift(f -> f.minus_pos[2], feat_R), text = "–", align = (:center, :center), font = :bold, fontsize = 22, color = "#2563eb")

    # Vertical reference lines (0° straight ahead)
    lines!(ax_eyes, [CL[1], CL[1]], [CL[2], CL[2] + 1.8 * R], color = :gray60, linestyle = :dash, linewidth = 1.2)
    lines!(ax_eyes, [CR[1], CR[1]], [CR[2], CR[2] + 1.8 * R], color = :gray60, linestyle = :dash, linewidth = 1.2)

    # Angle Arc on Right Eye
    arc_pts = lift(current_display_angle) do deg
        r_arc = 1.65
        deg_step = deg >= 0 ? 1.0 : -1.0
        angles = deg == 0 ? [0.0] : range(0.0, deg, step = deg_step)
        [CR .+ Point2f(r_arc * sin(a * π / 180), r_arc * cos(a * π / 180)) for a in angles]
    end
    lines!(ax_eyes, arc_pts, color = "#475569", linewidth = 1.8)

    angle_label_text = lift(current_display_angle) do deg
        s = deg >= 0 ? "+" : ""
        "$(s)$(round(deg, digits=1))°"
    end
    text!(
        ax_eyes,
        lift(deg -> CR[1] + 1.95 * sin(deg * π / 360), current_display_angle),
        lift(deg -> CR[2] + 1.95 * cos(deg * π / 360), current_display_angle),
        text = angle_label_text,
        align = (:center, :center),
        font = :bold,
        fontsize = 14,
        color = "#1e293b",
    )

    # Dipole Electric Field Lines (looping around Left Eye dipole)
    field_lines = lift(current_display_angle, show_fields) do deg, show
        !show && return Point2f[]
        pts = Point2f[]
        rad = deg * π / 180.0
        r0_list = [1.35, 1.85, 2.45]
        phi_range = range(0.12π, 0.88π, length = 35)

        for r0 in r0_list
            for sign in [1.0, -1.0]
                for phi in phi_range
                    r_val = r0 * sin(phi)^2
                    x_prime = sign * r_val * sin(phi)
                    y_prime = r_val * cos(phi)
                    xr = x_prime * cos(rad) + y_prime * sin(rad)
                    yr = -x_prime * sin(rad) + y_prime * cos(rad)
                    push!(pts, CL .+ Point2f(xr, yr))
                end
                push!(pts, Point2f(NaN, NaN))
            end
        end
        return pts
    end
    lines!(ax_eyes, field_lines, color = RGBAf(0.2, 0.5, 0.35, 0.45), linestyle = :dash, linewidth = 1.2)

    # Outer Canthi Temple Electrodes
    poly!(ax_eyes, Rect2f(-4.3, 0.2, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, -4.3, 1.0, text = "[-] Left Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")

    poly!(ax_eyes, Rect2f(4.08, 0.2, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, 4.19, 1.0, text = "[+] Right Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")

    # ── DIFFERENCE AMPLIFIER: Directly Below Both Eyes (Centered at x = 0.0) ──
    amp_tri = [Point2f(-0.7, -1.2), Point2f(-0.7, -2.2), Point2f(0.7, -1.7)]
    poly!(ax_eyes, amp_tri, color = "#f8fafc", strokewidth = 2.0, strokecolor = "#1e293b")
    text!(ax_eyes, -0.45, -1.40, text = "+", align = (:center, :center), font = :bold, fontsize = 18, color = "#dc2626")
    text!(ax_eyes, -0.45, -2.00, text = "–", align = (:center, :center), font = :bold, fontsize = 20, color = "#2563eb")
    text!(ax_eyes, 0.0, -0.55, text = "Difference Amplifier", align = (:center, :center), font = :bold, fontsize = 16, color = "#1e293b")

    # Output lead
    lines!(ax_eyes, [Point2f(0.7, -1.7), Point2f(1.2, -1.7)], color = "#1e293b", linewidth = 2.0)
    scatter!(ax_eyes, [Point2f(1.25, -1.7)], marker = :circle, markersize = 10, color = :white, strokecolor = "#1e293b", strokewidth = 2.0)

    # Lead wires routing cleanly into the centered amplifier:
    # Right wire to (+) input at (-0.7, -1.40), routing above the amplifier:
    wire_R = [
        Point2f(4.08, 0.5),
        Point2f(3.5, 0.5),
        Point2f(3.5, -0.85),
        Point2f(-1.0, -0.85),
        Point2f(-1.0, -1.40),
        Point2f(-0.70, -1.40),
    ]
    lines!(ax_eyes, wire_R, color = "#dc2626", linewidth = 2.0)

    # Left wire to (–) input at (-0.7, -2.00):
    wire_L = [
        Point2f(-4.08, 0.5),
        Point2f(-3.5, 0.5),
        Point2f(-3.5, -2.00),
        Point2f(-0.70, -2.00),
    ]
    lines!(ax_eyes, wire_L, color = "#2563eb", linewidth = 2.0)

    # Live Voltage Calculations and Readout Box (centered below the amplifier)
    v_calc = lift(current_display_angle) do deg
        v_diff = 5.0 * deg
        v_r = +(v_diff / 2.0)
        v_l = -(v_diff / 2.0)
        return (v_r, v_l, v_diff)
    end

    box_text = lift(v_calc) do (vr, vl, vd)
        sr = vr >= 0 ? "+" : ""
        sl = vl >= 0 ? "+" : ""
        sd = vd >= 0 ? "+" : ""
        "V(+) = $(sr)$(round(vr, digits=1)) µV   |   V(–) = $(sl)$(round(vl, digits=1)) µV   |   Vout = $(sd)$(round(vd, digits=1)) µV"
    end
    poly!(ax_eyes, Rect2f(-3.3, -2.88, 6.6, 0.56), color = "#f8fafc", strokewidth = 1.2, strokecolor = "#cbd5e1")
    text!(ax_eyes, 0.0, -2.60, text = box_text, align = (:center, :center), font = :bold, fontsize = 15.5, color = "#0f172a")

    # ── Left Column [2, 1]: Controls directly underneath both eyes ───────────
    ctrl_layout = fig[2, 1] = GridLayout(tellwidth = false, halign = :center)

    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)
    Label(slider_sub[1, 1], text = "Gaze:", font = :bold, fontsize = label_font)
    sl_angle = Slider(slider_sub[1, 2], range = -45.0:1.0:45.0, startvalue = 30.0, width = 75)
    lbl_angle = Label(slider_sub[1, 3], text = lift(v -> "$(round(v, digits=1))°", sl_angle.value), width = 58, fontsize = label_font)

    Label(slider_sub[1, 5], text = "Duration:", font = :bold, fontsize = label_font)
    sl_dur = Slider(slider_sub[1, 6], range = 100.0:25.0:1200.0, startvalue = 700.0, width = 70)
    lbl_dur = Label(slider_sub[1, 7], text = lift(v -> "$(round(Int, v)) ms", sl_dur.value), width = 68, fontsize = label_font)

    Label(slider_sub[1, 9], text = "Noise:", font = :bold, fontsize = label_font)
    sl_noise = Slider(slider_sub[1, 10], range = 0.0:1.0:100.0, startvalue = 2.0, width = 70)
    lbl_noise = Label(slider_sub[1, 11], text = lift(v -> "$(round(Int, v)) µV", sl_noise.value), width = 58, fontsize = label_font)

    colgap!(slider_sub, 6)
    colsize!(slider_sub, 4, Fixed(60))
    colsize!(slider_sub, 8, Fixed(60))

    rowgap!(ctrl_layout, 12)

    btn_play = Button(ctrl_layout[2, 1], label = "▶ Play Saccade", buttoncolor = "#dcfce7", fontsize = btn_font, width = 240, height = 38)

    # ── Right Column [1:2, 2]: Recording Signal Axis spanning full height ────
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

    t_vec = range(0.0, 2.0, length = 600)
    base_noise = 0.8 .* randn(length(t_vec))

    # Real-time waveform update
    waveform = lift(angle_deg, noise_amp, duration_ms) do deg, n_amp, dur_ms
        amp = 5.0 * deg
        t_on = 0.40
        t_off = t_on + dur_ms / 1000.0
        sig_on  = @. 1.0 / (1.0 + exp(-(t_vec - t_on) / 0.015))
        sig_off = @. 1.0 / (1.0 + exp((t_vec - t_off) / 0.015))
        pulse = sig_on .* sig_off
        return @. amp * pulse + n_amp * base_noise
    end
    lines!(ax_eog, t_vec, waveform, color = :black, linewidth = 2.2)

    # Time cursor line during playback
    cursor_x = lift(playback_time) do pt
        pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -300), Point2f(pt, 320)]
    end
    lines!(ax_eog, cursor_x, color = "#e11d48", linewidth = 2.0, linestyle = :dash)

    # Step guide line & amplitude label
    hlines!(ax_eog, lift(deg -> [5.0 * deg], angle_deg), color = :gray60, linestyle = :dot, linewidth = 1.2)
    v_text = lift(angle_deg) do deg
        val = round(5.0 * deg, digits = 1)
        sign_str = val >= 0 ? "+" : ""
        return "ΔV = $(sign_str)$(val) µV ($(round(deg, digits=1))°)"
    end
    text!(
        ax_eog,
        lift(d -> 0.40 + (d / 2000.0), duration_ms),
        lift(deg -> 5.0 * deg + (deg >= 0 ? 25.0 : -35.0), angle_deg),
        text = v_text,
        align = (:center, :center),
        font = :bold,
        fontsize = 14,
        color = :black,
    )

    # ── Callback Connections ─────────────────────────────────────────────────
    on(sl_angle.value) do val
        angle_deg[] = val
        current_display_angle[] = val
    end

    on(sl_dur.value) do val
        duration_ms[] = val
    end

    on(sl_noise.value) do val
        noise_amp[] = val
    end

    on(btn_play.clicks) do _
        @async begin
            target = angle_deg[]
            dur = duration_ms[] / 1000.0
            t_on = 0.40
            t_off = t_on + dur
            t_total = 2.0
            n_frames = 60
            fps = 30.0
            dt = t_total / n_frames
            for frame in 0:n_frames
                t_curr = frame * dt
                playback_time[] = t_curr
                sig_on  = 1.0 / (1.0 + exp(-(t_curr - t_on) / 0.015))
                sig_off = 1.0 / (1.0 + exp((t_curr - t_off) / 0.015))
                current_display_angle[] = target * sig_on * sig_off
                sleep(1.0 / fps)
            end
            playback_time[] = -1.0
            current_display_angle[] = target
        end
    end

    rowsize!(fig.layout, 1, Relative(0.80))
    rowsize!(fig.layout, 2, Relative(0.20))
    colsize!(fig.layout, 1, Relative(0.55))
    colsize!(fig.layout, 2, Relative(0.45))

    display(fig)

    return fig, ax_eyes, ax_eog
end
