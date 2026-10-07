"""
    signal_example_eog()

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
| **Noise Slider** | Injects realistic physiological baseline noise (0–5 µV). |
| **30° Right Preset** | Recreates the textbook +30° rightward saccade (+150 µV step). |
| **15° Left Preset** | Recreates the textbook -15° leftward saccade (-75 µV step). |
| **Center (0°) Preset** | Returns eyes to baseline resting position (0 µV). |
| **▶ Play Saccade** | Plays an animated real-time saccade sequence (0° → target angle → 0°) with live cursor sweep. |
| **Toggle Field Lines** | Shows/hides electric dipole field lines looping around the eye dipole. |

# Example
```julia
using EegFun
signal_example_eog()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI.
- `ax_eyes::Axis`: The anatomical dipole & circuit schematic axis.
- `ax_eog::Axis`: The potential vs. time signal recording axis.
"""
function signal_example_eog()
    fig = Figure(
        size = (1280, 800),
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
            scale = area.widths[1] / 1280
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
    noise_amp = Observable(0.5)              # Additive noise standard deviation [µV]
    show_fields = Observable(true)           # Electric field lines toggle
    playback_time = Observable(-1.0)         # < 0: static, >= 0: time cursor position [s]

    # ── Geometry Constants ───────────────────────────────────────────────────
    CL = Point2f(-2.1, 0.4) # Left Eye center
    CR = Point2f(2.1, 0.4)  # Right Eye center
    R = 1.15                # Eyeball radius

    n_pts = 64
    theta_circ = range(0, 2π, length = n_pts)

    # ── Left Axis: Anatomy & Circuit Schematic ──────────────────────────────
    ax_eyes = Axis(
        fig[1, 1],
        title = "Corneo-Retinal Dipoles & Bipolar Difference Circuit",
        titlesize = title_font,
        aspect = DataAspect(),
    )
    hidedecorations!(ax_eyes)
    hidespines!(ax_eyes)
    xlims!(ax_eyes, -5.2, 5.2)
    ylims!(ax_eyes, -2.8, 2.5)

    # Eyeball sclera bases
    poly!(
        ax_eyes,
        [CL .+ Point2f(R * cos(t), R * sin(t)) for t in theta_circ],
        color = "#f8fafc",
        strokewidth = 1.5,
        strokecolor = "#334155",
    )
    poly!(
        ax_eyes,
        [CR .+ Point2f(R * cos(t), R * sin(t)) for t in theta_circ],
        color = "#f8fafc",
        strokewidth = 1.5,
        strokecolor = "#334155",
    )

    # Feature generator function for each eye at rotation `deg`
    function get_eye_features(C, deg)
        rad = deg * π / 180.0
        u = Point2f(sin(rad), cos(rad)) # gaze direction vector

        # Retinal backing (rear semicircle opposite gaze direction)
        ret_angles = range(deg - 90, deg + 90, length = 32)
        ret_pts = [C .- Point2f(R * sin(a * π / 180), R * cos(a * π / 180)) for a in ret_angles]
        push!(ret_pts, C)

        # Cornea bulge (protruding convex dome along gaze vector)
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

        return (; ret_pts, cornea_pts, iris_center, pupil_center, dipole_start, dipole_end, plus_pos, minus_pos, u)
    end

    feat_L = lift(current_display_angle) do deg
        get_eye_features(CL, deg)
    end
    feat_R = lift(current_display_angle) do deg
        get_eye_features(CR, deg)
    end

    # Retinal backing (electronegative indigo/violet shell)
    poly!(ax_eyes, lift(f -> f.ret_pts, feat_L), color = "#4338ca", strokewidth = 0)
    poly!(ax_eyes, lift(f -> f.ret_pts, feat_R), color = "#4338ca", strokewidth = 0)

    # Cornea bulge (electropositive translucent cyan dome)
    poly!(ax_eyes, lift(f -> f.cornea_pts, feat_L), color = RGBAf(0.2, 0.8, 0.95, 0.7), strokecolor = "#0284c7", strokewidth = 2.0)
    poly!(ax_eyes, lift(f -> f.cornea_pts, feat_R), color = RGBAf(0.2, 0.8, 0.95, 0.7), strokecolor = "#0284c7", strokewidth = 2.0)

    # Irises and pupils
    poly!(ax_eyes, lift(f -> [f.iris_center .+ Point2f(0.38 * cos(t), 0.38 * sin(t)) for t in theta_circ], feat_L), color = "#374151")
    poly!(ax_eyes, lift(f -> [f.iris_center .+ Point2f(0.38 * cos(t), 0.38 * sin(t)) for t in theta_circ], feat_R), color = "#374151")
    poly!(ax_eyes, lift(f -> [f.pupil_center .+ Point2f(0.18 * cos(t), 0.18 * sin(t)) for t in theta_circ], feat_L), color = :black)
    poly!(ax_eyes, lift(f -> [f.pupil_center .+ Point2f(0.18 * cos(t), 0.18 * sin(t)) for t in theta_circ], feat_R), color = :black)

    # Dipole arrows (green vector pointing along gaze direction)
    linesegments!(ax_eyes, lift(f -> [f.dipole_start, f.dipole_end], feat_L), color = "#10b981", linewidth = 2.5)
    linesegments!(ax_eyes, lift(f -> [f.dipole_start, f.dipole_end], feat_R), color = "#10b981", linewidth = 2.5)

    # Dipole charge signs: (+) at cornea, (–) at retina
    text!(ax_eyes, lift(f -> f.plus_pos[1], feat_L), lift(f -> f.plus_pos[2], feat_L), text = "+", align = (:center, :center), font = :bold, fontsize = 18, color = "#e11d48")
    text!(ax_eyes, lift(f -> f.plus_pos[1], feat_R), lift(f -> f.plus_pos[2], feat_R), text = "+", align = (:center, :center), font = :bold, fontsize = 18, color = "#e11d48")
    text!(ax_eyes, lift(f -> f.minus_pos[1], feat_L), lift(f -> f.minus_pos[2], feat_L), text = "–", align = (:center, :center), font = :bold, fontsize = 20, color = "#2563eb")
    text!(ax_eyes, lift(f -> f.minus_pos[1], feat_R), lift(f -> f.minus_pos[2], feat_R), text = "–", align = (:center, :center), font = :bold, fontsize = 20, color = "#2563eb")

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
    poly!(ax_eyes, Rect2f(-4.3, 0.1, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, -4.3, 0.9, text = "[-] Left Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")

    poly!(ax_eyes, Rect2f(4.08, 0.1, 0.22, 0.6), color = "#334155", strokewidth = 1, strokecolor = :black)
    text!(ax_eyes, 4.19, 0.9, text = "[+] Right Temple", align = (:center, :center), font = :bold, fontsize = 12, color = "#334155")

    # Difference Amplifier Triangle
    amp_tri = [Point2f(0.6, -0.85), Point2f(0.6, -1.95), Point2f(2.1, -1.40)]
    poly!(ax_eyes, amp_tri, color = "#f8fafc", strokewidth = 2.0, strokecolor = "#1e293b")
    text!(ax_eyes, 0.85, -1.05, text = "+", align = (:center, :center), font = :bold, fontsize = 18, color = "#dc2626")
    text!(ax_eyes, 0.85, -1.75, text = "–", align = (:center, :center), font = :bold, fontsize = 20, color = "#2563eb")
    text!(ax_eyes, 1.35, -0.62, text = "Difference Amplifier", align = (:center, :center), font = :bold, fontsize = 13, color = "#1e293b")

    # Amplifier output lead & terminal
    lines!(ax_eyes, [Point2f(2.1, -1.40), Point2f(2.6, -1.40)], color = "#1e293b", linewidth = 2.0)
    scatter!(ax_eyes, [Point2f(2.65, -1.40)], marker = :circle, markersize = 10, color = :white, strokecolor = "#1e293b", strokewidth = 2.0)

    # Lead wires
    # Right temple lead (red wire to non-inverting + input)
    wire_R = [Point2f(4.08, 0.4), Point2f(3.6, 0.4), Point2f(3.6, -1.05), Point2f(0.6, -1.05)]
    lines!(ax_eyes, wire_R, color = "#dc2626", linewidth = 2.0)

    # Left temple lead (blue wire to inverting - input)
    wire_L = [Point2f(-4.08, 0.4), Point2f(-3.6, 0.4), Point2f(-3.6, -2.25), Point2f(0.6, -2.25), Point2f(0.6, -1.75)]
    lines!(ax_eyes, wire_L, color = "#2563eb", linewidth = 2.0)

    # Live Voltage Calculations and Readout Box
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
        "V(+) = $(sr)$(round(vr, digits=1)) µV\n" *
        "V(–) = $(sl)$(round(vl, digits=1)) µV\n" *
        "─────────────────────\n" *
        "Vout = V(+) – V(–) = $(sd)$(round(vd, digits=1)) µV"
    end
    text!(ax_eyes, -1.8, -1.50, text = box_text, align = (:center, :center), font = :bold, fontsize = 12, color = "#1e293b")

    # ── Right Axis: Signal Potential vs Time ─────────────────────────────────
    ax_eog = Axis(
        fig[1, 2],
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
    ylims!(ax_eog, -220, 240)
    hlines!(ax_eog, [0.0], color = :gray70, linestyle = :dash, linewidth = 1.2)

    t_vec = range(0.0, 2.0, length = 600)
    base_noise = 0.3 .* randn(length(t_vec))

    # Real-time waveform update
    waveform = lift(angle_deg, noise_amp) do deg, n_amp
        amp = 5.0 * deg
        sig_on  = @. 1.0 / (1.0 + exp(-(t_vec - 0.45) / 0.015))
        sig_off = @. 1.0 / (1.0 + exp((t_vec - 1.15) / 0.015))
        pulse = sig_on .* sig_off
        return @. amp * pulse + n_amp * base_noise
    end
    lines!(ax_eog, t_vec, waveform, color = "#2563eb", linewidth = 2.5)

    # Time cursor line during playback
    cursor_x = lift(playback_time) do pt
        pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -220), Point2f(pt, 240)]
    end
    lines!(ax_eog, cursor_x, color = "#e11d48", linewidth = 2.0, linestyle = :dash)

    # Step guide line & amplitude label
    hlines!(ax_eog, lift(deg -> [5.0 * deg], angle_deg), color = "#93c5fd", linestyle = :dot, linewidth = 1.5)
    v_text = lift(angle_deg) do deg
        val = round(5.0 * deg, digits = 1)
        sign_str = val >= 0 ? "+" : ""
        return "ΔV = $(sign_str)$(val) µV ($(round(deg, digits=1))°)"
    end
    text!(
        ax_eog,
        0.8,
        lift(deg -> 5.0 * deg + (deg >= 0 ? 25.0 : -35.0), angle_deg),
        text = v_text,
        align = (:center, :center),
        font = :bold,
        fontsize = 14,
        color = "#1e40af",
    )

    # ── Controls Layout ──────────────────────────────────────────────────────
    ctrl_layout = fig[2, 1:2] = GridLayout()

    # Preset Buttons Row
    btn_sub = ctrl_layout[1, 1] = GridLayout()
    btn_right30 = Button(btn_sub[1, 1], label = "30° Right (+150 µV)", buttoncolor = "#e0f2fe", fontsize = btn_font)
    btn_left15  = Button(btn_sub[1, 2], label = "15° Left (-75 µV)", buttoncolor = "#fce7f3", fontsize = btn_font)
    btn_center  = Button(btn_sub[1, 3], label = "Center (0°)", buttoncolor = "#f1f5f9", fontsize = btn_font)
    btn_play    = Button(btn_sub[1, 4], label = "▶ Play Saccade", buttoncolor = "#dcfce7", fontsize = btn_font)
    btn_fields  = Button(btn_sub[1, 5], label = "Toggle Field Lines", buttoncolor = "#f8fafc", fontsize = btn_font)

    # Sliders Row
    slider_sub = ctrl_layout[2, 1] = GridLayout()
    Label(slider_sub[1, 1], text = "Gaze Angle:", font = :bold, fontsize = label_font)
    sl_angle = Slider(slider_sub[1, 2], range = -45.0:1.0:45.0, startvalue = 30.0, width = 240)
    lbl_angle = Label(slider_sub[1, 3], text = lift(v -> "$(round(v, digits=1))°", sl_angle.value), width = 60, fontsize = label_font)

    Label(slider_sub[1, 4], text = "Noise (µV):", font = :bold, fontsize = label_font)
    sl_noise = Slider(slider_sub[1, 5], range = 0.0:0.1:5.0, startvalue = 0.5, width = 150)
    lbl_noise = Label(slider_sub[1, 6], text = lift(v -> "$(round(v, digits=1)) µV", sl_noise.value), width = 60, fontsize = label_font)

    # ── Callback Connections ─────────────────────────────────────────────────
    on(sl_angle.value) do val
        angle_deg[] = val
        current_display_angle[] = val
    end

    on(sl_noise.value) do val
        noise_amp[] = val
    end

    on(btn_right30.clicks) do _
        set_close_to!(sl_angle, 30.0)
        angle_deg[] = 30.0
        current_display_angle[] = 30.0
    end

    on(btn_left15.clicks) do _
        set_close_to!(sl_angle, -15.0)
        angle_deg[] = -15.0
        current_display_angle[] = -15.0
    end

    on(btn_center.clicks) do _
        set_close_to!(sl_angle, 0.0)
        angle_deg[] = 0.0
        current_display_angle[] = 0.0
    end

    on(btn_fields.clicks) do _
        show_fields[] = !show_fields[]
    end

    on(btn_play.clicks) do _
        @async begin
            target = angle_deg[]
            n_frames = 60
            fps = 30.0
            dt = 2.0 / n_frames
            for frame in 0:n_frames
                t_curr = frame * dt
                playback_time[] = t_curr
                sig_on  = 1.0 / (1.0 + exp(-(t_curr - 0.45) / 0.015))
                sig_off = 1.0 / (1.0 + exp((t_curr - 1.15) / 0.015))
                current_display_angle[] = target * sig_on * sig_off
                sleep(1.0 / fps)
            end
            playback_time[] = -1.0
            current_display_angle[] = target
        end
    end

    rowsize!(fig.layout, 1, Relative(0.80))
    rowsize!(fig.layout, 2, Relative(0.20))

    display(fig)

    return fig, ax_eyes, ax_eog
end

# Alias for convenience
const signal_example_heog = signal_example_eog
