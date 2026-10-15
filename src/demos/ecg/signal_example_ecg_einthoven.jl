"""
    signal_example_ecg_einthoven()

Interactive Electrocardiography (ECG) Simulator — Einthoven's Triangle
& Vectorcardiogram Demo.

Creates an interactive psychophysiology demonstration illustrating the cellular-dipole
origin of the electrocardiogram, Einthoven's triangle, and bipolar lead
derivations (Lead I, Lead II, Lead III).

## Educational Background

### 1. The Heart as an Electrical Dipole
During each cardiac cycle, waves of electrical depolarization and repolarization sweep
through myocardial tissue:
- **Atrial Depolarization (P wave)**: Initiated by the Sinoatrial (SA) node in the upper right
  atrium, spreading downward and leftward toward the Atrioventricular (AV) node (~+45° to +60°).
- **Septal Depolarization (Q wave)**: Depolarization crosses the interventricular septum from the
  left bundle branch to the right (~-125°), producing a small initial negative deflection in Lead II.
- **Ventricular Apical Depolarization (R wave)**: Depolarization spreads rapidly through the Purkinje
  network into the massive muscular left ventricular apex (~+60° normal axis), generating the prominent
  tall R wave.
- **Ventricular Basal Depolarization (S wave)**: Late depolarization of the posterobasal left ventricle
  and pulmonary conus pointing superiorly and rightward (~-105°), producing the S wave.
- **Ventricular Repolarization (T wave)**: Unlike depolarization, repolarization spreads from the
  epicardium to the endocardium (reversing both direction and electrical sign), yielding an upright
  T wave concordant with the QRS complex in normal physiological states.

### 2. Einthoven's Triangle & Bipolar Leads
Willem Einthoven modeled the human torso as a homogeneous volume conductor with the heart
at the center of an equilateral triangle formed by the right arm (RA), left arm (LA), and left leg (LL):
- **Lead I**: Potential difference between Left Arm (+) and Right Arm (–), axis 0°.
- **Lead II**: Potential difference between Left Leg (+) and Right Arm (–), axis +60°.
- **Lead III**: Potential difference between Left Leg (+) and Left Arm (–), axis +120°.

## Controls

| Control | Description |
|:---|:---|
| **Heart Axis (α) Slider** | Rotates the heart's mean electrical axis (–45° to +120°). Demonstrates axis shifts and waveform morphology. |
| **HR (BPM) Slider** | Adjusts heart rate over an expanded range (10–220 BPM). Demonstrates extreme bradycardia to tachycardia. |
| **Noise Slider** | Injects high-frequency somatic muscle EMG & line noise over an expanded range (0–200, unitless). |
| **Play / Pause** | Animates the cardiac cycle in real-time with sweeping dipole vector and synchronized ECG playhead. |
| **Lead Buttons** | Toggles active lead: `Lead I`, `Lead II (Std)`, `Lead III`, or `All 3 (I, II, III)`. |
| **Vectorcardiogram Checkbox** | Toggles display of the complete 3-lobed Vectorcardiogram (P, QRS, T loops). |
| **Projection Checkbox** | Toggles orthogonal projection line and scalar voltage bar on active lead axis. |

## Interactive Backend (GLMakie)

This interactive demo requires **GLMakie** for real-time reactivity, dragging sliders,
and smooth cardiac cycle playback animations:

```julia
using GLMakie
using EegFun

signal_example_ecg_einthoven()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI.
- `ax_torso::Axis`: Anatomical torso, Einthoven's triangle, cardiac dipole & Vectorcardiogram schematic axis.
- `ax_ecg::Axis`: High-resolution calibrated ECG strip axis.
"""
function signal_example_ecg_einthoven()
    # Explicitly activate GLMakie before constructing the Figure, in case another backend (e.g. CairoMakie) is active
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

    _set_window_title("Interactive ECG Simulator — Einthoven's Triangle & Vectorcardiogram")

    fig = Figure(
        size = (1400, 860),
        title = "Interactive ECG Simulator — Einthoven's Triangle & Vectorcardiogram",
        backgroundcolor = :white,
    )

    # ── Adaptive Font and UI Sizing ──────────────────────────────────────────
    function setup_adaptive_sizing(fig)
        title_font = Observable(16)
        label_font = Observable(13.0)
        tick_font  = Observable(11.0)
        btn_font   = Observable(12.5)

        on(fig.scene.viewport) do area
            scale = area.widths[1] / 1400
            title_font[] = max(13, round(Int, 16 * scale))
            label_font[] = max(11, round(Int, 13.0 * scale))
            tick_font[]  = max(10, round(Int, 11.0 * scale))
            btn_font[]   = max(11, round(Int, 12.5 * scale))
        end

        return title_font, label_font, tick_font, btn_font
    end

    title_font, label_font, tick_font, btn_font = setup_adaptive_sizing(fig)

    # ── Observables ──────────────────────────────────────────────────────────
    axis_alpha_deg   = Observable(60.0)      # Mean QRS electrical axis in degrees (-45° to +120°)
    heart_rate_bpm   = Observable(72.0)      # Heart rate in beats per minute
    noise_level      = Observable(15.0)      # High-frequency noise level (0-200)
    selected_lead    = Observable(:lead2)    # :lead1, :lead2, :lead3, :all3
    show_vcg_loop    = Observable(true)      # Show Vectorcardiogram loop trajectory
    show_projection  = Observable(true)      # Show orthogonal projection onto active lead
    playback_time    = Observable(-1.0)      # < 0: static R-peak inspection, >= 0: live time cursor [s]
    is_playing       = Observable(false)     # Continuous heartbeat playback flag

    # ── Geometry & Anatomical Constants ──────────────────────────────────────
    # Equilateral Einthoven Triangle positions (Standard limb placement)
    # RA: (-2.60, 1.80), LA: (+2.60, 1.80), LL: (0.00, -2.70), side ≈ 5.20
    pos_ra = Point2f(-2.60, 1.80)
    pos_la = Point2f(2.60, 1.80)
    pos_ll = Point2f(0.00, -2.70)
    pos_rl = Point2f(-1.10, -2.70) # Right leg ground

    # Anatomical Heart Center (dipole origin, slightly left of thoracic midline)
    heart_center = Point2f(0.25, 0.35)
    dipole_scale = 1.35 # Visual vector scale on torso plot

    # ── Left Column [1, 1]: Torso, Einthoven Triangle & Vector Schematic ─────
    ax_torso = Axis(
        fig[1, 1],
        title = "Einthoven's Triangle & Vectorcardiogram",
        titlesize = title_font,
        aspect = DataAspect(),
    )
    hidedecorations!(ax_torso)
    hidespines!(ax_torso)
    xlims!(ax_torso, -3.85, 3.85)
    ylims!(ax_torso, -3.65, 3.45)

    # ── Stylized Human Torso Silhouette ──────────────────────────────────────
    torso_half = [
        Point2f(0.55, 3.10),   # neck top
        Point2f(0.55, 2.70),   # neck base
        Point2f(1.60, 2.60),   # clavicle slope
        Point2f(2.90, 2.30),   # shoulder
        Point2f(3.40, 1.80),   # outer deltoid
        Point2f(3.20, 0.90),   # armpit outer
        Point2f(2.55, 0.70),   # armpit inner
        Point2f(2.40, -0.40),  # ribcage lateral
        Point2f(2.15, -1.50),  # waist
        Point2f(2.45, -2.60),  # iliac crest / hip
        Point2f(2.20, -3.30),  # groin lateral
        Point2f(0.00, -3.45),  # pubis / perineum
    ]
    torso_left = [Point2f(-p[1], p[2]) for p in reverse(torso_half)]
    torso_pts  = [torso_left; torso_half]

    poly!(ax_torso, torso_pts, color = RGBAf(0.96, 0.97, 0.98, 1.0), strokecolor = "#94a3b8", strokewidth = 1.8)

    # Clavicles and Sternum anatomical reference lines
    lines!(ax_torso, [Point2f(0.2, 2.50), Point2f(1.2, 2.45), Point2f(2.4, 2.20)], color = "#cbd5e1", linewidth = 2.0)
    lines!(ax_torso, [Point2f(-0.2, 2.50), Point2f(-1.2, 2.45), Point2f(-2.4, 2.20)], color = "#cbd5e1", linewidth = 2.0)
    poly!(ax_torso, Rect2f(-0.18, 0.80, 0.36, 1.55), color = "#e2e8f0", strokecolor = "#cbd5e1", strokewidth = 1.0)
    text!(ax_torso, 0.0, 1.55, text = "Sternum", align = (:center, :center), font = :bold, fontsize = 8.5, color = "#94a3b8")

    # ── Stylized Anatomical Heart at (0.25, 0.35) ────────────────────────────
    # Great vessels (Aorta arch in crimson, Pulmonary trunk in blue)
    poly!(ax_torso, [Point2f(0.12, 1.00), Point2f(0.35, 1.35), Point2f(0.60, 1.15), Point2f(0.40, 0.85)], color = "#dc2626", strokecolor = "#991b1b", strokewidth = 1.0)
    poly!(ax_torso, [Point2f(-0.05, 0.90), Point2f(0.15, 1.15), Point2f(0.30, 0.80), Point2f(0.10, 0.65)], color = "#2563eb", strokecolor = "#1d4ed8", strokewidth = 1.0)

    # Atria (upper chambers)
    poly!(ax_torso, [Point2f(-0.45, 0.85), Point2f(-0.05, 0.90), Point2f(0.00, 0.55), Point2f(-0.35, 0.50)], color = "#fecaca", strokecolor = "#f87171", strokewidth = 1.2) # RA
    poly!(ax_torso, [Point2f(0.20, 0.85), Point2f(0.70, 0.80), Point2f(0.65, 0.50), Point2f(0.15, 0.55)], color = "#fecaca", strokecolor = "#f87171", strokewidth = 1.2) # LA

    # Ventricles (conical muscular myocardium pointing down-left toward apex)
    ventricle_pts = [
        Point2f(-0.35, 0.50),
        Point2f(0.00, 0.55),
        Point2f(0.65, 0.50),
        Point2f(0.75, 0.15),
        Point2f(0.68, -0.30), # Left ventricular apex
        Point2f(0.40, -0.38), # Apex tip
        Point2f(0.05, -0.20),
        Point2f(-0.25, 0.15),
    ]
    poly!(ax_torso, ventricle_pts, color = "#ef4444", strokecolor = "#b91c1c", strokewidth = 1.8)

    # Interventricular septum line
    lines!(ax_torso, [Point2f(0.12, 0.52), Point2f(0.45, -0.35)], color = "#991b1b", linewidth = 1.5, linestyle = :dash)

    # Conduction nodes: SA Node (yellow glowing dot) and AV Node (emerald dot)
    scatter!(ax_torso, [-0.28], [0.82], color = "#facc15", strokecolor = "#ca8a04", strokewidth = 1.5, markersize = 12)
    text!(ax_torso, -0.32, 0.92, text = "SA Node", align = (:right, :center), font = :bold, fontsize = 8.5, color = "#ca8a04")

    scatter!(ax_torso, [0.08], [0.48], color = "#10b981", strokecolor = "#047857", strokewidth = 1.5, markersize = 10)
    text!(ax_torso, 0.12, 0.38, text = "AV Node", align = (:left, :center), font = :bold, fontsize = 8.5, color = "#047857")

    # ── Dipole Physics & Vector Function ─────────────────────────────────────
    # Computes 2D dipole components (Px, Py) for an elapsed time t_in_beat in [0, T]
    # Coordinate convention: Px points patient's left (+x); Py points inferiorly/feet (-y)
    function compute_cardiac_dipole(t_in_beat::Float64, bpm_val::Float64, alpha_val::Float64)
        T = 60.0 / bpm_val
        T0 = 60.0 / 72.0 # 0.8333 s reference period

        # For tachycardia (BPM > 72), scale wave timings to fit inside T
        k = bpm_val > 72.0 ? (T / T0) : 1.0

        # QT interval slightly lengthens with bradycardia (Bazett-like, capped)
        qt_scale = bpm_val > 72.0 ? k : min(1.35, sqrt(T / T0))

        # Physiological wave timings (in seconds):
        # P wave: ~80-100 ms duration
        tP = 0.12 * k; wP = 0.022 * k
        ampP = 0.18; thetaP = 48.0

        # Q wave: septal depolarization (~30-40 ms)
        tQ = 0.22 * k; wQ = 0.010 * k
        ampQ = 0.16; thetaQ = -125.0 + (alpha_val - 60.0) * 0.25

        # R wave: ventricular apical depolarization (~70-90 ms QRS total)
        tR = 0.25 * k; wR = 0.014 * k
        ampR = 1.35; thetaR = alpha_val

        # S wave: ventricular basal depolarization
        tS = 0.28 * k; wS = 0.012 * k
        ampS = 0.32; thetaS = -105.0 + (alpha_val - 60.0) * 0.35

        # T wave: ventricular repolarization
        tT = (0.25 + 0.18 * qt_scale) * (bpm_val > 72.0 ? k : 1.0)
        wT = 0.045 * qt_scale * (bpm_val > 72.0 ? k : 1.0)
        ampT = 0.38; thetaT = 45.0 + (alpha_val - 60.0) * 0.70

        # Beyond the end of the active beat (diastolic TP pause), dipole is 0
        t_active_end = tT + 3.0 * wT
        if t_in_beat > t_active_end && bpm_val <= 72.0
            return 0.0, 0.0
        end

        gP = exp(-((t_in_beat - tP) / wP)^2)
        gQ = exp(-((t_in_beat - tQ) / wQ)^2)
        gR = exp(-((t_in_beat - tR) / wR)^2)
        gS = exp(-((t_in_beat - tS) / wS)^2)
        gT = exp(-((t_in_beat - tT) / wT)^2)

        Px = ampP*gP*cosd(thetaP) + ampQ*gQ*cosd(thetaQ) + ampR*gR*cosd(thetaR) + ampS*gS*cosd(thetaS) + ampT*gT*cosd(thetaT)
        Py = -(ampP*gP*sind(thetaP) + ampQ*gQ*sind(thetaQ) + ampR*gR*sind(thetaR) + ampS*gS*sind(thetaS) + ampT*gT*sind(thetaT))

        return Px, Py
    end

    # Current cardiac phase for vector display
    # When paused (playback_time < 0), display the R-peak
    current_dipole = lift(playback_time, heart_rate_bpm, axis_alpha_deg) do pt, bpm, alpha
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        t_in = pt < 0.0 ? (0.25 * k) : mod(pt, T)
        compute_cardiac_dipole(t_in, bpm, alpha)
    end

    # ── Precomputed Vectorcardiogram ─────────────────────────────────────────
    vcg_pts = lift(axis_alpha_deg, heart_rate_bpm) do alpha, bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        qt_scale = bpm > 72.0 ? k : min(1.35, sqrt(T / (60.0 / 72.0)))
        tT = (0.25 + 0.18 * qt_scale) * (bpm > 72.0 ? k : 1.0)
        wT = 0.045 * qt_scale * (bpm > 72.0 ? k : 1.0)
        t_end = min(T, tT + 2.8 * wT)

        n_vcg = 160
        pts = Point2f[]
        for i in 0:n_vcg
            t_sample = (i / n_vcg) * t_end
            px, py = compute_cardiac_dipole(t_sample, bpm, alpha)
            push!(pts, heart_center .+ Point2f(dipole_scale * px, dipole_scale * py))
        end
        push!(pts, heart_center)
        return pts
    end

    lines!(
        ax_torso,
        vcg_pts,
        color = "#06b6d4",
        linewidth = 2.0,
        linestyle = :solid,
        visible = show_vcg_loop,
        label = "Vectorcardiogram",
    )

    # ── Active Dipole Vector Arrow ───────────────────────────────────────────
    dipole_tip = lift(current_dipole) do (px, py)
        heart_center .+ Point2f(dipole_scale * px, dipole_scale * py)
    end

    # Dipole vector line from heart center to dipole tip
    linesegments!(
        ax_torso,
        lift(tip -> [heart_center, tip], dipole_tip),
        color = "#e11d48",
        linewidth = 3.2,
    )
    scatter!(
        ax_torso,
        lift(tip -> [tip[1]], dipole_tip),
        lift(tip -> [tip[2]], dipole_tip),
        color = "#e11d48",
        strokecolor = "#9f1239",
        strokewidth = 1.5,
        markersize = 11,
    )

    # ── Einthoven Triangle Lead Lines & Active Highlight ─────────────────────
    line_lead1 = [pos_ra, pos_la]
    line_lead2 = [pos_ra, pos_ll]
    line_lead3 = [pos_la, pos_ll]

    # Lead line colors depending on selection
    col_lead1 = lift(s -> (s == :lead1 || s == :all3) ? "#0284c7" : "#cbd5e1", selected_lead)
    col_lead2 = lift(s -> (s == :lead2 || s == :all3) ? "#e11d48" : "#cbd5e1", selected_lead)
    col_lead3 = lift(s -> (s == :lead3 || s == :all3) ? "#10b981" : "#cbd5e1", selected_lead)

    lw_lead1 = lift(s -> (s == :lead1 || s == :all3) ? 3.0 : 1.4, selected_lead)
    lw_lead2 = lift(s -> (s == :lead2 || s == :all3) ? 3.0 : 1.4, selected_lead)
    lw_lead3 = lift(s -> (s == :lead3 || s == :all3) ? 3.0 : 1.4, selected_lead)

    lines!(ax_torso, line_lead1, color = col_lead1, linewidth = lw_lead1)
    lines!(ax_torso, line_lead2, color = col_lead2, linewidth = lw_lead2)
    lines!(ax_torso, line_lead3, color = col_lead3, linewidth = lw_lead3)

    # Lead direction indicator arrowheads
    text!(ax_torso, 0.0, pos_la[2] + 0.16, text = "Lead I (0°) ➔", align = (:center, :bottom), font = :bold, fontsize = 10, color = col_lead1)
    text!(ax_torso, (pos_ra[1] + pos_ll[1])/2 - 0.22, (pos_ra[2] + pos_ll[2])/2, text = "Lead II (+60°) ➔", align = (:right, :center), font = :bold, fontsize = 10, color = col_lead2)
    text!(ax_torso, (pos_la[1] + pos_ll[1])/2 + 0.22, (pos_la[2] + pos_ll[2])/2, text = "➔ Lead III (+120°)", align = (:left, :center), font = :bold, fontsize = 10, color = col_lead3)

    # ── Orthogonal Vector Projection onto Active Lead Axis ───────────────────
    # Active unit lead vector
    lead_unit_vec = lift(selected_lead) do s
        if s == :lead1
            Point2f(1.0, 0.0)
        elseif s == :lead2
            Point2f(cosd(60.0), -sind(60.0))
        elseif s == :lead3
            Point2f(cosd(120.0), -sind(120.0))
        else # :all3 default to Lead II
            Point2f(cosd(60.0), -sind(60.0))
        end
    end

    # Projected scalar and point along the lead line passing through heart_center
    proj_point = lift(current_dipole, lead_unit_vec) do (px, py), u
        # scalar dot product: Px * ux + Py * uy
        v_proj = px * u[1] + py * u[2]
        return heart_center .+ Point2f(dipole_scale * v_proj * u[1], dipole_scale * v_proj * u[2])
    end

    # Dashed perpendicular projection line from dipole tip to active lead axis
    linesegments!(
        ax_torso,
        lift((tip, pproj) -> [tip, pproj], dipole_tip, proj_point),
        color = "#64748b",
        linewidth = 1.8,
        linestyle = :dash,
        visible = show_projection,
    )

    # Projected scalar vector along the active lead axis
    linesegments!(
        ax_torso,
        lift(pproj -> [heart_center, pproj], proj_point),
        color = "#2563eb",
        linewidth = 4.0,
        visible = show_projection,
    )
    scatter!(
        ax_torso,
        lift(p -> [p[1]], proj_point),
        lift(p -> [p[2]], proj_point),
        color = "#2563eb",
        markersize = 9,
        visible = show_projection,
    )

    # ── Electrode Pads on Torso (Concentric Rings & Polarity Badges) ─────────
    # Outer adhesive rings
    scatter!(ax_torso, [pos_ra, pos_la, pos_ll], color = "#cbd5e1", markersize = 24)
    scatter!(ax_torso, [pos_rl], color = "#cbd5e1", markersize = 20)

    # Inner IEC standard colored snap buttons (RA: Red, LA: Yellow, LL: Green, RL: Black)
    scatter!(ax_torso, [pos_ra, pos_la, pos_ll], color = ["#ef4444", "#eab308", "#22c55e"], markersize = 14)
    scatter!(ax_torso, [pos_rl], color = "#1e293b", markersize = 12)

    # Dynamic Polarity Labels (+ / – / GND) based on active lead
    lbl_ra_pol = lift(selected_lead) do s
        s in (:lead1, :lead2, :all3) ? "RA (–)" : "RA"
    end
    lbl_la_pol = lift(selected_lead) do s
        if s == :lead1
            "LA (+)"
        elseif s == :lead3
            "LA (–)"
        elseif s == :all3
            "LA (+/–)"
        else
            "LA"
        end
    end
    lbl_ll_pol = lift(selected_lead) do s
        s in (:lead2, :lead3, :all3) ? "LL (+)" : "LL"
    end

    text!(ax_torso, pos_ra[1] - 0.28, pos_ra[2], text = lbl_ra_pol, align = (:right, :center), font = :bold, fontsize = 11, color = "#991b1b")
    text!(ax_torso, pos_la[1] + 0.28, pos_la[2], text = lbl_la_pol, align = (:left, :center), font = :bold, fontsize = 11, color = "#854d0e")
    text!(ax_torso, pos_ll[1] + 0.28, pos_ll[2], text = lbl_ll_pol, align = (:left, :center), font = :bold, fontsize = 11, color = "#166534")
    text!(ax_torso, pos_rl[1] - 0.25, pos_rl[2], text = "RL (GND)", align = (:right, :center), font = :bold, fontsize = 9.5, color = "#475569")

    # Standard clinical monitor sweep window (5.0 s)
    t_window = 5.0

    # ── Right Column [1, 2]: ECG Strip ───────────────────────────────────────
    ax_ecg = Axis(
        fig[1, 2],
        title = "ECG Recording",
        xlabel = "Time [s]",
        ylabel = "Potential (a.u.)",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = :white,
        xgridcolor = "#fecdd3",
        xgridwidth = 1.0,
        ygridcolor = "#fecdd3",
        ygridwidth = 1.0,
        xminorgridcolor = "#fff1f2",
        xminorgridwidth = 0.5,
        yminorgridcolor = "#fff1f2",
        yminorgridwidth = 0.5,
        xminorgridvisible = true,
        yminorgridvisible = true,
    )
    xlims!(ax_ecg, 0.0, t_window)
    ylims!(ax_ecg, -0.75, 1.75)

    hlines!(ax_ecg, [0.0], color = "#f43f5e", linestyle = :dash, linewidth = 1.2)

    # ── ECG Signal Synthesis Pipeline ────────────────────────────────────────
    function generate_ecg_traces(alpha_val, bpm_val, n_lvl)
        T = 60.0 / bpm_val
        n_samples = 1500
        t_arr = collect(range(0.0, t_window, length = n_samples))

        pts_I   = Vector{Point2f}(undef, n_samples)
        pts_II  = Vector{Point2f}(undef, n_samples)
        pts_III = Vector{Point2f}(undef, n_samples)

        n_scale = (n_lvl / 100.0) * 0.30

        for i in 1:n_samples
            t = t_arr[i]
            t_in = mod(t, T)
            px, py = compute_cardiac_dipole(t_in, bpm_val, alpha_val)

            # Projections
            v1 = px * 1.0 + py * 0.0
            v2 = px * cosd(60.0) - py * sind(60.0)
            v3 = px * cosd(120.0) - py * sind(120.0)

            # Noise
            n = n_scale * (0.6 * sin(i * 12.9898) + 0.35 * sin(2π * 50.0 * t))

            pts_I[i]   = Point2f(t, v1 + n)
            pts_II[i]  = Point2f(t, v2 + n)
            pts_III[i] = Point2f(t, v3 + n)
        end

        return (; pts_I, pts_II, pts_III)
    end

    ecg_data = lift(axis_alpha_deg, heart_rate_bpm, noise_level) do alpha, bpm, n_lvl
        generate_ecg_traces(alpha, bpm, n_lvl)
    end

    # Active single lead signal
    active_ecg_pts = lift(ecg_data, selected_lead) do d, s
        if s == :lead1
            d.pts_I
        elseif s == :lead2
            d.pts_II
        elseif s == :lead3
            d.pts_III
        else
            d.pts_II
        end
    end

    active_ecg_color = lift(selected_lead) do s
        if s == :lead1
            "#0284c7"
        elseif s == :lead2
            "#e11d48"
        elseif s == :lead3
            "#10b981"
        else
            "#e11d48"
        end
    end

    # Single active lead line
    lines!(
        ax_ecg,
        active_ecg_pts,
        color = active_ecg_color,
        linewidth = 2.4,
        visible = lift(s -> s != :all3, selected_lead),
    )

    # Multi-lead superimposed lines for "All 3 (I, II, III)" mode
    lines!(
        ax_ecg,
        lift(d -> d.pts_I, ecg_data),
        color = "#0284c7",
        linewidth = 2.0,
        label = "Lead I (0°)",
        visible = lift(s -> s == :all3, selected_lead),
    )
    lines!(
        ax_ecg,
        lift(d -> d.pts_II, ecg_data),
        color = "#e11d48",
        linewidth = 2.4,
        label = "Lead II (+60°)",
        visible = lift(s -> s == :all3, selected_lead),
    )
    lines!(
        ax_ecg,
        lift(d -> d.pts_III, ecg_data),
        color = "#10b981",
        linewidth = 2.0,
        label = "Lead III (+120°)",
        visible = lift(s -> s == :all3, selected_lead),
    )

    leg_ecg = axislegend(ax_ecg, position = :rt, labelsize = 12, backgroundcolor = RGBAf(1, 1, 1, 0.88))
    leg_ecg.blockscene.visible[] = false
    on(selected_lead) do s
        leg_ecg.blockscene.visible[] = (s == :all3)
    end

    # Sweeping time cursor during playback
    cursor_x = lift(playback_time) do pt
        pt < 0.0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -0.75), Point2f(pt, 1.75)]
    end
    lines!(ax_ecg, cursor_x, color = "#e11d48", linewidth = 2.0, linestyle = :dash)

    # Wave markers (P, Q, R, S, T) on first beat
    t_P = lift(heart_rate_bpm) do bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        0.12 * k
    end
    t_Q = lift(heart_rate_bpm) do bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        0.22 * k
    end
    t_R = lift(heart_rate_bpm) do bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        0.25 * k
    end
    t_S = lift(heart_rate_bpm) do bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        0.28 * k
    end
    t_T = lift(heart_rate_bpm) do bpm
        T = 60.0 / bpm
        k = bpm > 72.0 ? (T / (60.0 / 72.0)) : 1.0
        qt_scale = bpm > 72.0 ? k : min(1.35, sqrt(T / (60.0 / 72.0)))
        (0.25 + 0.18 * qt_scale) * (bpm > 72.0 ? k : 1.0)
    end

    # Wave marker label styling: consistent crisp dark tone, responsive font sizing
    wave_font = lift(f -> f * 1.20, label_font)
    wave_col  = "#1e293b"

    text!(ax_ecg, t_P, 0.35, text = "P", align = (:center, :bottom), font = :bold, fontsize = wave_font, color = wave_col)
    text!(ax_ecg, t_Q, -0.28, text = "Q", align = (:center, :top),    font = :bold, fontsize = wave_font, color = wave_col)
    text!(ax_ecg, t_R, 1.52, text = "R", align = (:center, :bottom), font = :bold, fontsize = lift(f -> f + 1.5, wave_font), color = wave_col)
    text!(ax_ecg, t_S, -0.45, text = "S", align = (:center, :top),    font = :bold, fontsize = wave_font, color = wave_col)
    text!(ax_ecg, t_T, 0.55, text = "T", align = (:center, :bottom), font = :bold, fontsize = wave_font, color = wave_col)

    # ── Bottom Control Panel [2, 1:2]: Sliders, Buttons & Toggles ────────────
    ctrl_layout = fig[2, 1:2] = GridLayout(tellwidth = false, halign = :center)

    # Sliders Sub-Grid (Row 1)
    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(slider_sub[1, 1], text = "Heart Axis (α):", font = :bold, fontsize = label_font)
    sl_axis = Slider(slider_sub[1, 2], range = -45.0:1.0:120.0, startvalue = 60.0, width = 110)
    lbl_axis = Label(slider_sub[1, 3], text = lift(v -> "$(round(v, digits=1))°", sl_axis.value), width = 60, fontsize = label_font)

    Label(slider_sub[1, 5], text = "HR (BPM):", font = :bold, fontsize = label_font)
    sl_hr = Slider(slider_sub[1, 6], range = 10.0:1.0:220.0, startvalue = 72.0, width = 110)
    lbl_hr = Label(slider_sub[1, 7], text = lift(v -> "$(round(Int, v)) BPM", sl_hr.value), width = 75, fontsize = label_font)

    Label(slider_sub[1, 9], text = "Noise:", font = :bold, fontsize = label_font)
    sl_noise = Slider(slider_sub[1, 10], range = 0.0:1.0:200.0, startvalue = 15.0, width = 120)

    colgap!(slider_sub, 8)
    colsize!(slider_sub, 4, Fixed(50))
    colsize!(slider_sub, 8, Fixed(50))

    rowgap!(ctrl_layout, 12)

    # Buttons & Checkboxes Sub-Grid (Row 2)
    btn_sub = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)

    btn_play = Button(btn_sub[1, 1], label = "Play Heartbeat", buttoncolor = "#dcfce7", fontsize = btn_font, width = 160, height = 36)

    Label(btn_sub[1, 3], text = "Lead:", font = :bold, fontsize = label_font)
    btn_l1   = Button(btn_sub[1, 4], label = "Lead I", fontsize = btn_font, width = 85, height = 36)
    btn_l2   = Button(btn_sub[1, 5], label = "Lead II (Std)", buttoncolor = "#fee2e2", fontsize = btn_font, width = 115, height = 36)
    btn_l3   = Button(btn_sub[1, 6], label = "Lead III", fontsize = btn_font, width = 90, height = 36)
    btn_all  = Button(btn_sub[1, 7], label = "All 3 (I, II, III)", fontsize = btn_font, width = 140, height = 36)

    chk_vcg = Checkbox(btn_sub[1, 9], checked = true)
    Label(btn_sub[1, 10], text = "Show Vectorcardiogram", fontsize = label_font)

    chk_proj = Checkbox(btn_sub[1, 12], checked = true)
    Label(btn_sub[1, 13], text = "Show Projection", fontsize = label_font)

    colgap!(btn_sub, 8)
    colsize!(btn_sub, 2, Fixed(24))
    colsize!(btn_sub, 8, Fixed(28))
    colsize!(btn_sub, 11, Fixed(20))

    # ── Interactive Callback Connections ─────────────────────────────────────
    on(sl_axis.value) do val
        axis_alpha_deg[] = val
    end

    on(sl_hr.value) do val
        heart_rate_bpm[] = val
    end

    on(sl_noise.value) do val
        noise_level[] = val
    end

    on(chk_vcg.checked) do val
        show_vcg_loop[] = val
    end

    on(chk_proj.checked) do val
        show_projection[] = val
    end

    # Lead Selection Button Callbacks with Visual Highlighting
    function set_lead_button_colors(active_sym)
        btn_l1.buttoncolor[]  = active_sym == :lead1 ? "#bae6fd" : "#f1f5f9"
        btn_l2.buttoncolor[]  = active_sym == :lead2 ? "#fecdd3" : "#f1f5f9"
        btn_l3.buttoncolor[]  = active_sym == :lead3 ? "#a7f3d0" : "#f1f5f9"
        btn_all.buttoncolor[] = active_sym == :all3  ? "#fbcfe8" : "#f1f5f9"
    end
    set_lead_button_colors(:lead2)

    on(btn_l1.clicks) do _
        selected_lead[] = :lead1
        set_lead_button_colors(:lead1)
    end
    on(btn_l2.clicks) do _
        selected_lead[] = :lead2
        set_lead_button_colors(:lead2)
    end
    on(btn_l3.clicks) do _
        selected_lead[] = :lead3
        set_lead_button_colors(:lead3)
    end
    on(btn_all.clicks) do _
        selected_lead[] = :all3
        set_lead_button_colors(:all3)
    end

    # Play/Pause Heartbeat Animation Loop
    on(btn_play.clicks) do _
        is_playing[] = !is_playing[]
        if is_playing[]
            btn_play.label[] = "Pause Heartbeat"
            btn_play.buttoncolor[] = "#fee2e2"

            @async begin
                fps = 30.0
                dt = 1.0 / fps
                curr_t = playback_time[] < 0.0 ? 0.0 : playback_time[]

                try
                    while is_playing[]
                        curr_t += dt
                        if curr_t > t_window
                            curr_t = 0.0
                        end
                        playback_time[] = curr_t
                        sleep(dt)
                    end
                finally
                    if !is_playing[]
                        btn_play.label[] = "Play Heartbeat"
                        btn_play.buttoncolor[] = "#dcfce7"
                    end
                end
            end
        else
            btn_play.label[] = "Play Heartbeat"
            btn_play.buttoncolor[] = "#dcfce7"
        end
    end

    # Relative Grid Layout Sizing
    rowsize!(fig.layout, 1, Relative(0.80))
    rowsize!(fig.layout, 2, Relative(0.20))
    colsize!(fig.layout, 1, Relative(0.53))
    colsize!(fig.layout, 2, Relative(0.47))

    display(fig)

    return fig, ax_torso, ax_ecg
end
