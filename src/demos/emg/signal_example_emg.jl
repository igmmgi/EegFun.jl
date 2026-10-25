"""
    signal_example_emg()

Interactive Surface Electromyography (sEMG) Simulator —
Muscle Activation & Signal Processing Pipeline Demo (Raw → Rectify → Filter → Integrate).

Creates an interactive psychophysiology demonstration illustrating the cellular origin of
surface electromyography (sEMG) and its canonical 4-stage signal processing cascade:
**Raw → Rectify → Filter → Integrate**.

## Educational Background

### 1. Biophysical Source of Surface EMG
- **Motor Unit Pool (100s of Units)**: A skeletal muscle is innervated by a motor neuron pool
  containing hundreds of alpha motor neurons in the spinal cord. Each motor neuron and all the
  individual muscle fibers it innervates constitute a single Motor Unit (MU). The schematic
  depicts this dense motor neuron pool with representative units (`MU 1` through `MU n`) and
  clustered somas.
- **Motor Nerve & Asynchronous Rate Coding**: Axons exit the spinal cord and travel together in a
  peripheral motor nerve trunk before arborizing into neuromuscular junctions across the muscle.
  During voluntary contraction, hundreds of motor units fire action potential trains
  asynchronously at frequencies scaling with effort (rate coding).
- **Surface Interference Pattern**: Macro bipolar electrodes placed over the muscle record
  the collective volume-conducted superposition of thousands of depolarization waves. Because
  positive and negative phases overlap and balance out across fibers and space, the raw surface
  EMG is a zero-mean, stochastic-like interference pattern (10–500 Hz). Individual motor units
  cannot be isolated with macro surface electrodes.

### 2. The 4-Stage Signal Processing Pipeline
1. **Raw Surface EMG [µV]**:
   Bipolar difference potential \$V(t) = V_{(+)} - V_{(-)}\$. Because positive and negative phases
   balance out, the arithmetic mean of raw EMG is approximately 0 µV regardless of contraction force.
2. **Rectify [µV]**:
   Calculates the full-wave rectified absolute value \$|V(t)|\$. Converts all negative excursions
   into positive deflections, preserving signal energy for subsequent averaging.
3. **Filter [µV]**:
   Extracts the smoothed linear envelope / Root Mean Square (RMS) over a sliding time window \$\\tau\$:
   ```math
   \\text{RMS}(t) = \\sqrt{\\frac{1}{\\tau} \\int_{t - \\tau/2}^{t + \\tau/2} V^2(s) \\, ds}
   ```
   - **Window Size Trade-off**: Short windows (e.g. 25 ms) retain rapid fluctuations but remain noisy;
     long windows (e.g. 200 ms) produce a smooth envelope but introduce temporal lag / smoothing.
   - **Electromechanical Delay (EMD)**: Physiological delay (~30–60 ms) between electrical onset
     and mechanical tension development.
4. **Integrate [µV·s] (iEMG)**:
   Calculates the cumulative integral of the rectified signal over time:
   ```math
   \\text{iEMG}(t) = \\int_0^t |V(s)| \\, ds
   ```
   - **Teaching Note**: If one integrates raw EMG, the positive and negative phases cancel out to ~0 µV·s.
     Rectification is the crucial non-linear prerequisite that allows integration to measure cumulative
     muscular work and total neural drive over an epoch.

## Controls

| Control | Description |
|:---|:---|
| **Force Slider** | Commanded contraction force (0–100% Maximum Voluntary Contraction / MVC). |
| **Duration Slider** | Contraction burst hold duration (200–1400 ms). |
| **RMS Window Slider** | Sliding RMS smoothing window (20–250 ms) to demonstrate filtering trade-offs. |
| **Noise Slider** | Injects baseline physiological & amplifier noise (0–100 µV). |
| **▶ Play Contraction** | Plays an animated real-time contraction sequence with live cursor sweep & active muscle bulge. |

## Interactive Backend (GLMakie)

This interactive demo requires **GLMakie** for real-time reactivity, dragging sliders,
and smooth muscle contraction animations:

```julia
using GLMakie
using EegFun

signal_example_emg()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI.
- `ax_anatomy::Axis`: The anatomical muscle, motor units & bipolar amplifier schematic axis.
- `ax_raw::Axis`: The raw surface EMG interference pattern axis.
- `ax_rect::Axis`: The full-wave rectified EMG axis.
- `ax_filt::Axis`: The linear envelope / filtered EMG axis.
- `ax_integ::Axis`: The cumulative integrated EMG (iEMG) axis.
"""
function signal_example_emg()
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

    _set_window_title("Interactive Surface EMG Simulator — Raw → Rectify → Filter → Integrate")

    fig = Figure(
        size = (1440, 920),
        title = "Interactive Surface EMG Simulator — Raw → Rectify → Filter → Integrate",
        backgroundcolor = :white,
    )

    # ── Adaptive Font and UI Sizing ──────────────────────────────────────────
    function setup_adaptive_sizing(fig)
        title_font = Observable(16)
        label_font = Observable(12.5)
        tick_font  = Observable(10.5)
        btn_font   = Observable(12.5)

        on(fig.scene.viewport) do area
            scale        = area.widths[1] / 1440
            title_font[] = max(13, round(Int, 16 * scale))
            label_font[] = max(11, round(Int, 12.5 * scale))
            tick_font[]  = max(10, round(Int, 10.5 * scale))
            btn_font[]   = max(11, round(Int, 12.5 * scale))
        end

        return title_font, label_font, tick_font, btn_font
    end

    title_font, label_font, tick_font, btn_font = setup_adaptive_sizing(fig)

    # ── Observables ──────────────────────────────────────────────────────────
    target_force       = Observable(45.0)  # Commanded % MVC (0–100%)
    current_disp_force = Observable(45.0)  # Currently rendered % MVC (for animation)
    duration_ms        = Observable(700.0) # Contraction duration in ms
    rms_window_ms      = Observable(75.0)  # Moving RMS window width in ms
    noise_amp          = Observable(4.0)   # Baseline noise amplitude in µV
    playback_time      = Observable(-1.0)  # < 0: static, >= 0: time cursor position [s]

    # ── Geometry & Anatomical Constants ──────────────────────────────────────
    color_mu_active   = "#1e293b" # dark slate / black
    color_mu_inactive = "#94a3b8" # muted grey

    n_pts = 48
    theta_circ = range(0, 2π, length = n_pts)

    # ── Left Column [1, 1]: Anatomy, Motor Units & Circuit Schematic ─────────
    ax_anatomy = Axis(fig[1, 1], title = "Motor Unit Pool, Muscle & Bipolar Recording Setup", titlesize = title_font, aspect = DataAspect())
    hidedecorations!(ax_anatomy)
    hidespines!(ax_anatomy)
    xlims!(ax_anatomy, -5.90, 4.10)
    ylims!(ax_anatomy, -3.35, 2.95)

    # ── Tendons & Bone Attachments ───────────────────────────────────────────
    y_center = 0.50
    h_muscle = 0.86

    # Left Tendon (Origin)
    poly!(
        ax_anatomy,
        [Point2f(-3.5, y_center + 0.30), Point2f(-2.2, y_center + 0.12), Point2f(-2.2, y_center - 0.12), Point2f(-3.5, y_center - 0.30)],
        color = "#cbd5e1",
        strokecolor = "#64748b",
        strokewidth = 1.5,
    )
    text!(
        ax_anatomy,
        -2.85,
        y_center - 0.38,
        text = "Tendon (Origin)",
        align = (:center, :top),
        font = :bold,
        fontsize = 9.0,
        color = "#475569",
    )

    # Right Tendon (Insertion)
    poly!(
        ax_anatomy,
        [Point2f(2.2, y_center + 0.12), Point2f(3.5, y_center + 0.30), Point2f(3.5, y_center - 0.30), Point2f(2.2, y_center - 0.12)],
        color = "#cbd5e1",
        strokecolor = "#64748b",
        strokewidth = 1.5,
    )
    text!(
        ax_anatomy,
        2.85,
        y_center - 0.38,
        text = "Tendon (Insertion)",
        align = (:center, :top),
        font = :bold,
        fontsize = 9.0,
        color = "#475569",
    )

    # ── Dynamic Muscle (Contour Bulges with Force) ───────────────────────────
    function make_muscle_poly(force)
        bulge = 1.0 + 0.30 * (force / 100.0)
        xs_top = range(-2.2, 2.2, length = 32)
        pts_top = [Point2f(x, y_center + (h_muscle * bulge) * (1.0 - (x / 2.2)^2)) for x in xs_top]
        xs_bot = range(2.2, -2.2, length = 32)
        pts_bot = [Point2f(x, y_center - (h_muscle * bulge) * (1.0 - (x / 2.2)^2)) for x in xs_bot]
        return [pts_top; pts_bot]
    end

    muscle_poly = lift(make_muscle_poly, current_disp_force)
    poly!(ax_anatomy, muscle_poly, color = "#991b1b", strokecolor = "#450a0a", strokewidth = 2.0)

    # ── Dense Muscle Fiber Architecture (Conveys 1000s of fibers in overlapping MUs) ──
    striation_factors = [-0.75, -0.55, -0.35, -0.18, 0.0, 0.18, 0.35, 0.55, 0.75]
    for s_factor in striation_factors
        s_pts = lift(current_disp_force) do force
            bulge = 1.0 + 0.30 * (force / 100.0)
            xs = range(-2.1, 2.1, length = 40)
            [Point2f(x, y_center + s_factor * (0.80 * bulge) * (1.0 - (x / 2.2)^2)) for x in xs]
        end
        lines!(ax_anatomy, s_pts, color = RGBAf(0.3, 0.04, 0.04, 0.45), linewidth = 1.5, linestyle = :dot)
    end

    text!(
        ax_anatomy,
        0.5,
        y_center,
        text = "Muscle\n(1000s of fibers in\noverlapping MUs)",
        align = (:center, :center),
        font = :bold,
        fontsize = 8.0,
        color = "#fef2f2",
    )

    # ── Spinal Motor Neuron Pool (Dense cluster conveying 100s of units) ────
    poly!(ax_anatomy, Rect2f(-5.70, -0.75, 1.25, 3.10), color = RGBAf(0.94, 0.96, 0.98, 0.85), strokecolor = "#cbd5e1", strokewidth = 1.5)
    text!(
        ax_anatomy,
        -5.07,
        2.40,
        text = "Motor Neuron Pool\n(~100–1000+ MUs)",
        align = (:center, :bottom),
        font = :bold,
        fontsize = 8.5,
        color = "#0f172a",
    )

    # Clustered background somas representing the larger pool population
    bg_somas = [
        Point2f(-5.48, 2.10),
        Point2f(-5.28, 1.85),
        Point2f(-5.48, 1.60),
        Point2f(-5.30, 1.35),
        Point2f(-5.50, 1.05),
        Point2f(-5.30, 0.80),
        Point2f(-5.50, 0.50),
        Point2f(-5.30, 0.20),
        Point2f(-5.48, -0.10),
        Point2f(-5.30, -0.35),
        Point2f(-5.48, -0.58),
    ]
    bg_soma_r = 0.10
    for pos in bg_somas
        poly!(
            ax_anatomy,
            [pos .+ Point2f(bg_soma_r * cos(t), bg_soma_r * sin(t)) for t in theta_circ],
            color = lift(f -> f > 0.0 ? "#475569" : "#cbd5e1", current_disp_force),
            strokecolor = "#334155",
            strokewidth = 1.0,
        )
    end

    # Representative foreground somas with designated axon outputs
    soma1_pos = Point2f(-4.88, 2.05) # MU 1
    soma2_pos = Point2f(-4.88, 1.55) # MU 2
    soma3_pos = Point2f(-4.88, 1.05) # MU 3
    soma4_pos = Point2f(-4.88, 0.55) # MU 4
    somaN_pos = Point2f(-4.88, -0.25) # MU n
    rep_somas = (soma1_pos, soma2_pos, soma3_pos, soma4_pos, somaN_pos)

    r_soma = 0.17

    # Active glow halos around representative somas scaling with force
    for pos in rep_somas
        poly!(
            ax_anatomy,
            lift(f -> [pos .+ Point2f((r_soma + 0.12) * cos(t), (r_soma + 0.12) * sin(t)) for t in theta_circ], current_disp_force),
            color = lift(f -> f > 0.0 ? RGBAf(0.12, 0.16, 0.23, 0.15 + 0.30 * (f / 100.0)) : RGBAf(0, 0, 0, 0), current_disp_force),
        )
    end

    # Representative soma bodies (black/grey)
    for pos in rep_somas
        poly!(
            ax_anatomy,
            [pos .+ Point2f(r_soma * cos(t), r_soma * sin(t)) for t in theta_circ],
            color = lift(f -> f > 0.0 ? color_mu_active : color_mu_inactive, current_disp_force),
            strokecolor = "#334155",
            strokewidth = 1.5,
        )
    end

    # Motor unit pool labels
    text!(ax_anatomy, -4.50, soma1_pos[2], text = "MU 1", align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")
    text!(ax_anatomy, -4.50, soma2_pos[2], text = "MU 2", align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")
    text!(ax_anatomy, -4.50, soma3_pos[2], text = "MU 3", align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")
    text!(ax_anatomy, -4.50, soma4_pos[2], text = "MU 4", align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")
    text!(ax_anatomy, -4.50, 0.15, text = "⋮", align = (:left, :center), font = :bold, fontsize = 12.0, color = "#64748b")
    text!(ax_anatomy, -4.50, somaN_pos[2], text = "MU n", align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")

    # ── Peripheral Motor Nerve Conduit (Axon Bundle) ─────────────────────────
    nerve_poly = [
        Point2f(-4.25, 2.00),
        Point2f(-3.50, 1.55),
        Point2f(-2.65, 1.20),
        Point2f(-2.65, 0.15),
        Point2f(-3.50, 0.35),
        Point2f(-4.25, -0.45),
    ]
    poly!(ax_anatomy, nerve_poly, color = RGBAf(0.92, 0.95, 0.98, 0.65), strokecolor = "#94a3b8", strokewidth = 1.2, linestyle = :dash)
    text!(
        ax_anatomy,
        -3.50,
        1.72,
        text = "Motor Nerve\n(Axon Bundle)",
        align = (:center, :bottom),
        font = :bold,
        fontsize = 8.0,
        color = "#475569",
    )

    # Axons routing through the nerve bundle and arborizing into the muscle
    axon1 = [soma1_pos, Point2f(-3.7, 1.45), Point2f(-2.6, 1.10), Point2f(-1.3, 1.00)]
    axon2 = [soma2_pos, Point2f(-3.7, 1.20), Point2f(-2.6, 0.90), Point2f(-0.9, 0.78)]
    axon3 = [soma3_pos, Point2f(-3.7, 0.95), Point2f(-2.6, 0.70), Point2f(-1.5, 0.58)]
    axon4 = [soma4_pos, Point2f(-3.7, 0.70), Point2f(-2.6, 0.50), Point2f(-1.1, 0.38)]
    axon5 = [Point2f(-5.30, 0.20), Point2f(-3.7, 0.45), Point2f(-2.6, 0.30), Point2f(-1.4, 0.18)]
    axon6 = [somaN_pos, Point2f(-3.7, 0.20), Point2f(-2.6, 0.10), Point2f(-1.0, -0.05)]

    axon_col = lift(f -> f > 0.0 ? "#334155" : "#cbd5e1", current_disp_force)
    lines!(ax_anatomy, axon1, color = axon_col, linewidth = 2.0, linestyle = :dot)
    lines!(ax_anatomy, axon2, color = axon_col, linewidth = 2.0, linestyle = :dot)
    lines!(ax_anatomy, axon3, color = axon_col, linewidth = 2.0, linestyle = :dot)
    lines!(ax_anatomy, axon4, color = axon_col, linewidth = 2.0, linestyle = :dot)
    lines!(ax_anatomy, axon5, color = axon_col, linewidth = 2.0, linestyle = :dot)
    lines!(ax_anatomy, axon6, color = axon_col, linewidth = 2.0, linestyle = :dot)

    # Neuromuscular Junctions (Motor Endplates)
    endplate_col = lift(f -> f > 0.0 ? color_mu_active : color_mu_inactive, current_disp_force)
    endplate_pts =
        [Point2f(-1.3, 1.00), Point2f(-0.9, 0.78), Point2f(-1.5, 0.58), Point2f(-1.1, 0.38), Point2f(-1.4, 0.18), Point2f(-1.0, -0.05)]
    scatter!(ax_anatomy, endplate_pts, marker = :utriangle, markersize = 11, color = endplate_col)
    text!(
        ax_anatomy,
        -1.65,
        -0.28,
        text = "Neuromuscular Junctions\n(Innervation Zone)",
        align = (:right, :center),
        font = :bold,
        fontsize = 8.0,
        color = "#475569",
    )

    # ── Bipolar Surface Electrodes & Difference Amplifier Circuit ────────────
    y_skin = 2.15
    lines!(ax_anatomy, [-2.2, 2.2], [y_skin, y_skin], color = "#fcd34d", linewidth = 3.0)
    text!(ax_anatomy, 2.30, y_skin, text = "Skin Surface", align = (:left, :center), font = :italic, fontsize = 9.0, color = "#78350f")

    # Bipolar Electrodes placed along fiber axis
    poly!(ax_anatomy, Rect2f(-0.95, y_skin, 0.40, 0.28), color = "#2563eb", strokewidth = 1.2, strokecolor = "#1e293b")
    text!(ax_anatomy, -0.75, y_skin + 0.38, text = "[-]", align = (:center, :center), font = :bold, fontsize = 11.5, color = "#2563eb")

    poly!(ax_anatomy, Rect2f(0.55, y_skin, 0.40, 0.28), color = "#dc2626", strokewidth = 1.2, strokecolor = "#1e293b")
    text!(ax_anatomy, 0.75, y_skin + 0.38, text = "[+]", align = (:center, :center), font = :bold, fontsize = 11.5, color = "#dc2626")

    # Difference Amplifier (Positioned with generous clearance below muscle at y = -2.10)
    y_amp = -2.10
    amp_tri = [Point2f(-0.70, y_amp + 0.45), Point2f(-0.70, y_amp - 0.45), Point2f(0.70, y_amp)]
    poly!(ax_anatomy, amp_tri, color = "#f8fafc", strokewidth = 2.0, strokecolor = "#1e293b")
    text!(ax_anatomy, -0.48, y_amp + 0.22, text = "+", align = (:center, :center), font = :bold, fontsize = 17, color = "#dc2626")
    text!(ax_anatomy, -0.48, y_amp - 0.22, text = "–", align = (:center, :center), font = :bold, fontsize = 19, color = "#2563eb")
    text!(
        ax_anatomy,
        0.0,
        y_amp - 0.58,
        text = "Difference Amplifier",
        align = (:center, :top),
        font = :bold,
        fontsize = 11.0,
        color = "#1e293b",
    )
    text!(
        ax_anatomy,
        0.0,
        y_amp - 0.78,
        text = "Vout = V(+) – V(-)",
        align = (:center, :top),
        font = :italic,
        fontsize = 8.5,
        color = "#64748b",
    )

    # Lead wires routing cleanly around muscle with ample clearance
    wire_pos = [
        Point2f(0.75, y_skin),
        Point2f(1.60, y_skin),
        Point2f(1.60, -1.35),
        Point2f(-1.05, -1.35),
        Point2f(-1.05, y_amp + 0.22),
        Point2f(-0.70, y_amp + 0.22),
    ]
    lines!(ax_anatomy, wire_pos, color = "#dc2626", linewidth = 2.0)

    wire_neg = [Point2f(-0.75, y_skin), Point2f(-1.60, y_skin), Point2f(-1.60, y_amp - 0.22), Point2f(-0.70, y_amp - 0.22)]
    lines!(ax_anatomy, wire_neg, color = "#2563eb", linewidth = 2.0)

    # Output lead
    lines!(ax_anatomy, [Point2f(0.70, y_amp), Point2f(1.50, y_amp)], color = "#1e293b", linewidth = 2.2)
    scatter!(
        ax_anatomy,
        [Point2f(1.55, y_amp)],
        marker = :circle,
        markersize = 9,
        color = :white,
        strokecolor = "#1e293b",
        strokewidth = 2.0,
    )
    text!(ax_anatomy, 1.70, y_amp, text = "Vout", align = (:left, :center), font = :bold, fontsize = 12, color = "#0f172a")

    # ── Left Column [2, 1]: 1x4 Slider Grid & Play Button ────────────────────
    ctrl_layout = fig[2, 1] = GridLayout(tellwidth = false, halign = :center)

    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    # Row 1: Force Slider
    Label(slider_sub[1, 1], text = "Force:", font = :bold, fontsize = label_font, halign = :right)
    sl_force = Slider(slider_sub[1, 2], range = 0.0:1.0:100.0, startvalue = 45.0, width = 210)
    lbl_force = Label(
        slider_sub[1, 3],
        text = lift(v -> "$(round(Int, v))% MVC", sl_force.value),
        width = 76,
        halign = :left,
        fontsize = label_font,
    )

    # Row 2: Duration Slider
    Label(slider_sub[2, 1], text = "Duration:", font = :bold, fontsize = label_font, halign = :right)
    sl_dur = Slider(slider_sub[2, 2], range = 200.0:50.0:1400.0, startvalue = 700.0, width = 210)
    lbl_dur =
        Label(slider_sub[2, 3], text = lift(v -> "$(round(Int, v)) ms", sl_dur.value), width = 76, halign = :left, fontsize = label_font)

    # Row 3: RMS Window Slider
    Label(slider_sub[3, 1], text = "RMS Win:", font = :bold, fontsize = label_font, halign = :right)
    sl_rms = Slider(slider_sub[3, 2], range = 20.0:5.0:250.0, startvalue = 75.0, width = 210)
    lbl_rms =
        Label(slider_sub[3, 3], text = lift(v -> "$(round(Int, v)) ms", sl_rms.value), width = 76, halign = :left, fontsize = label_font)

    # Row 4: Noise Slider (0–100 µV)
    Label(slider_sub[4, 1], text = "Noise:", font = :bold, fontsize = label_font, halign = :right)
    sl_noise = Slider(slider_sub[4, 2], range = 0.0:1.0:100.0, startvalue = 4.0, width = 210)
    lbl_noise =
        Label(slider_sub[4, 3], text = lift(v -> "$(round(Int, v)) µV", sl_noise.value), width = 76, halign = :left, fontsize = label_font)

    colgap!(slider_sub, 10)
    rowgap!(slider_sub, 4)

    # Row 2: Play Contraction Button (Centered below sliders)
    btn_sub = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    btn_play = Button(btn_sub[1, 1], label = "▶ Play Contraction", buttoncolor = "#dcfce7", fontsize = btn_font, width = 200, height = 34)

    rowgap!(ctrl_layout, 8)

    # ── Synthetic EMG & Processing Pipeline Calculation ──────────────────────
    fs = 1000.0
    dt = 1.0 / fs
    t_vec = collect(range(0.0, 2.0, step = dt))
    N_pts = length(t_vec)

    # Deterministic background physiological noise
    base_noise =
        [0.65 * sin(2π * 23.4 * tt) * cos(2π * 61.2 * tt) + 0.35 * sin(2π * 117.8 * tt) + 0.25 * sin(2π * 243.1 * tt) for tt in t_vec]

    # Pre-computed pseudo-random jitter phase sequences for 24 asynchronous motor units
    coprimes = [
        31.0,
        37.0,
        41.0,
        43.0,
        47.0,
        53.0,
        59.0,
        61.0,
        67.0,
        71.0,
        73.0,
        79.0,
        83.0,
        89.0,
        97.0,
        101.0,
        103.0,
        107.0,
        109.0,
        113.0,
        127.0,
        131.0,
        137.0,
        139.0,
    ]
    N_units = length(coprimes)
    jitter_phases = [0.20 * sin.(2π .* (1:N_pts) ./ coprimes[idx] .+ (1.3 * idx)) for idx = 1:N_units]

    # Motor unit pool: Henneman's size principle
    # Units have progressive recruitment thresholds (from 2% up to 80% MVC).
    # Small units recruited first fire at moderate rates; larger units recruited at higher force
    # fire at higher rates and generate larger action potentials.
    fiber_pools = [
        (
            thresh = 2.0 + 78.0 * ((i - 1) / (N_units - 1))^1.2,
            f_base = 8.0 + 4.0 * (i / N_units),
            f_max  = 24.0 + 32.0 * (i / N_units),
            amp    = 30.0 + 220.0 * ((i - 1) / (N_units - 1))^1.4,
            sigma  = 0.0018 + 0.0018 * (i / N_units),
        ) for i = 1:N_units
    ]

    function compute_emg_signals(force_val, dur_val, rms_win_val, noise_val)
        t_on = 0.40
        t_off = t_on + dur_val / 1000.0

        force_t = zeros(Float64, N_pts)
        for i = 1:N_pts
            tt = t_vec[i]
            s_on = 1.0 / (1.0 + exp(-(tt - t_on) / 0.025))
            s_off = 1.0 / (1.0 + exp((tt - t_off) / 0.025))
            force_t[i] = force_val * s_on * s_off
        end

        raw_sig = zeros(Float64, N_pts)

        for u_idx = 1:N_units
            p = fiber_pools[u_idx]
            j_phase = jitter_phases[u_idx]
            phi = 0.0
            w_samples = round(Int, 5.0 * p.sigma * fs)

            for i = 2:N_pts
                F = force_t[i]
                if F >= p.thresh
                    drive = (F - p.thresh) / (100.0 - p.thresh)
                    freq = p.f_base + (p.f_max - p.f_base) * drive
                    phi += freq * dt
                    phi_eff = phi + j_phase[i]
                    if floor(phi_eff) > floor(phi_eff - freq * dt)
                        amp_unit = p.amp * (0.60 + 0.40 * sqrt(drive))
                        for k = (-w_samples):w_samples
                            idx = i + k
                            if 1 <= idx <= N_pts
                                tau = k * dt
                                u = tau / p.sigma
                                val = amp_unit * (-u + 0.6 * u^3) * exp(-0.5 * u^2)
                                raw_sig[idx] += val
                            end
                        end
                    end
                end
            end
        end

        raw_sig .+= (noise_val .* base_noise)
        rect_sig = abs.(raw_sig)

        # Running RMS via cumulative sum
        w_rms = max(3, round(Int, (rms_win_val / 1000.0) * fs))
        half_w = div(w_rms, 2)
        raw_sq = raw_sig .^ 2
        csum = [0.0; cumsum(raw_sq)]
        rms_sig = zeros(Float64, N_pts)
        for i = 1:N_pts
            i1 = max(1, i - half_w)
            i2 = min(N_pts, i + half_w)
            diff_sq = max(0.0, csum[i2+1] - csum[i1])
            rms_sig[i] = sqrt(diff_sq / (i2 - i1 + 1))
        end

        # Stage 4: Cumulative Integral (iEMG = ∫ |V| dt) in µV·s
        integ_sig = cumsum(rect_sig) .* dt
        max_integ = integ_sig[end]

        return (; raw = raw_sig, rect = rect_sig, rms = rms_sig, integ = integ_sig, max_integ = max_integ, force_curve = force_t)
    end

    signals = lift(target_force, duration_ms, rms_window_ms, noise_amp) do f, d, r, n
        compute_emg_signals(f, d, r, n)
    end

    # ── Right Column [1:2, 2]: 4-Stage Signal Processing Pipeline Stack ──────
    pipe_layout = fig[1:2, 2] = GridLayout()

    # Stage 1: Raw
    ax_raw = Axis(
        pipe_layout[1, 1],
        title = "1. Raw Surface EMG [µV]",
        ylabel = "Potential [µV]",
        titlesize = title_font,
        ylabelsize = label_font,
        yticklabelsize = tick_font,
        ylabelpadding = 16.0,
        yticklabelspace = 48.0,
    )
    xlims!(ax_raw, 0.0, 2.0)
    ylims!(ax_raw, -700, 700)
    hidexdecorations!(ax_raw, grid = true, ticks = true)
    hlines!(ax_raw, [0.0], color = :gray75, linestyle = :dash, linewidth = 1.0)

    # Stage 2: Rectify
    ax_rect = Axis(
        pipe_layout[2, 1],
        title = "2. Rectify [µV] — |V(t)|",
        ylabel = "|V| [µV]",
        titlesize = title_font,
        ylabelsize = label_font,
        yticklabelsize = tick_font,
        ylabelpadding = 16.0,
        yticklabelspace = 48.0,
    )
    xlims!(ax_rect, 0.0, 2.0)
    ylims!(ax_rect, 0, 700)
    hidexdecorations!(ax_rect, grid = true, ticks = true)

    # Stage 3: Filter (Linear Envelope / RMS)
    ax_filt = Axis(
        pipe_layout[3, 1],
        title = "3. Filter [µV] — Linear Envelope",
        ylabel = "Envelope [µV]",
        titlesize = title_font,
        ylabelsize = label_font,
        yticklabelsize = tick_font,
        ylabelpadding = 16.0,
        yticklabelspace = 48.0,
    )
    xlims!(ax_filt, 0.0, 2.0)
    ylims!(ax_filt, 0, 320)
    hidexdecorations!(ax_filt, grid = true, ticks = true)

    # Stage 4: Integrate (Cumulative iEMG = ∫ |V| dt)
    ax_integ = Axis(
        pipe_layout[4, 1],
        title = "4. Integrate [µV·s] — iEMG = ∫ |V| dt",
        xlabel = "Time [s]",
        ylabel = "iEMG [µV·s]",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        ylabelpadding = 16.0,
        yticklabelspace = 48.0,
    )
    xlims!(ax_integ, 0.0, 2.0)
    ylims!(ax_integ, 0, 240)

    # Render Stage 1 (Raw)
    lines!(ax_raw, t_vec, lift(s -> s.raw, signals), color = "#0f172a", linewidth = 1.3)

    # Render Stage 2 (Rectify)
    band!(ax_rect, t_vec, fill(0.0, N_pts), lift(s -> s.rect, signals), color = RGBAf(0.0, 0.0, 0.0, 0.08))
    lines!(ax_rect, t_vec, lift(s -> s.rect, signals), color = "#0f172a", linewidth = 1.3)

    # Render Stage 3 (Filter)
    band!(ax_filt, t_vec, fill(0.0, N_pts), lift(s -> s.rms, signals), color = RGBAf(0.0, 0.0, 0.0, 0.10))
    lines!(ax_filt, t_vec, lift(s -> s.rms, signals), color = "#0f172a", linewidth = 2.2, label = "Filtered Envelope")

    # Force Curve overlay on Filter axis
    force_scaled = lift(signals) do s
        max_rms = max(10.0, maximum(s.rms))
        return (s.force_curve ./ 100.0) .* (0.92 * max_rms)
    end
    line_force = lines!(ax_filt, t_vec, force_scaled, color = "#64748b", linestyle = :dash, linewidth = 2.0, label = "Muscle Force")
    axislegend(ax_filt, position = :lt, labelsize = 10.5, backgroundcolor = RGBAf(1, 1, 1, 0.85))

    # Render Stage 4 (Integrate)
    band!(ax_integ, t_vec, fill(0.0, N_pts), lift(s -> s.integ, signals), color = RGBAf(0.0, 0.0, 0.0, 0.08))
    lines!(ax_integ, t_vec, lift(s -> s.integ, signals), color = "#0f172a", linewidth = 2.4, label = "iEMG")

    # Live iEMG Value readout badge directly on Integrate plot (mirrors slide readout)
    iemg_text = lift(signals) do s
        "iEMG = $(round(s.max_integ, digits=2)) µV·s"
    end
    text!(
        ax_integ,
        0.35,
        lift(s -> min(200.0, max(35.0, 0.80 * s.max_integ)), signals),
        text = iemg_text,
        font = :bold,
        fontsize = 12.5,
        color = "#0f172a",
    )

    # Sweeping Time Cursors during Contraction Playback across all 4 axes (slate grey dash)
    cursor_x_raw = lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, -700), Point2f(pt, 700)], playback_time)
    cursor_x_rec = lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, 0), Point2f(pt, 700)], playback_time)
    cursor_x_flt = lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, 0), Point2f(pt, 320)], playback_time)
    cursor_x_int = lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, 0), Point2f(pt, 240)], playback_time)

    lines!(ax_raw, cursor_x_raw, color = "#475569", linewidth = 2.0, linestyle = :dash)
    lines!(ax_rect, cursor_x_rec, color = "#475569", linewidth = 2.0, linestyle = :dash)
    lines!(ax_filt, cursor_x_flt, color = "#475569", linewidth = 2.0, linestyle = :dash)
    lines!(ax_integ, cursor_x_int, color = "#475569", linewidth = 2.0, linestyle = :dash)

    # Link time x-axes across all 4 pipeline stages
    linkxaxes!(ax_raw, ax_rect, ax_filt, ax_integ)

    # ── Interactive Callback Connections ─────────────────────────────────────
    is_animating = Ref(false)

    on(sl_force.value) do val
        target_force[] = val
        if !is_animating[]
            current_disp_force[] = val
        end
    end

    on(sl_dur.value) do val
        duration_ms[] = val
    end

    on(sl_rms.value) do val
        rms_window_ms[] = val
    end

    on(sl_noise.value) do val
        noise_amp[] = val
    end

    # Play Contraction Animation Callback
    on(btn_play.clicks) do _
        is_animating[] && return
        is_animating[] = true

        @async begin
            try
                target = target_force[]
                dur = duration_ms[] / 1000.0
                t_on = 0.40
                t_off = t_on + dur
                t_total = 2.0
                n_frames = 65
                fps = 30.0
                dt_frame = t_total / n_frames

                for frame = 0:n_frames
                    t_curr = frame * dt_frame
                    playback_time[] = t_curr

                    # Instantaneous force profile at this moment
                    s_on = 1.0 / (1.0 + exp(-(t_curr - t_on) / 0.025))
                    s_off = 1.0 / (1.0 + exp((t_curr - t_off) / 0.025))
                    f_moment = target * s_on * s_off
                    current_disp_force[] = f_moment

                    sleep(1.0 / fps)
                end

                playback_time[] = -1.0
                current_disp_force[] = target
            finally
                is_animating[] = false
            end
        end
    end

    # Relative Grid Layout Sizing
    rowsize!(fig.layout, 1, Relative(0.75))
    rowsize!(fig.layout, 2, Relative(0.25))
    colsize!(fig.layout, 1, Relative(0.53))
    colsize!(fig.layout, 2, Relative(0.47))

    rowsize!(pipe_layout, 1, Relative(0.26))
    rowsize!(pipe_layout, 2, Relative(0.24))
    rowsize!(pipe_layout, 3, Relative(0.25))
    rowsize!(pipe_layout, 4, Relative(0.25))

    display(fig)

    return fig, ax_anatomy, ax_raw, ax_rect, ax_filt, ax_integ
end
