"""
    signal_example_emg()

Interactive Surface Electromyography (sEMG) simulator illustrating motor unit recruitment
(Henneman's size principle), bipolar interference recording, and the 4-stage processing cascade:
**Raw → Rectify → Filter (RMS Envelope) → Integrate (iEMG)**.

```julia
using GLMakie, EegFun
signal_example_emg()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax_anatomy::Axis`: The anatomical muscle, motor units & bipolar amplifier schematic axis
- `ax_raw::Axis`: The raw surface EMG interference pattern axis
- `ax_rect::Axis`: The full-wave rectified EMG axis
- `ax_filt::Axis`: The linear envelope / filtered EMG axis
- `ax_integ::Axis`: The cumulative integrated EMG (iEMG) axis
"""
function signal_example_emg()
    _check_interactive_demo_backend("Interactive Surface EMG demo") || return nothing
    _set_window_title("Interactive Surface EMG Simulator — Raw → Rectify → Filter → Integrate")

    fig = _eegfun_figure(
        size = (1440, 920),
        title = "Interactive Surface EMG Simulator — Raw → Rectify → Filter → Integrate",
        backgroundcolor = :white,
    )
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1440, (16, 13), (12.5, 11), (10.5, 10), (12.5, 11))

    target_force = Observable(45.0)
    current_disp_force = Observable(45.0)
    duration_ms = Observable(700.0)
    rms_window_ms = Observable(75.0)
    noise_amp = Observable(4.0)
    playback_time = Observable(-1.0)

    color_mu_active = "#1e293b"
    color_mu_inactive = "#94a3b8"
    unit_circ = [Point2f(cos(t), sin(t)) for t in range(0, 2π, length = 48)]

    ax_anatomy = Axis(fig[1, 1], title = "Motor Unit Pool, Muscle & Bipolar Recording Setup", titlesize = title_font, aspect = DataAspect())
    hidedecorations!(ax_anatomy)
    hidespines!(ax_anatomy)
    xlims!(ax_anatomy, -5.90, 4.10)
    ylims!(ax_anatomy, -3.35, 2.95)

    y_center = 0.50f0
    h_muscle = 0.86f0

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

    xs_muscle = collect(Float32, range(-2.2f0, 2.2f0, length = 32))
    muscle_arch = [1.0f0 - (x / 2.2f0)^2 for x in xs_muscle]
    function make_muscle_poly(force)
        scale_h = h_muscle * Float32(1.0 + 0.30 * (force / 100.0))
        pts = Vector{Point2f}(undef, 64)
        @inbounds for i = 1:32
            dy = scale_h * muscle_arch[i]
            pts[i] = Point2f(xs_muscle[i], y_center + dy)
            pts[65-i] = Point2f(xs_muscle[i], y_center - dy)
        end
        return pts
    end

    poly!(ax_anatomy, lift(make_muscle_poly, current_disp_force), color = "#991b1b", strokecolor = "#450a0a", strokewidth = 2.0)

    xs_stri = collect(Float32, range(-2.1f0, 2.1f0, length = 40))
    stri_arch = [1.0f0 - (x / 2.2f0)^2 for x in xs_stri]
    for s_factor in (-0.75f0, -0.55f0, -0.35f0, -0.18f0, 0.0f0, 0.18f0, 0.35f0, 0.55f0, 0.75f0)
        let sf = s_factor
            s_pts = lift(current_disp_force) do force
                amp = sf * 0.80f0 * Float32(1.0 + 0.30 * (force / 100.0))
                [Point2f(xs_stri[i], y_center + amp * stri_arch[i]) for i in eachindex(xs_stri)]
            end
            lines!(ax_anatomy, s_pts, color = RGBAf(0.3, 0.04, 0.04, 0.45), linewidth = 1.5, linestyle = :dot)
        end
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

    bg_soma_col = lift(f -> f > 0.0 ? "#475569" : "#cbd5e1", current_disp_force)
    for pos in (
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
    )
        poly!(ax_anatomy, [pos + 0.10f0 * p for p in unit_circ], color = bg_soma_col, strokecolor = "#334155", strokewidth = 1.0)
    end

    soma1_pos = Point2f(-4.88, 2.05)
    soma2_pos = Point2f(-4.88, 1.55)
    soma3_pos = Point2f(-4.88, 1.05)
    soma4_pos = Point2f(-4.88, 0.55)
    somaN_pos = Point2f(-4.88, -0.25)
    rep_somas = (soma1_pos, soma2_pos, soma3_pos, soma4_pos, somaN_pos)
    r_soma = 0.17f0

    halo_col = lift(f -> f > 0.0 ? RGBAf(0.12, 0.16, 0.23, 0.15 + 0.30 * (f / 100.0)) : RGBAf(0, 0, 0, 0), current_disp_force)
    soma_col = lift(f -> f > 0.0 ? color_mu_active : color_mu_inactive, current_disp_force)
    for pos in rep_somas
        poly!(ax_anatomy, [pos + (r_soma + 0.12f0) * p for p in unit_circ], color = halo_col)
        poly!(ax_anatomy, [pos + r_soma * p for p in unit_circ], color = soma_col, strokecolor = "#334155", strokewidth = 1.5)
    end

    for (pos, lbl) in ((soma1_pos, "MU 1"), (soma2_pos, "MU 2"), (soma3_pos, "MU 3"), (soma4_pos, "MU 4"), (somaN_pos, "MU n"))
        text!(ax_anatomy, -4.50, pos[2], text = lbl, align = (:left, :center), font = :bold, fontsize = 8.0, color = "#0f172a")
    end
    text!(ax_anatomy, -4.50, 0.15, text = "⋮", align = (:left, :center), font = :bold, fontsize = 12.0, color = "#64748b")

    poly!(
        ax_anatomy,
        [
            Point2f(-4.25, 2.00),
            Point2f(-3.50, 1.55),
            Point2f(-2.65, 1.20),
            Point2f(-2.65, 0.15),
            Point2f(-3.50, 0.35),
            Point2f(-4.25, -0.45),
        ],
        color = RGBAf(0.92, 0.95, 0.98, 0.65),
        strokecolor = "#94a3b8",
        strokewidth = 1.2,
        linestyle = :dash,
    )
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

    axon_col = lift(f -> f > 0.0 ? "#334155" : "#cbd5e1", current_disp_force)
    for axon in (
        [soma1_pos, Point2f(-3.7, 1.45), Point2f(-2.6, 1.10), Point2f(-1.3, 1.00)],
        [soma2_pos, Point2f(-3.7, 1.20), Point2f(-2.6, 0.90), Point2f(-0.9, 0.78)],
        [soma3_pos, Point2f(-3.7, 0.95), Point2f(-2.6, 0.70), Point2f(-1.5, 0.58)],
        [soma4_pos, Point2f(-3.7, 0.70), Point2f(-2.6, 0.50), Point2f(-1.1, 0.38)],
        [Point2f(-5.30, 0.20), Point2f(-3.7, 0.45), Point2f(-2.6, 0.30), Point2f(-1.4, 0.18)],
        [somaN_pos, Point2f(-3.7, 0.20), Point2f(-2.6, 0.10), Point2f(-1.0, -0.05)],
    )
        lines!(ax_anatomy, axon, color = axon_col, linewidth = 2.0, linestyle = :dot)
    end

    scatter!(
        ax_anatomy,
        [Point2f(-1.3, 1.00), Point2f(-0.9, 0.78), Point2f(-1.5, 0.58), Point2f(-1.1, 0.38), Point2f(-1.4, 0.18), Point2f(-1.0, -0.05)],
        marker = :utriangle,
        markersize = 11,
        color = soma_col,
    )
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

    y_skin = 2.15f0
    lines!(ax_anatomy, [-2.2, 2.2], [y_skin, y_skin], color = "#fcd34d", linewidth = 3.0)
    text!(ax_anatomy, 2.30, y_skin, text = "Skin Surface", align = (:left, :center), font = :italic, fontsize = 9.0, color = "#78350f")

    poly!(ax_anatomy, Rect2f(-0.95, y_skin, 0.40, 0.28), color = "#2563eb", strokewidth = 1.2, strokecolor = "#1e293b")
    text!(ax_anatomy, -0.75, y_skin + 0.38, text = "[-]", align = (:center, :center), font = :bold, fontsize = 11.5, color = "#2563eb")
    poly!(ax_anatomy, Rect2f(0.55, y_skin, 0.40, 0.28), color = "#dc2626", strokewidth = 1.2, strokecolor = "#1e293b")
    text!(ax_anatomy, 0.75, y_skin + 0.38, text = "[+]", align = (:center, :center), font = :bold, fontsize = 11.5, color = "#dc2626")

    y_amp = -2.10f0
    poly!(
        ax_anatomy,
        [Point2f(-0.70, y_amp + 0.45), Point2f(-0.70, y_amp - 0.45), Point2f(0.70, y_amp)],
        color = "#f8fafc",
        strokewidth = 2.0,
        strokecolor = "#1e293b",
    )
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

    lines!(
        ax_anatomy,
        [
            Point2f(0.75, y_skin),
            Point2f(1.60, y_skin),
            Point2f(1.60, -1.35),
            Point2f(-1.05, -1.35),
            Point2f(-1.05, y_amp + 0.22),
            Point2f(-0.70, y_amp + 0.22),
        ],
        color = "#dc2626",
        linewidth = 2.0,
    )
    lines!(
        ax_anatomy,
        [Point2f(-0.75, y_skin), Point2f(-1.60, y_skin), Point2f(-1.60, y_amp - 0.22), Point2f(-0.70, y_amp - 0.22)],
        color = "#2563eb",
        linewidth = 2.0,
    )
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

    ctrl_layout = fig[2, 1] = GridLayout(tellwidth = false, halign = :center)
    slider_sub = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    function _emg_slider(row::Int, lbl_str::String, rng, startv::Float64, unit_str::String)
        Label(slider_sub[row, 1], text = lbl_str, font = :bold, fontsize = label_font, halign = :right)
        sl = Slider(slider_sub[row, 2], range = rng, startvalue = startv, width = 210)
        Label(
            slider_sub[row, 3],
            text = lift(v -> "$(round(Int, v)) $unit_str", sl.value),
            width = 76,
            halign = :left,
            fontsize = label_font,
        )
        return sl
    end

    sl_force = _emg_slider(1, "Force:", 0.0:1.0:100.0, 45.0, "% MVC")
    sl_dur = _emg_slider(2, "Duration:", 200.0:50.0:1400.0, 700.0, "ms")
    sl_rms = _emg_slider(3, "RMS Win:", 20.0:5.0:250.0, 75.0, "ms")
    sl_noise = _emg_slider(4, "Noise:", 0.0:1.0:100.0, 4.0, "µV")

    colgap!(slider_sub, 10)
    rowgap!(slider_sub, 4)

    btn_sub = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    btn_play = Button(btn_sub[1, 1], label = "▶ Play Contraction", buttoncolor = "#dcfce7", fontsize = btn_font, width = 200, height = 34)
    rowgap!(ctrl_layout, 8)

    fs = 1000.0
    dt = 1.0 / fs
    t_vec = collect(range(0.0, 2.0, step = dt))
    N_pts = length(t_vec)
    zero_vec = zeros(Float64, N_pts)

    base_noise =
        [0.65 * sin(2π * 23.4 * tt) * cos(2π * 61.2 * tt) + 0.35 * sin(2π * 117.8 * tt) + 0.25 * sin(2π * 243.1 * tt) for tt in t_vec]

    coprimes = (
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
    )
    N_units = length(coprimes)
    jitter_phases = [[0.20 * sin(2π * k / coprimes[idx] + 1.3 * idx) for k = 1:N_pts] for idx = 1:N_units]

    fiber_pools = [
        begin
            sigma = 0.0018 + 0.0018 * (i / N_units)
            w_samples = round(Int, 5.0 * sigma * fs)
            kernel = [
                begin
                    u = (k * dt) / sigma
                    (-u + 0.6 * u^3) * exp(-0.5 * u^2)
                end for k = (-w_samples):w_samples
            ]
            (
                thresh = 2.0 + 78.0 * ((i - 1) / (N_units - 1))^1.2,
                f_base = 8.0 + 4.0 * (i / N_units),
                f_max = 24.0 + 32.0 * (i / N_units),
                amp = 30.0 + 220.0 * ((i - 1) / (N_units - 1))^1.4,
                w_samples = w_samples,
                kernel = kernel,
            )
        end for i = 1:N_units
    ]

    csum_buf = Vector{Float64}(undef, N_pts + 1)

    function compute_emg_signals(force_val, dur_val, rms_win_val, noise_val)
        t_on = 0.40
        t_off = t_on + dur_val / 1000.0

        force_t = Vector{Float64}(undef, N_pts)
        @inbounds for i = 1:N_pts
            tt = t_vec[i]
            s_on = 1.0 / (1.0 + exp(-(tt - t_on) / 0.025))
            s_off = 1.0 / (1.0 + exp((tt - t_off) / 0.025))
            force_t[i] = force_val * s_on * s_off
        end

        raw_sig = zeros(Float64, N_pts)
        @inbounds for u_idx = 1:N_units
            p = fiber_pools[u_idx]
            j_phase = jitter_phases[u_idx]
            kern = p.kernel
            w_s = p.w_samples
            phi = 0.0
            inv_span = 1.0 / (100.0 - p.thresh)

            for i = 2:N_pts
                F = force_t[i]
                if F >= p.thresh
                    drive = (F - p.thresh) * inv_span
                    freq = muladd(p.f_max - p.f_base, drive, p.f_base)
                    dphi = freq * dt
                    phi += dphi
                    phi_eff = phi + j_phase[i]
                    if floor(phi_eff) > floor(phi_eff - dphi)
                        amp_unit = p.amp * (0.60 + 0.40 * sqrt(drive))
                        k_lo = max(-w_s, 1 - i)
                        k_hi = min(w_s, N_pts - i)
                        for k = k_lo:k_hi
                            raw_sig[i+k] = muladd(amp_unit, kern[k+w_s+1], raw_sig[i+k])
                        end
                    end
                end
            end
        end

        rect_sig = Vector{Float64}(undef, N_pts)
        integ_sig = Vector{Float64}(undef, N_pts)
        csum_buf[1] = 0.0
        acc_sq = 0.0
        acc_int = 0.0
        @inbounds for i = 1:N_pts
            v = muladd(noise_val, base_noise[i], raw_sig[i])
            raw_sig[i] = v
            rv = abs(v)
            rect_sig[i] = rv
            acc_sq = muladd(v, v, acc_sq)
            csum_buf[i+1] = acc_sq
            acc_int = muladd(rv, dt, acc_int)
            integ_sig[i] = acc_int
        end

        half_w = div(max(3, round(Int, (rms_win_val / 1000.0) * fs)), 2)
        rms_sig = Vector{Float64}(undef, N_pts)
        @inbounds for i = 1:N_pts
            i1 = max(1, i - half_w)
            i2 = min(N_pts, i + half_w)
            rms_sig[i] = sqrt(max(0.0, csum_buf[i2+1] - csum_buf[i1]) / (i2 - i1 + 1))
        end

        return (; raw = raw_sig, rect = rect_sig, rms = rms_sig, integ = integ_sig, max_integ = integ_sig[end], force_curve = force_t)
    end

    signals = lift(compute_emg_signals, target_force, duration_ms, rms_window_ms, noise_amp)

    pipe_layout = fig[1:2, 2] = GridLayout()
    pipe_ax_kw = (titlesize = title_font, ylabelsize = label_font, yticklabelsize = tick_font, ylabelpadding = 16.0, yticklabelspace = 48.0)

    ax_raw = Axis(pipe_layout[1, 1]; title = "1. Raw Surface EMG [µV]", ylabel = "Potential [µV]", pipe_ax_kw...)
    xlims!(ax_raw, 0.0, 2.0)
    ylims!(ax_raw, -700, 700)
    hidexdecorations!(ax_raw, grid = true, ticks = true)
    hlines!(ax_raw, [0.0], color = :gray75, linestyle = :dash, linewidth = 1.0)

    ax_rect = Axis(pipe_layout[2, 1]; title = "2. Rectify [µV] — |V(t)|", ylabel = "|V| [µV]", pipe_ax_kw...)
    xlims!(ax_rect, 0.0, 2.0)
    ylims!(ax_rect, 0, 700)
    hidexdecorations!(ax_rect, grid = true, ticks = true)

    ax_filt = Axis(pipe_layout[3, 1]; title = "3. Filter [µV] — Linear Envelope", ylabel = "Envelope [µV]", pipe_ax_kw...)
    xlims!(ax_filt, 0.0, 2.0)
    ylims!(ax_filt, 0, 320)
    hidexdecorations!(ax_filt, grid = true, ticks = true)

    ax_integ = Axis(
        pipe_layout[4, 1];
        title = "4. Integrate [µV·s] — iEMG = ∫ |V| dt",
        xlabel = "Time [s]",
        ylabel = "iEMG [µV·s]",
        xlabelsize = label_font,
        xticklabelsize = tick_font,
        pipe_ax_kw...,
    )
    xlims!(ax_integ, 0.0, 2.0)
    ylims!(ax_integ, 0, 240)

    lines!(ax_raw, t_vec, lift(s -> s.raw, signals), color = "#0f172a", linewidth = 1.3)

    band!(ax_rect, t_vec, zero_vec, lift(s -> s.rect, signals), color = RGBAf(0.0, 0.0, 0.0, 0.08))
    lines!(ax_rect, t_vec, lift(s -> s.rect, signals), color = "#0f172a", linewidth = 1.3)

    band!(ax_filt, t_vec, zero_vec, lift(s -> s.rms, signals), color = RGBAf(0.0, 0.0, 0.0, 0.10))
    lines!(ax_filt, t_vec, lift(s -> s.rms, signals), color = "#0f172a", linewidth = 2.2, label = "Filtered Envelope")

    force_scaled = lift(s -> (s.force_curve ./ 100.0) .* (0.92 * max(10.0, maximum(s.rms))), signals)
    lines!(ax_filt, t_vec, force_scaled, color = "#64748b", linestyle = :dash, linewidth = 2.0, label = "Muscle Force")
    axislegend(ax_filt, position = :lt, labelsize = 10.5, backgroundcolor = RGBAf(1, 1, 1, 0.85))

    band!(ax_integ, t_vec, zero_vec, lift(s -> s.integ, signals), color = RGBAf(0.0, 0.0, 0.0, 0.08))
    lines!(ax_integ, t_vec, lift(s -> s.integ, signals), color = "#0f172a", linewidth = 2.4, label = "iEMG")
    text!(
        ax_integ,
        0.35,
        lift(s -> clamp(0.80 * s.max_integ, 35.0, 200.0), signals),
        text = lift(s -> "iEMG = $(round(s.max_integ, digits=2)) µV·s", signals),
        font = :bold,
        fontsize = 12.5,
        color = "#0f172a",
    )

    for (ax, y_lo, y_hi) in ((ax_raw, -700, 700), (ax_rect, 0, 700), (ax_filt, 0, 320), (ax_integ, 0, 240))
        let lo = Float32(y_lo), hi = Float32(y_hi)
            lines!(
                ax,
                lift(pt -> pt < 0 ? [Point2f(NaN, NaN)] : [Point2f(pt, lo), Point2f(pt, hi)], playback_time),
                color = "#475569",
                linewidth = 2.0,
                linestyle = :dash,
            )
        end
    end

    linkxaxes!(ax_raw, ax_rect, ax_filt, ax_integ)

    is_animating = Ref(false)

    on(sl_force.value) do val
        target_force[] = val
        is_animating[] || (current_disp_force[] = val)
    end
    on(val -> (duration_ms[] = val), sl_dur.value)
    on(val -> (rms_window_ms[] = val), sl_rms.value)
    on(val -> (noise_amp[] = val), sl_noise.value)

    on(btn_play.clicks) do _
        is_animating[] && return
        is_animating[] = true
        @async begin
            try
                target = target_force[]
                t_on = 0.40
                t_off = t_on + duration_ms[] / 1000.0
                n_frames = 65
                dt_frame = 2.0 / n_frames
                for frame = 0:n_frames
                    (!is_animating[] || !events(fig).window_open[]) && break
                    t_curr = frame * dt_frame
                    playback_time[] = t_curr
                    s_on = 1.0 / (1.0 + exp(-(t_curr - t_on) / 0.025))
                    s_off = 1.0 / (1.0 + exp((t_curr - t_off) / 0.025))
                    current_disp_force[] = target * s_on * s_off
                    sleep(1.0 / 30.0)
                end
                if events(fig).window_open[]
                    playback_time[] = -1.0
                    current_disp_force[] = target
                end
            finally
                is_animating[] = false
            end
        end
    end

    rowsize!(fig.layout, 1, Relative(0.75))
    rowsize!(fig.layout, 2, Relative(0.25))
    colsize!(fig.layout, 1, Relative(0.53))
    colsize!(fig.layout, 2, Relative(0.47))

    rowsize!(pipe_layout, 1, Relative(0.26))
    rowsize!(pipe_layout, 2, Relative(0.24))
    rowsize!(pipe_layout, 3, Relative(0.25))
    rowsize!(pipe_layout, 4, Relative(0.25))

    _display_demo_figure(fig)
    return fig, ax_anatomy, ax_raw, ax_rect, ax_filt, ax_integ
end
