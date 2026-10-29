"""
    signal_example_ecg_hrv()

Interactive Heart Rate Variability (HRV) Simulator.

Creates a clean, intuitive psychophysiology demonstration illustrating how breathing
dynamically modulates heart rate.

## Educational Concepts

1. **The Heart is Not a Metronome**: In healthy individuals, the heart does not beat at a
   constant rhythm. The interval between beats changes continuously. This is **Heart Rate Variability (HRV)**.
2. **Breathing Changes Heart Rate**:
   - **Inhale**: Heart rate speeds up (beats get closer together).
   - **Exhale**: Heart rate slows down (beats spread farther apart).
3. **Relaxation vs. Stress**:
   - **Relaxed**: Heart rate oscillates with large, flexible waves (**High HRV**).
   - **Stressed**: Heart rate flattens into a rigid, invariant rhythm (**Low HRV**).

## Interactive Backend (GLMakie)

```julia
using GLMakie
using EegFun

signal_example_ecg_hrv()
```

# Returns
- `fig::Figure`: Makie figure object.
- `ax_ecg::Axis`: ECG strip with beat-to-beat intervals.
- `ax_hr::Axis`: Heart rate and breathing synchrony plot.
"""
function signal_example_ecg_hrv()
    _check_interactive_demo_backend("Interactive ECG HRV demo") || return nothing

    _set_window_title("Interactive HRV Simulator")

    fig =
        _eegfun_figure(size = (1280, 860), title = "Interactive HRV Simulator", backgroundcolor = :white, figure_padding = (24, 24, 24, 24))

    # ── Adaptive Font Sizing ─────────────────────────────────────────────────
    function setup_adaptive_sizing(fig)
        title_font = Observable(15)
        label_font = Observable(12.5)
        tick_font  = Observable(10.5)
        btn_font   = Observable(12.0)

        on(fig.scene.viewport) do area
            scale        = area.widths[1] / 1280
            title_font[] = max(12, round(Int, 15 * scale))
            label_font[] = max(10, round(Int, 12.5 * scale))
            tick_font[]  = max(9, round(Int, 10.5 * scale))
            btn_font[]   = max(10, round(Int, 12.0 * scale))
        end

        return title_font, label_font, tick_font, btn_font
    end

    title_font, label_font, tick_font, btn_font = setup_adaptive_sizing(fig)

    # ── Observables ──────────────────────────────────────────────────────────
    obs_mean_hr   = Observable(72.0)   # Baseline heart rate (BPM)
    obs_resp_rate = Observable(12.0)   # Breaths per minute
    obs_vagal_pct = Observable(75.0)   # Relaxation / Vagal tone (0 to 100%)

    # ── Simulation Engine ────────────────────────────────────────────────────
    function compute_hrv_data(mean_hr, resp_bpm, vagal_pct)
        f_resp = resp_bpm / 60.0  # Breathing frequency (Hz)

        # Heart rate swing amplitude in BPM:
        # At 100% relaxation (high HRV): +/- 10 BPM swing with each breath
        # At 0% relaxation (high stress): +/- 0.8 BPM swing (rigid metronome)
        amp_bpm = 0.8 + (vagal_pct / 100.0) * 9.5

        # 1. Continuous Time Series for Synchrony Plot (20 seconds)
        t_span = 20.0
        n_hr_pts = 400
        t_hr = collect(range(0.0, t_span, length = n_hr_pts))

        # Heart Rate curve: peaks during inhalation (with small natural 0.35s delay)
        hr_curve = mean_hr .+ amp_bpm .* sin.(2π .* f_resp .* (t_hr .- 0.35))

        # Breathing curve: normalized to visually match the heart rate swing
        # so students can immediately see the two curves move together
        resp_curve = mean_hr .+ amp_bpm .* sin.(2π .* f_resp .* t_hr)

        # 2. Discrete Beat Generation for ECG Strip (6 seconds)
        t_ecg_win = 6.0
        t_beats = Float64[]
        curr_t = 0.10
        while curr_t < t_ecg_win + 0.8
            inst_hr = mean_hr + amp_bpm * sin(2π * f_resp * (curr_t - 0.35))
            ibi_s = 60.0 / max(35.0, inst_hr)
            push!(t_beats, curr_t)
            curr_t += ibi_s
        end

        ibis_ms = length(t_beats) > 1 ? diff(t_beats) .* 1000.0 : Float64[]

        # Synthesize Lead II ECG waveform for the 6-second strip
        fs_ecg = 250.0
        n_ecg = round(Int, t_ecg_win * fs_ecg)
        t_ecg = collect(range(0.0, t_ecg_win, length = n_ecg))
        ecg_sig = zeros(Float64, n_ecg)
        r_peaks = Point2f[]
        bracket_lines = Point2f[]
        ibi_labels = Tuple{Point2f,String}[]

        for i = 1:length(t_beats)
            tb = t_beats[i]
            tR = tb + 0.25
            if 0.0 <= tR <= t_ecg_win
                push!(r_peaks, Point2f(tR, 1.35))

                if i < length(t_beats)
                    tR_next = t_beats[i+1] + 0.25
                    if tR_next <= t_ecg_win && i <= length(ibis_ms)
                        # Measurement bracket between consecutive R-peaks
                        y_bar = 1.46
                        # Horizontal bar
                        push!(bracket_lines, Point2f(tR, y_bar))
                        push!(bracket_lines, Point2f(tR_next, y_bar))
                        # Downward end-ticks
                        push!(bracket_lines, Point2f(tR, y_bar))
                        push!(bracket_lines, Point2f(tR, y_bar - 0.05))
                        push!(bracket_lines, Point2f(tR_next, y_bar))
                        push!(bracket_lines, Point2f(tR_next, y_bar - 0.05))

                        # Text label centered above the bracket
                        t_mid = (tR + tR_next) / 2.0
                        push!(ibi_labels, (Point2f(t_mid, y_bar + 0.06), "$(round(Int, ibis_ms[i])) ms"))
                    end
                end
            end

            # ECG wave synthesis (P, Q, R, S, T)
            idx_start = max(1, round(Int, tb * fs_ecg) + 1)
            idx_end   = min(n_ecg, round(Int, (tb + 0.65) * fs_ecg))
            for idx = idx_start:idx_end
                dt = t_ecg[idx] - tb
                if 0.0 <= dt <= 0.65
                    gP = 0.18 * exp(-((dt - 0.12) / 0.022)^2)
                    gQ = -0.16 * exp(-((dt - 0.22) / 0.010)^2)
                    gR = 1.35 * exp(-((dt - 0.25) / 0.014)^2)
                    gS = -0.32 * exp(-((dt - 0.28) / 0.012)^2)
                    gT = 0.38 * exp(-((dt - 0.45) / 0.048)^2)
                    ecg_sig[idx] += gP + gQ + gR + gS + gT
                end
            end
        end

        # Subtle clean noise
        for idx = 1:n_ecg
            ecg_sig[idx] += 0.012 * sin(idx * 13.7)
        end

        ecg_pts  = Point2f.(t_ecg, ecg_sig)
        hr_pts   = Point2f.(t_hr, hr_curve)
        resp_pts = Point2f.(t_hr, resp_curve)

        status_text =
            vagal_pct >= 50 ? "Relaxed State: Heart rate flexibly speeds up & slows down with breathing (High HRV)" :
            "Stressed State: Heart rate is rigid and metronome-like with minimal variation (Low HRV)"

        status_color = vagal_pct >= 50 ? "#047857" : "#b91c1c"

        return (;
            mean_hr,
            resp_bpm,
            vagal_pct,
            amp_bpm,
            t_ecg_win,
            ecg_pts,
            r_peaks,
            bracket_lines,
            ibi_labels,
            t_hr,
            hr_curve,
            hr_pts,
            resp_pts,
            status_text,
            status_color,
        )
    end

    sim_data = lift(obs_mean_hr, obs_resp_rate, obs_vagal_pct) do hr, resp, vag
        compute_hrv_data(hr, resp, vag)
    end

    # ── Top Status Banner ────────────────────────────────────────────────────
    status_lbl = Label(
        fig[1, 1],
        text = lift(d -> d.status_text, sim_data),
        font = :bold,
        fontsize = lift(f -> f * 1.05, label_font),
        color = lift(d -> d.status_color, sim_data),
        halign = :center,
    )

    # ── Top Plot [2, 1]: Clean ECG Strip with Beat Intervals ─────────────────
    ax_ecg = Axis(
        fig[2, 1],
        title = "ECG Strip: Time Between Consecutive Heartbeats",
        xlabel = "Time [seconds]",
        ylabel = "ECG Signal",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        titlegap = 10,
        xlabelpadding = 8,
        yticklabelsvisible = false,
        yticksvisible = false,
        backgroundcolor = :white,
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
        xgridwidth = 0.8,
        ygridwidth = 0.8,
    )
    xlims!(ax_ecg, 0.0, 6.0)
    ylims!(ax_ecg, -0.55, 1.90)

    # Baseline zero line
    hlines!(ax_ecg, [0.0], color = "#cbd5e1", linestyle = :dash, linewidth = 1.0)

    # ECG trace
    lines!(ax_ecg, lift(d -> d.ecg_pts, sim_data), color = "#e11d48", linewidth = 2.2)

    # R-peak markers
    scatter!(ax_ecg, lift(d -> d.r_peaks, sim_data), color = "#e11d48", strokecolor = "#991b1b", strokewidth = 1.2, markersize = 10)

    # Interval brackets between peaks
    linesegments!(ax_ecg, lift(d -> d.bracket_lines, sim_data), color = "#0284c7", linewidth = 1.5)

    # Interval text labels (up to 12 intervals)
    for k = 1:12
        txt_obs = lift(sim_data) do d
            k <= length(d.ibi_labels) ? d.ibi_labels[k][2] : ""
        end
        pos_obs = lift(sim_data) do d
            k <= length(d.ibi_labels) ? d.ibi_labels[k][1] : Point2f(NaN, NaN)
        end
        text!(
            ax_ecg,
            pos_obs,
            text = txt_obs,
            align = (:center, :bottom),
            font = :bold,
            fontsize = lift(f -> f * 0.95, label_font),
            color = "#0284c7",
        )
    end

    # ── Bottom Plot [3, 1]: Heart Rate & Breathing Synchrony ─────────────────
    ax_hr = Axis(
        fig[3, 1],
        title = "Breathing & Heart Rate Synchrony",
        xlabel = "Time [seconds]",
        ylabel = "Heart Rate [BPM]",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        titlegap = 10,
        xlabelpadding = 8,
        ylabelpadding = 10,
        backgroundcolor = "#f8fafc",
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
        xgridwidth = 0.8,
        ygridwidth = 0.8,
    )
    xlims!(ax_hr, 0.0, 20.0)

    # Respiration curve (smooth dashed teal line)
    line_resp = lines!(ax_hr, lift(d -> d.resp_pts, sim_data), color = "#0891b2", linewidth = 2.4, linestyle = :dash)

    # Instantaneous Heart Rate curve (bold blue line)
    line_hr = lines!(ax_hr, lift(d -> d.hr_pts, sim_data), color = "#2563eb", linewidth = 3.0)

    # Dynamic limits for HR axis: generous headroom (+13 BPM) so curves never collide with legend
    on(sim_data) do d
        hr_min = minimum(d.hr_curve)
        hr_max = maximum(d.hr_curve)
        ylims!(ax_hr, hr_min - 5.0, hr_max + 13.0)
    end

    # Clean, uncrowded Legend with comfortable margins and internal padding
    axislegend(
        ax_hr,
        [line_hr, line_resp],
        ["Heart Rate (BPM)", "Breathing (Inhale / Exhale)"],
        position = :rt,
        labelsize = lift(f -> min(Float64(f) * 0.90, 11.5), label_font),
        padding = (12.0, 12.0, 8.0, 8.0),
        margin = (14.0, 14.0, 10.0, 10.0),
        backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.92),
        framecolor = "#cbd5e1",
    )

    # ── Controls Panel [4, 1]: 3 Intuitive Sliders & 2 Scenario Buttons ─────
    ctrl_layout = fig[4, 1] = GridLayout(tellwidth = false, halign = :center)

    # Row 1: Sliders
    sl_grid = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(sl_grid[1, 1], text = "Heart Rate:", font = :bold, fontsize = label_font)
    sl_hr = Slider(sl_grid[1, 2], range = 50.0:1.0:110.0, startvalue = 72.0, width = 120)
    lbl_hr = Label(sl_grid[1, 3], text = lift(v -> "$(round(Int, v)) BPM", sl_hr.value), width = 75, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 5], text = "Breathing Speed:", font = :bold, fontsize = label_font)
    sl_resp = Slider(sl_grid[1, 6], range = 6.0:1.0:20.0, startvalue = 12.0, width = 120)
    lbl_resp =
        Label(sl_grid[1, 7], text = lift(v -> "$(round(Int, v)) br/min", sl_resp.value), width = 85, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 9], text = "Relaxation (HRV):", font = :bold, fontsize = label_font)
    sl_vag = Slider(sl_grid[1, 10], range = 0.0:1.0:100.0, startvalue = 75.0, width = 120)
    lbl_vag = Label(sl_grid[1, 11], text = lift(v -> "$(round(Int, v))%", sl_vag.value), width = 60, fontsize = label_font, halign = :left)

    colgap!(sl_grid, 10)
    colsize!(sl_grid, 4, Fixed(45))
    colsize!(sl_grid, 8, Fixed(45))

    rowgap!(ctrl_layout, 16)

    # Row 2: Scenario Buttons
    btn_grid = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(btn_grid[1, 1], text = "Presets:", font = :bold, fontsize = label_font)
    btn_calm   = Button(btn_grid[1, 2], label = "Relaxed (High HRV)", buttoncolor = "#dcfce7", fontsize = btn_font, width = 230, height = 36)
    btn_stress = Button(btn_grid[1, 3], label = "Stressed (Low HRV)", buttoncolor = "#fee2e2", fontsize = btn_font, width = 230, height = 36)

    colgap!(btn_grid, 18)

    # ── Callbacks ────────────────────────────────────────────────────────────
    on(sl_hr.value) do val
        obs_mean_hr[] = val
    end

    on(sl_resp.value) do val
        obs_resp_rate[] = val
    end

    on(sl_vag.value) do val
        obs_vagal_pct[] = val
    end

    # Scenario Button Callbacks
    on(btn_calm.clicks) do _
        set_close_to!(sl_hr, 68.0)
        set_close_to!(sl_resp, 10.0)
        set_close_to!(sl_vag, 85.0)
    end

    on(btn_stress.clicks) do _
        set_close_to!(sl_hr, 95.0)
        set_close_to!(sl_resp, 18.0)
        set_close_to!(sl_vag, 8.0)
    end

    # Layout Sizing & Spacing
    rowsize!(fig.layout, 1, Fixed(28))       # Status banner
    rowsize!(fig.layout, 2, Relative(0.38))  # ECG strip
    rowsize!(fig.layout, 3, Relative(0.38))  # HR & breathing synchrony
    rowsize!(fig.layout, 4, Fixed(130))      # Controls panel (guaranteed space, no cutoff)

    rowgap!(fig.layout, 12)
    rowgap!(fig.layout, 3, Fixed(26))        # Generous separation between HR axis and controls

    _display_demo_figure(fig)

    return fig, ax_ecg, ax_hr
end
