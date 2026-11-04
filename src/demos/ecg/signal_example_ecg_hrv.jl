"""
    signal_example_ecg_hrv()

Interactive Heart Rate Variability (HRV) simulator illustrating respiratory sinus arrhythmia
(breathing-coupled R–R interval modulation) across relaxed (high vagal tone) and stressed states.

```julia
using GLMakie, EegFun
EegFun.signal_example_ecg_hrv()
```

# Returns
- `fig::Figure`: Makie figure object
- `ax_ecg::Axis`: ECG strip with beat-to-beat intervals
- `ax_hr::Axis`: Heart rate and breathing synchrony plot
"""
function signal_example_ecg_hrv()
    _check_interactive_demo_backend("Interactive ECG HRV demo") || return nothing
    _set_window_title("Interactive HRV Simulator")

    fig =
        _eegfun_figure(size = (1280, 860), title = "Interactive HRV Simulator", backgroundcolor = :white, figure_padding = (24, 24, 24, 24))
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1280, (15.0, 12.0), (12.5, 10.0), (10.5, 9.0), (12.0, 10.0))

    obs_mean_hr = Observable(72.0)
    obs_resp_rate = Observable(12.0)
    obs_vagal_pct = Observable(75.0)

    t_hr = collect(range(0.0, 20.0, length = 400))
    t_ecg_win = 6.0
    fs_ecg = 250.0
    n_ecg = round(Int, t_ecg_win * fs_ecg)
    t_ecg = collect(range(0.0, t_ecg_win, length = n_ecg))
    base_ecg_noise = [0.012 * sin(idx * 13.7) for idx = 1:n_ecg]

    function compute_hrv_data(mean_hr, resp_bpm, vagal_pct)
        f_resp = resp_bpm / 60.0
        amp_bpm = 0.8 + (vagal_pct / 100.0) * 9.5

        hr_curve = @. mean_hr + amp_bpm * sin(2π * f_resp * (t_hr - 0.35))
        resp_curve = @. mean_hr + amp_bpm * sin(2π * f_resp * t_hr)

        t_beats = Float64[]
        curr_t = 0.10
        while curr_t < t_ecg_win + 0.8
            inst_hr = mean_hr + amp_bpm * sin(2π * f_resp * (curr_t - 0.35))
            push!(t_beats, curr_t)
            curr_t += 60.0 / max(35.0, inst_hr)
        end

        ecg_sig = copy(base_ecg_noise)
        r_peaks = Point2f[]
        bracket_lines = Point2f[]
        ibi_pos = Point2f[]
        ibi_strs = String[]

        @inbounds for i in eachindex(t_beats)
            tb = t_beats[i]
            tR = tb + 0.25
            if 0.0 <= tR <= t_ecg_win
                push!(r_peaks, Point2f(tR, 1.35))
                if i < length(t_beats)
                    tR_next = t_beats[i+1] + 0.25
                    if tR_next <= t_ecg_win
                        y_bar = 1.46f0
                        append!(
                            bracket_lines,
                            (
                                Point2f(tR, y_bar),
                                Point2f(tR_next, y_bar),
                                Point2f(tR, y_bar),
                                Point2f(tR, y_bar - 0.05f0),
                                Point2f(tR_next, y_bar),
                                Point2f(tR_next, y_bar - 0.05f0),
                            ),
                        )
                        ibi_ms = round(Int, (t_beats[i+1] - tb) * 1000.0)
                        push!(ibi_pos, Point2f(0.5 * (tR + tR_next), y_bar + 0.06f0))
                        push!(ibi_strs, "$ibi_ms ms")
                    end
                end
            end

            idx_start = max(1, round(Int, tb * fs_ecg) + 1)
            idx_end = min(n_ecg, round(Int, (tb + 0.65) * fs_ecg))
            for idx = idx_start:idx_end
                dt = t_ecg[idx] - tb
                if 0.0 <= dt <= 0.65
                    ecg_sig[idx] +=
                        0.18 * exp(-((dt - 0.12) / 0.022)^2) - 0.16 * exp(-((dt - 0.22) / 0.010)^2) + 1.35 * exp(-((dt - 0.25) / 0.014)^2) -
                        0.32 * exp(-((dt - 0.28) / 0.012)^2) + 0.38 * exp(-((dt - 0.45) / 0.048)^2)
                end
            end
        end

        status_text =
            vagal_pct >= 50 ? "Relaxed State: Heart rate flexibly speeds up & slows down with breathing (High HRV)" :
            "Stressed State: Heart rate is rigid and metronome-like with minimal variation (Low HRV)"
        status_color = vagal_pct >= 50 ? "#047857" : "#b91c1c"

        return (;
            ecg_pts = Point2f.(t_ecg, ecg_sig),
            r_peaks,
            bracket_lines,
            ibi_pos,
            ibi_strs,
            hr_curve,
            hr_pts = Point2f.(t_hr, hr_curve),
            resp_pts = Point2f.(t_hr, resp_curve),
            status_text,
            status_color,
        )
    end

    sim_data = lift(compute_hrv_data, obs_mean_hr, obs_resp_rate, obs_vagal_pct)

    Label(
        fig[1, 1],
        text = lift(d -> d.status_text, sim_data),
        font = :bold,
        fontsize = lift(f -> f * 1.05, label_font),
        color = lift(d -> d.status_color, sim_data),
        halign = :center,
    )

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

    hlines!(ax_ecg, [0.0], color = "#cbd5e1", linestyle = :dash, linewidth = 1.0)
    lines!(ax_ecg, lift(d -> d.ecg_pts, sim_data), color = "#e11d48", linewidth = 2.2)
    scatter!(ax_ecg, lift(d -> d.r_peaks, sim_data), color = "#e11d48", strokecolor = "#991b1b", strokewidth = 1.2, markersize = 10)
    linesegments!(ax_ecg, lift(d -> d.bracket_lines, sim_data), color = "#0284c7", linewidth = 1.5)

    ibi_font = lift(f -> f * 0.95, label_font)
    text!(
        ax_ecg,
        lift(d -> d.ibi_pos, sim_data),
        text = lift(d -> d.ibi_strs, sim_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = ibi_font,
        color = "#0284c7",
    )

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

    line_resp = lines!(ax_hr, lift(d -> d.resp_pts, sim_data), color = "#0891b2", linewidth = 2.4, linestyle = :dash)
    line_hr = lines!(ax_hr, lift(d -> d.hr_pts, sim_data), color = "#2563eb", linewidth = 3.0)

    on(sim_data) do d
        hr_min, hr_max = extrema(d.hr_curve)
        ylims!(ax_hr, hr_min - 5.0, hr_max + 13.0)
    end
    let (hr_min0, hr_max0) = extrema(sim_data[].hr_curve)
        ylims!(ax_hr, hr_min0 - 5.0, hr_max0 + 13.0)
    end

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

    ctrl_layout = fig[4, 1] = GridLayout(tellwidth = false, halign = :center)
    sl_grid = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(sl_grid[1, 1], text = "Heart Rate:", font = :bold, fontsize = label_font)
    sl_hr = Slider(sl_grid[1, 2], range = 50.0:1.0:110.0, startvalue = 72.0, width = 120)
    Label(sl_grid[1, 3], text = lift(v -> "$(round(Int, v)) BPM", sl_hr.value), width = 75, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 5], text = "Breathing Speed:", font = :bold, fontsize = label_font)
    sl_resp = Slider(sl_grid[1, 6], range = 6.0:1.0:20.0, startvalue = 12.0, width = 120)
    Label(sl_grid[1, 7], text = lift(v -> "$(round(Int, v)) br/min", sl_resp.value), width = 85, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 9], text = "Relaxation (HRV):", font = :bold, fontsize = label_font)
    sl_vag = Slider(sl_grid[1, 10], range = 0.0:1.0:100.0, startvalue = 75.0, width = 120)
    Label(sl_grid[1, 11], text = lift(v -> "$(round(Int, v))%", sl_vag.value), width = 60, fontsize = label_font, halign = :left)

    colgap!(sl_grid, 10)
    colsize!(sl_grid, 4, Fixed(45))
    colsize!(sl_grid, 8, Fixed(45))
    rowgap!(ctrl_layout, 16)

    btn_grid = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    Label(btn_grid[1, 1], text = "Presets:", font = :bold, fontsize = label_font)
    btn_calm = Button(btn_grid[1, 2], label = "Relaxed (High HRV)", buttoncolor = "#dcfce7", fontsize = btn_font, width = 230, height = 36)
    btn_stress =
        Button(btn_grid[1, 3], label = "Stressed (Low HRV)", buttoncolor = "#fee2e2", fontsize = btn_font, width = 230, height = 36)
    colgap!(btn_grid, 18)

    on(val -> (obs_mean_hr[] = val), sl_hr.value)
    on(val -> (obs_resp_rate[] = val), sl_resp.value)
    on(val -> (obs_vagal_pct[] = val), sl_vag.value)

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

    rowsize!(fig.layout, 1, Fixed(28))
    rowsize!(fig.layout, 2, Relative(0.38))
    rowsize!(fig.layout, 3, Relative(0.38))
    rowsize!(fig.layout, 4, Fixed(130))
    rowgap!(fig.layout, 12)
    rowgap!(fig.layout, 3, Fixed(26))

    _display_demo_figure(fig)
    return fig, ax_ecg, ax_hr
end
