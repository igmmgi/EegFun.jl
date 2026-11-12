"""
    signal_example_eda_habituation()

Interactive Electrodermal Activity (EDA) multi-trial simulator illustrating:
- **Tonic vs. Phasic Decomposition**: Skin Conductance Level (SCL), slow baseline drift, and event-related Skin Conductance Responses (ER-SCRs)
- **Spontaneous Fluctuations (NS-SCRs / SFs)**: Non-specific responses in the absence of stimuli (1–3 / min)
- **Orienting & Habituation (Barry & Sokolov, 1993)**: Exponential decline of \$\\sqrt{\\text{SCR}}\$ and tonic SCL across repeated stimulus presentations (\$S_1 \\dots S_6\$), plus dishabituation to a novel stimulus
- **Temporal Superposition (Overlap)**: How short inter-stimulus intervals (ISIs) cause successive SCRs to ride on the recovery tail of previous responses, distorting trough-to-peak amplitude measurements.

```julia
using GLMakie, EegFun
EegFun.signal_example_eda_habituation()
```

# Returns
- `fig::Figure`: The Makie figure object containing the interactive GUI
- `ax_eda::Axis`: Continuous multi-trial EDA recording strip (Total SC, Tonic SCL, SFs, and stimulus markers)
- `ax_hab_scr::Axis`: Trial-by-trial phasic SCR habituation plot (true underlying vs. measured trough-to-peak)
- `ax_hab_scl::Axis`: Trial-by-trial tonic SCL trajectory plot
"""
function signal_example_eda_habituation()
    _check_interactive_demo_backend("Interactive EDA Habituation demo") || return nothing
    _set_window_title("Interactive EDA Simulator — Tonic/Phasic Decomposition, Habituation & Overlap")

    fig = _eegfun_figure(
        size = (1420, 900),
        title = "Interactive EDA Simulator — Tonic/Phasic Decomposition, Habituation & Overlap",
        backgroundcolor = :white,
        figure_padding = (20, 20, 18, 18),
    )
    title_font, label_font, tick_font, btn_font = _demo_adaptive_fonts(fig, 1420, (15.5, 12.5), (12.5, 10.0), (10.5, 9.5), (12.0, 10.0))

    obs_scr1 = Observable(1.20)
    obs_hab_rate = Observable(45.0)
    obs_isi = Observable(12.0)
    obs_drift = Observable(-1.0)
    obs_sf_rate = Observable(1.0)
    obs_novel_s6 = Observable(false)

    t_max = 90.0
    t_vec = collect(range(0.0, t_max, length = 900))
    n_trials = 6
    trials_x = collect(1:n_trials)
    a_half = 1.6783469900166605

    function eval_scr_unit(t::Real, t_on::Real, t_r::Real, t_h::Real)
        u = t - t_on
        if u <= 0.0
            return 0.0
        elseif u <= t_r
            return 0.5 * (1.0 - cos(π * u / t_r))
        else
            x = (u - t_r) / t_h
            return (1.0 + a_half * x) * exp(-a_half * x)
        end
    end

    # Fixed candidate onset times and amplitudes for Spontaneous Fluctuations (NS-SCRs / SFs)
    sf_specs = [(6.5, 0.28), (28.0, 0.24), (52.5, 0.30), (74.0, 0.25)]

    function compute_hab_data(scr1, hab_pct, isi, drift_per_min, sf_per_min, novel_s6)
        base_scl0 = 8.0
        lat = 1.6
        rise = 2.0
        rec_half = 4.5
        decay_k = (hab_pct / 100.0) * 1.15

        stim_times = [14.0 + (k - 1) * isi for k = 1:n_trials]
        true_amps = Vector{Float64}(undef, n_trials)
        for k = 1:n_trials
            if novel_s6 && k == n_trials && scr1 > 0.02
                true_amps[k] = 0.92 * scr1
            else
                true_amps[k] = scr1 * exp(-decay_k * (k - 1))
            end
        end

        # Tonic SCL with slow drift + mild autonomic arousal bump at start of stimulation
        tonic_vec = [base_scl0 + (drift_per_min / 60.0) * t + 0.35 * scr1 * exp(-((t - 16.0) / 14.0)^2) for t in t_vec]

        # Active Spontaneous Fluctuations (NS-SCRs)
        n_sf = round(Int, clamp(sf_per_min, 0.0, 4.0))
        active_sfs = sf_specs[1:n_sf]

        phasic_er_vec = zeros(Float64, length(t_vec))
        sf_vec = zeros(Float64, length(t_vec))

        for (k, ts) in enumerate(stim_times)
            amp_k = true_amps[k]
            t_on = ts + lat
            @inbounds for i in eachindex(t_vec)
                phasic_er_vec[i] += amp_k * eval_scr_unit(t_vec[i], t_on, rise, rec_half)
            end
        end

        sf_peaks = Point2f[]
        for (t_sf, amp_sf) in active_sfs
            @inbounds for i in eachindex(t_vec)
                sf_vec[i] += amp_sf * eval_scr_unit(t_vec[i], t_sf, 1.6, 3.2)
            end
        end

        total_sc = tonic_vec .+ phasic_er_vec .+ sf_vec

        for (t_sf, _) in active_sfs
            t_pk_sf = t_sf + 1.6
            idx_sf = clamp(round(Int, (t_pk_sf / t_max) * (length(t_vec) - 1)) + 1, 1, length(t_vec))
            push!(sf_peaks, Point2f(t_pk_sf, total_sc[idx_sf]))
        end

        # Compute per-trial measured trough-to-peak SCR and pre-stimulus SCL
        meas_amps = Vector{Float64}(undef, n_trials)
        trial_scl = Vector{Float64}(undef, n_trials)
        peak_pts = Vector{Point2f}(undef, n_trials)
        trough_pts = Vector{Point2f}(undef, n_trials)
        bracket_segs = Point2f[]

        for k = 1:n_trials
            ts = stim_times[k]
            t_on = ts + lat
            t_pk = t_on + rise
            idx_on = clamp(round(Int, (t_on / t_max) * (length(t_vec) - 1)) + 1, 1, length(t_vec))
            idx_pk = clamp(round(Int, (t_pk / t_max) * (length(t_vec) - 1)) + 1, 1, length(t_vec))

            y_on = total_sc[idx_on]
            y_pk = total_sc[idx_pk]
            # Measured trough-to-peak amplitude (clipped at 0 if steep recovery tail dominates)
            meas_amps[k] = max(0.0, y_pk - y_on)
            trial_scl[k] = y_on
            trough_pts[k] = Point2f(t_on, y_on)
            peak_pts[k] = Point2f(t_pk, y_pk)
            append!(bracket_segs, (Point2f(t_pk, y_on), Point2f(t_pk, y_pk)))
        end

        status_text, status_color = if scr1 < 0.05
            ("Non-Responder (~10% of population): Absence of event-related SCRs despite stimulus presentation", "#64748b")
        elseif isi < 7.5
            (
                "Temporal Superposition (Overlap): Short ISI ($(round(isi, digits=1)) s) causes new SCRs to pile onto the previous recovery tail!",
                "#b91c1c",
            )
        elseif novel_s6
            ("Orienting & Dishabituation: SCRs habituate across S1–S5, then recover strongly to Novel Stimulus S6", "#0369a1")
        else
            ("Standard Habituation (Barry & Sokolov, 1993): Exponential SCR decline across repeated stimuli S1–S6", "#047857")
        end

        return (;
            tonic_vec,
            total_sc,
            tonic_pts = Point2f.(t_vec, tonic_vec),
            total_pts = Point2f.(t_vec, total_sc),
            stim_times,
            true_amps,
            meas_amps,
            sqrt_true = sqrt.(max.(0.0, true_amps)),
            sqrt_meas = sqrt.(max.(0.0, meas_amps)),
            trial_scl,
            peak_pts,
            trough_pts,
            bracket_segs,
            sf_peaks,
            status_text,
            status_color,
        )
    end

    sim_data = lift(compute_hab_data, obs_scr1, obs_hab_rate, obs_isi, obs_drift, obs_sf_rate, obs_novel_s6)

    Label(
        fig[1, 1:2],
        text = lift(d -> d.status_text, sim_data),
        font = :bold,
        fontsize = lift(f -> f * 1.05, label_font),
        color = lift(d -> d.status_color, sim_data),
        halign = :center,
    )

    # =========================================================================
    # TOP STRIP: Continuous Multi-Trial EDA Recording (Slides 29, 30, 32, 36)
    # =========================================================================
    ax_eda = Axis(
        fig[2, 1:2],
        title = "Continuous Skin Conductance Recording — Tonic SCL, Phasic ER-SCRs & Spontaneous Fluctuations (SF)",
        xlabel = "Time [s]",
        ylabel = "Skin Conductance [µS]",
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = :white,
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
    )
    xlims!(ax_eda, 0.0, t_max)

    # Shaded phasic area between Tonic SCL and Total SC
    band!(ax_eda, t_vec, lift(d -> d.tonic_vec, sim_data), lift(d -> d.total_sc, sim_data), color = RGBAf(0.02, 0.52, 0.78, 0.14))

    line_tonic = lines!(ax_eda, lift(d -> d.tonic_pts, sim_data), color = "#0284c7", linestyle = :dash, linewidth = 2.2)
    line_total = lines!(ax_eda, lift(d -> d.total_pts, sim_data), color = "#0f172a", linewidth = 2.8)

    # Trough-to-peak vertical measurement brackets
    linesegments!(ax_eda, lift(d -> d.bracket_segs, sim_data), color = "#e11d48", linewidth = 2.0)
    scatter!(ax_eda, lift(d -> d.trough_pts, sim_data), color = "#0284c7", markersize = 8)
    scatter_er =
        scatter!(ax_eda, lift(d -> d.peak_pts, sim_data), color = "#e11d48", strokecolor = :white, strokewidth = 1.2, markersize = 10)

    # Spontaneous Fluctuation (NS-SCR / SF) markers
    scatter_sf = scatter!(
        ax_eda,
        lift(d -> d.sf_peaks, sim_data),
        color = "#d97706",
        marker = :diamond,
        strokecolor = :white,
        strokewidth = 1.2,
        markersize = 12,
    )
    text!(
        ax_eda,
        lift(d -> [Point2f(p[1], p[2] + 0.12) for p in d.sf_peaks], sim_data),
        text = lift(d -> fill("SF", length(d.sf_peaks)), sim_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11,
        color = "#b45309",
    )

    # Stimulus onset vertical tick lines and S1..S6 labels
    stim_segs = lift(sim_data) do d
        ymin_v = minimum(d.tonic_vec) - 0.65
        segs = Point2f[]
        for ts in d.stim_times
            append!(segs, (Point2f(ts, ymin_v), Point2f(ts, ymin_v + 0.32)))
        end
        return segs
    end
    linesegments!(ax_eda, stim_segs, color = "#1e293b", linewidth = 3.2)

    text!(
        ax_eda,
        lift(d -> [Point2f(ts, minimum(d.tonic_vec) - 0.30) for ts in d.stim_times], sim_data),
        text = lift(d -> [obs_novel_s6[] && k == n_trials ? "S6 (Novel)" : "S$k" for k = 1:n_trials], sim_data),
        align = (:center, :bottom),
        font = :bold,
        fontsize = 11.5,
        color = "#1e293b",
    )

    on(sim_data) do d
        ymin_v = minimum(d.tonic_vec) - 0.75
        ymax_v = maximum(d.total_sc) + 0.65
        ylims!(ax_eda, ymin_v, max(ymin_v + 2.2, ymax_v))
    end
    let d0 = sim_data[]
        ymin_v = minimum(d0.tonic_vec) - 0.75
        ymax_v = maximum(d0.total_sc) + 0.65
        ylims!(ax_eda, ymin_v, max(ymin_v + 2.2, ymax_v))
    end

    axislegend(
        ax_eda,
        [line_total, line_tonic, scatter_er, scatter_sf],
        ["Total SC [µS]", "Tonic SCL + Drift", "ER-SCR Peak", "Spontaneous Fluctuation (SF / NS-SCR)"],
        position = :rt,
        orientation = :horizontal,
        labelsize = lift(f -> min(Float64(f) * 0.90, 11.5), label_font),
        backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.92),
        framecolor = "#cbd5e1",
    )

    # =========================================================================
    # BOTTOM-LEFT: Phasic SCR Habituation across Trials (Barry & Sokolov, 1993)
    # =========================================================================
    ax_hab_scr = Axis(
        fig[3, 1],
        title = "Phasic Habituation — Relative √SCR across Trials (Barry & Sokolov, 1993)",
        xlabel = "Stimulus Trial (S1–S6)",
        ylabel = "√SCR Amplitude [√µS]",
        xticks = (1:n_trials, ["S$k" for k = 1:n_trials]),
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = "#f8fafc",
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
    )
    xlims!(ax_hab_scr, 0.6, n_trials + 0.4)
    ylims!(ax_hab_scr, -0.05, 1.75)

    l_true = lines!(ax_hab_scr, trials_x, lift(d -> d.sqrt_true, sim_data), color = "#0284c7", linewidth = 2.6)
    s_true = scatter!(ax_hab_scr, trials_x, lift(d -> d.sqrt_true, sim_data), color = "#0284c7", markersize = 10)
    l_meas = lines!(ax_hab_scr, trials_x, lift(d -> d.sqrt_meas, sim_data), color = "#e11d48", linestyle = :dash, linewidth = 2.4)
    s_meas = scatter!(ax_hab_scr, trials_x, lift(d -> d.sqrt_meas, sim_data), color = "#e11d48", marker = :rect, markersize = 10)

    axislegend(
        ax_hab_scr,
        [[l_true, s_true], [l_meas, s_meas]],
        ["True Underlying SCR", "Measured (Trough-to-Peak)"],
        position = :rt,
        labelsize = lift(f -> min(Float64(f) * 0.88, 11.0), label_font),
        backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.92),
        framecolor = "#cbd5e1",
    )

    # =========================================================================
    # BOTTOM-RIGHT: Tonic SCL across Trials (Barry & Sokolov, 1993)
    # =========================================================================
    ax_hab_scl = Axis(
        fig[3, 2],
        title = "Tonic Adaptation — Pre-Stimulus SCL across Trials",
        xlabel = "Stimulus Trial (S1–S6)",
        ylabel = "Pre-Stimulus SCL [µS]",
        xticks = (1:n_trials, ["S$k" for k = 1:n_trials]),
        titlesize = title_font,
        xlabelsize = label_font,
        ylabelsize = label_font,
        xticklabelsize = tick_font,
        yticklabelsize = tick_font,
        backgroundcolor = "#f8fafc",
        xgridcolor = "#e2e8f0",
        ygridcolor = "#e2e8f0",
    )
    xlims!(ax_hab_scl, 0.6, n_trials + 0.4)

    lines!(ax_hab_scl, trials_x, lift(d -> d.trial_scl, sim_data), color = "#0f172a", linewidth = 2.6)
    scatter!(ax_hab_scl, trials_x, lift(d -> d.trial_scl, sim_data), color = "#0f172a", markersize = 10)

    on(sim_data) do d
        s_lo, s_hi = extrema(d.trial_scl)
        ylims!(ax_hab_scl, s_lo - 0.6, max(s_lo + 1.5, s_hi + 0.6))
    end
    let (s_lo0, s_hi0) = extrema(sim_data[].trial_scl)
        ylims!(ax_hab_scl, s_lo0 - 0.6, max(s_lo0 + 1.5, s_hi0 + 0.6))
    end

    # =========================================================================
    # CONTROLS PANEL (Row 4)
    # =========================================================================
    ctrl_layout = fig[4, 1:2] = GridLayout(tellwidth = false, halign = :center)
    sl_grid = ctrl_layout[1, 1] = GridLayout(tellwidth = false, halign = :center)

    Label(sl_grid[1, 1], text = "Initial SCR (S1):", font = :bold, fontsize = label_font)
    sl_scr1 = Slider(sl_grid[1, 2], range = 0.0:0.05:2.5, startvalue = 1.20, width = 105)
    Label(sl_grid[1, 3], text = lift(v -> "$(round(v, digits=2)) µS", sl_scr1.value), width = 65, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 4], text = "Habituation:", font = :bold, fontsize = label_font)
    sl_hab = Slider(sl_grid[1, 5], range = 0.0:5.0:100.0, startvalue = 45.0, width = 105)
    Label(sl_grid[1, 6], text = lift(v -> "$(round(Int, v))%", sl_hab.value), width = 48, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 7], text = "Stimulus ISI:", font = :bold, fontsize = label_font)
    sl_isi = Slider(sl_grid[1, 8], range = 4.0:0.5:14.0, startvalue = 12.0, width = 105)
    Label(sl_grid[1, 9], text = lift(v -> "$(round(v, digits=1)) s", sl_isi.value), width = 52, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 10], text = "Tonic Drift:", font = :bold, fontsize = label_font)
    sl_drift = Slider(sl_grid[1, 11], range = -3.0:0.25:2.0, startvalue = -1.0, width = 95)
    Label(sl_grid[1, 12], text = lift(v -> "$(round(v, digits=1)) µS/m", sl_drift.value), width = 72, fontsize = label_font, halign = :left)

    Label(sl_grid[1, 13], text = "SFs/min:", font = :bold, fontsize = label_font)
    sl_sf = Slider(sl_grid[1, 14], range = 0.0:1.0:4.0, startvalue = 1.0, width = 75)
    Label(sl_grid[1, 15], text = lift(v -> "$(round(Int, v))", sl_sf.value), width = 28, fontsize = label_font, halign = :left)

    colgap!(sl_grid, 7)
    rowgap!(ctrl_layout, 12)

    btn_grid = ctrl_layout[2, 1] = GridLayout(tellwidth = false, halign = :center)
    Label(btn_grid[1, 1], text = "Presets:", font = :bold, fontsize = label_font)
    btn_std = Button(btn_grid[1, 2], label = "Standard Habituation", buttoncolor = "#dcfce7", fontsize = btn_font, width = 175, height = 34)
    btn_ovl = Button(
        btn_grid[1, 3],
        label = "Fast ISI Overlap (Superposition)",
        buttoncolor = "#fee2e2",
        fontsize = btn_font,
        width = 235,
        height = 34,
    )
    btn_non = Button(btn_grid[1, 4], label = "Non-Responder (~10%)", buttoncolor = "#e2e8f0", fontsize = btn_font, width = 175, height = 34)
    chk_novel = Checkbox(btn_grid[1, 5], checked = false)
    Label(btn_grid[1, 6], text = "Novel Stimulus at S6 (Dishabituation)", font = :bold, fontsize = label_font)
    colgap!(btn_grid, 10)

    on(val -> (obs_scr1[] = val), sl_scr1.value)
    on(val -> (obs_hab_rate[] = val), sl_hab.value)
    on(val -> (obs_isi[] = val), sl_isi.value)
    on(val -> (obs_drift[] = val), sl_drift.value)
    on(val -> (obs_sf_rate[] = val), sl_sf.value)
    on(val -> (obs_novel_s6[] = val), chk_novel.checked)

    on(btn_std.clicks) do _
        set_close_to!(sl_scr1, 1.20)
        set_close_to!(sl_hab, 45.0)
        set_close_to!(sl_isi, 12.0)
        set_close_to!(sl_drift, -1.0)
        set_close_to!(sl_sf, 1.0)
    end
    on(btn_ovl.clicks) do _
        set_close_to!(sl_scr1, 1.35)
        set_close_to!(sl_hab, 15.0)
        set_close_to!(sl_isi, 4.5)
        set_close_to!(sl_drift, 0.0)
        set_close_to!(sl_sf, 0.0)
    end
    on(btn_non.clicks) do _
        set_close_to!(sl_scr1, 0.0)
        set_close_to!(sl_hab, 50.0)
        set_close_to!(sl_isi, 12.0)
        set_close_to!(sl_drift, -1.25)
        set_close_to!(sl_sf, 0.0)
    end

    rowsize!(fig.layout, 1, Fixed(26))
    rowsize!(fig.layout, 2, Relative(0.45))
    rowsize!(fig.layout, 3, Relative(0.36))
    rowsize!(fig.layout, 4, Fixed(105))
    rowgap!(fig.layout, 10)

    _display_demo_figure(fig)
    return fig, ax_eda, ax_hab_scr, ax_hab_scl
end
