"""
    signal_example_tf()

Interactive Time-Frequency Analysis Demo — Method & Resolution Comparison

Demonstrates the **time-frequency resolution trade-off** (Heisenberg uncertainty
principle) and compares three TF analysis methods on a synthetic signal composed
of up to three independent frequency components.

## Methods

| Method | Window shape | Parameter | Key property |
|--------|-------------|-----------|--------------|
| Morlet Wavelet | Gaussian | Cycles | Adaptive: window length = cycles / frequency |
| STFT (Hanning) | Hanning | Window (ms) | Fixed: same window length for all frequencies |
| Multitaper (DPSS) | Slepian tapers | Cycles | Adaptive: multiple tapers reduce variance |

**Note on the parameter slider:** For Morlet and Multitaper the slider sets the
number of wavelet cycles. For STFT it sets a fixed window length using a 30 Hz
reference (e.g. 7 cycles → 233 ms). This makes the comparison explicit: Morlet
adapts the window per frequency, STFT uses one window for all frequencies.

## Plots

| Row | Plot | Description |
|-----|------|-------------|
| 1 | Time Domain | Synthetic signal averaged across trials. Dashed line = t = 0. |
| 2 | Time-Frequency Map | Power heatmap computed with the selected method. |

## Controls

| Control | Range | Description |
|---------|-------|-------------|
| Method | Morlet / STFT / Multitaper | TF decomposition method |
| Cycles / Window | 3–12 | Cycles (Morlet/Multitaper) or window = N/30 Hz (STFT) |
| Noise | 0.5–5 | Amplitude of additive white noise per trial |
| Baseline | None / dB / % Change / … | Baseline correction type (uses `tf_baseline`) |
| Window | −1.0–2.0 s | Baseline time window (IntervalSlider) |
| Freq scale | Linear / Log | Toggle log₁₀ y-axis on the TF map |
| ☐ Comp 1/2/3 | on/off | Enable/disable each component |
| Freq (×3) | 1–40 Hz | Component frequency |
| Amp (×3) | 0–5 | Component amplitude |
| Interval (×3) | −1.0–2.0 s | Component onset–offset window |

## Teaching Notes

The **Heisenberg uncertainty principle** for TF analysis: you cannot have
perfect time *and* frequency resolution simultaneously.

- **Morlet, low cycles**: narrow Gaussian → good time localisation, blurry frequency.
- **Morlet, high cycles**: wide Gaussian → sharp frequency, smeared in time.
- **STFT vs Morlet (same slider value)**: STFT uses the same window for all
  frequencies, so low-frequency bands are blurry (few cycles fit the window).
  Morlet adapts — resolution improves for high frequencies automatically.
- **Multitaper vs Morlet**: similar resolution trade-off but lower variance,
  because multiple orthogonal tapers provide independent spectral estimates.

## See Also

- Cohen, M. X. (2014). *Analyzing Neural Time Series Data*. MIT Press. — Chapters 12–16

# Examples
```julia
using EegFun
EegFun.signal_example_tf()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::NamedTuple`: Named tuple with `:time` and `:tf` axes
"""
function signal_example_tf()

    # ── Fixed simulation parameters ──────────────────────────────────────────
    const_samp_rate = 250
    const_n_trials  = 20
    const_time_win  = [-1.0, 2.0]
    const_freqs_out = collect(1.0:1.0:50.0)
    const_n_freqs   = length(const_freqs_out)
    const_n_samples = round(Int, (const_time_win[2] - const_time_win[1]) * const_samp_rate)
    const_times_vec = range(const_time_win[1], step = 1 / const_samp_rate, length = const_n_samples)

    const_step_out = 0.01
    const_sub_step = round(Int, const_step_out * const_samp_rate)
    const_t_sub    = collect(const_times_vec)[1:const_sub_step:end]  # pre-computed subsampled time axis

    # Reference frequency for converting STFT cycles → window_length:
    # window_length_s = cycles / stft_ref_freq  (e.g. 7 cycles → 233 ms at 30 Hz)
    const_stft_ref_hz = 30.0

    # ── Observables ───────────────────────────────────────────────────────────
    method_val            = Observable(:morlet)
    param_val             = Observable(7.0)   # cycles (Morlet/MT) or window = param/30 Hz in s (STFT)
    noise_val             = Observable(0.5)
    baseline_method_val   = Observable(:none)   # :none, :db, :percent, :zscore, :absolute, :relative, :relchange, :normchange
    baseline_interval_val = Observable((-1.0, 0.0))

    comp1_active  = Observable(true)
    freq1_val     = Observable(10.0)
    amp1_val      = Observable(3.0)
    interval1_val = Observable((0.0, 0.8))

    comp2_active  = Observable(false)
    freq2_val     = Observable(25.0)
    amp2_val      = Observable(3.0)
    interval2_val = Observable((0.5, 1.5))

    comp3_active  = Observable(false)
    freq3_val     = Observable(40.0)
    amp3_val      = Observable(3.0)
    interval3_val = Observable((1.0, 2.0))

    # Dynamic parameter label (updates with method)
    param_lbl_obs = Observable("Cycles: 7")

    # ── Figure & axes ─────────────────────────────────────────────────────────
    fig = Figure(size = (1200, 850))

    ax_time = Axis(
        fig[1, 1],
        title        = "Time-Domain Signal (trial average)",
        xlabel       = "Time (s)",
        ylabel       = "Amplitude",
        xgridvisible = true,
        ygridvisible = true,
    )
    ax_tf = Axis(
        fig[2, 1],
        title        = "Time-Frequency Map",
        xlabel       = "Time (s)",
        ylabel       = "Frequency (Hz)",
        xgridvisible = false,
        ygridvisible = false,
    )

    # Stimulus onset markers
    vlines!(ax_time, [0.0], color = :black, linestyle = :dash, linewidth = 1.5)
    hlines!(ax_time, [0.0], color = (:black, 0.25), linewidth = 1.0)
    vlines!(ax_tf, [0.0], color = :white, linestyle = :dash, linewidth = 1.5)

    rowsize!(fig.layout, 1, Relative(0.20))
    rowsize!(fig.layout, 2, Relative(0.55))

    # ── Plot data ─────────────────────────────────────────────────────────────
    time_signal  = Observable(zeros(const_n_samples))
    hm_times     = Observable(const_t_sub)
    tf_power_mat = Observable(zeros(length(hm_times[]), const_n_freqs))

    lines!(ax_time, collect(const_times_vec), time_signal, color = :blue, linewidth = 2)

    hm = heatmap!(ax_tf, hm_times, const_freqs_out, tf_power_mat, colormap = :jet, interpolate = true)
    colorbar_label = Observable("Power (a.u.)")
    Colorbar(fig[2, 2], hm, label = colorbar_label)

    # ── Helper: label string for the colorbar based on baseline method ────────
    """Return a descriptive colorbar label string for the given baseline method."""
    function _colorbar_label(m::Symbol)
        m == :db && return "Power (dB re baseline)"
        m == :percent && return "Power (% change)"
        m == :zscore && return "Power (z-score)"
        m == :absolute && return "Power (absolute change)"
        m == :relative && return "Power (relative power)"
        m == :relchange && return "Power (relative change)"
        (m == :normchange || m == :vssum) && return "Power (norm. change)"
        return "Power (a.u.)"
    end

    # ── Helper: apply tf_baseline if method is not :none ─────────────────────
    """Apply baseline correction to TF data if the current method is not `:none`."""
    function _maybe_baseline(tf_dat::TimeFreqData)
        m  = baseline_method_val[]
        iv = baseline_interval_val[]
        m == :none && return tf_dat
        return Logging.with_logger(Logging.NullLogger()) do
            tf_baseline(tf_dat, iv; method = m)
        end
    end

    # ── Helper: (n_t × n_f) matrix from STFT/MT TF result ────────────────────
    """Extract a (n_t × n_f) power matrix and time vector from an STFT/Multitaper TF result."""
    function _extract_power(tf_dat, n_f)
        pwr_vec      = tf_dat.data[!, :Ch1]
        times_unique = unique(tf_dat.data.time)
        n_t          = length(times_unique)
        pwr_mat      = reshape(pwr_vec, n_f, n_t)   # (n_f × n_t)
        return times_unique, Matrix(pwr_mat')         # → (n_t × n_f)
    end

    # ── Simulation ────────────────────────────────────────────────────────────
    """Generate synthetic signal, run TF decomposition, and update heatmap."""
    function run_simulation()

        freqs_list = Float64[]
        amps_list  = Float64[]
        times_list = Vector{Float64}[]

        if comp1_active[]
            push!(freqs_list, freq1_val[])
            push!(amps_list, amp1_val[])
            iv = interval1_val[]
            push!(times_list, [iv[1], iv[2]])
        end
        if comp2_active[]
            push!(freqs_list, freq2_val[])
            push!(amps_list, amp2_val[])
            iv = interval2_val[]
            push!(times_list, [iv[1], iv[2]])
        end
        if comp3_active[]
            push!(freqs_list, freq3_val[])
            push!(amps_list, amp3_val[])
            iv = interval3_val[]
            push!(times_list, [iv[1], iv[2]])
        end

        if isempty(freqs_list)
            time_signal[]  = zeros(const_n_samples)
            hm_times[]     = const_t_sub
            tf_power_mat[] = zeros(length(const_t_sub), const_n_freqs)
            return
        end

        times_gen, signal_mat =
            generate_signal(const_n_trials, const_time_win, Float64(const_samp_rate), freqs_list, amps_list, times_list, noise_val[])

        time_signal[] = vec(mean(signal_mat, dims = 2))

        epoch_dat = signal_to_data(times_gen, signal_mat, :Ch1, const_samp_rate)
        m         = method_val[]
        p         = param_val[]

        if m == :morlet
            tf_dat = tf_morlet(epoch_dat; frequencies = const_freqs_out, cycles = p)
            tf_dat = _maybe_baseline(tf_dat)
            pwr_vec = tf_dat.data[!, :Ch1]
            pwr_full = reshape(pwr_vec, const_n_freqs, const_n_samples)  # (n_f × n_t_full)
            pwr_sub = pwr_full[:, 1:const_sub_step:end]
            hm_times[] = const_t_sub
            tf_power_mat[] = Matrix(pwr_sub')

        elseif m == :stft
            win_s = p / const_stft_ref_hz
            tf_dat = tf_stft(epoch_dat; frequencies = const_freqs_out, window_length = win_s, time_steps = const_step_out)
            tf_dat = _maybe_baseline(tf_dat)
            t_out, pwr_mat = _extract_power(tf_dat, const_n_freqs)
            hm_times[] = t_out
            tf_power_mat[] = pwr_mat

        elseif m == :multitaper
            tf_dat = tf_multitaper(epoch_dat; frequencies = const_freqs_out, cycles = p, time_steps = const_step_out)
            tf_dat = _maybe_baseline(tf_dat)
            t_out, pwr_mat = _extract_power(tf_dat, const_n_freqs)
            hm_times[] = t_out
            tf_power_mat[] = pwr_mat
        end

        colorbar_label[] = _colorbar_label(baseline_method_val[])

        autolimits!(ax_time)
        xlims!(ax_tf, const_time_win[1], const_time_win[2])
        ylims!(ax_tf, 1.0, 50.0)
    end

    # ── Connect observables → simulation ──────────────────────────────────────
    for obs in [
        method_val,
        param_val,
        noise_val,
        baseline_method_val,
        baseline_interval_val,
        comp1_active,
        freq1_val,
        amp1_val,
        interval1_val,
        comp2_active,
        freq2_val,
        amp2_val,
        interval2_val,
        comp3_active,
        freq3_val,
        amp3_val,
        interval3_val,
    ]
        on(obs) do _
            run_simulation()
        end
    end

    # ── Control panel ──────────────────────────────────────────────────────────
    ctrl  = fig[3, 1:2] = GridLayout(tellheight = false)
    fsize = 16

    # Row 1: method + param + noise
    global_row = ctrl[1, 1:10] = GridLayout()

    Label(global_row[1, 1], "Method:"; fontsize = fsize, halign = :right)
    menu = Menu(global_row[1, 2], options = ["Morlet Wavelet", "STFT (Hanning)", "Multitaper (DPSS)"], default = "Morlet Wavelet")
    on(menu.selection) do sel
        if sel == "Morlet Wavelet"
            method_val[] = :morlet
        elseif sel == "STFT (Hanning)"
            method_val[] = :stft
        else
            method_val[] = :multitaper
        end
        # Update label
        p = param_val[]
        if method_val[] == :stft
            param_lbl_obs[] = "Window: $(round(Int, p / const_stft_ref_hz * 1000)) ms"
        else
            param_lbl_obs[] = "Cycles: $(Int(p))"
        end
    end

    # Single parameter slider — semantics depend on method
    Label(global_row[1, 3], "  Param:"; fontsize = fsize, halign = :right)
    sl_param = Slider(global_row[1, 4], range = 3.0:1.0:12.0, startvalue = 7.0)
    Label(global_row[1, 5], param_lbl_obs; fontsize = fsize)
    on(sl_param.value) do v
        param_val[] = v
        if method_val[] == :stft
            param_lbl_obs[] = "Window: $(round(Int, v / const_stft_ref_hz * 1000)) ms"
        else
            param_lbl_obs[] = "Cycles: $(Int(v))"
        end
    end

    Label(global_row[1, 6], "  Noise:"; fontsize = fsize, halign = :right)
    sl_noise  = Slider(global_row[1, 7], range = 0.5:0.5:5.0, startvalue = 0.5)
    lbl_noise = Label(global_row[1, 8], "0.5"; fontsize = fsize)
    on(sl_noise.value) do v
        noise_val[] = v
        lbl_noise.text = string(round(v, digits = 1))
    end

    # Row 2: baseline controls
    baseline_row = ctrl[2, 1:10] = GridLayout()
    Label(baseline_row[1, 1], "Baseline:"; fontsize = fsize, halign = :right)
    baseline_menu = Menu(
        baseline_row[1, 2],
        options = ["None", "dB", "% Change", "Z-score", "Absolute", "Relative", "Rel. Change", "Norm. Change"],
        default = "None",
    )
    on(baseline_menu.selection) do sel
        baseline_method_val[] = Dict(
            "None"         => :none,
            "dB"           => :db,
            "% Change"     => :percent,
            "Z-score"      => :zscore,
            "Absolute"     => :absolute,
            "Relative"     => :relative,
            "Rel. Change"  => :relchange,
            "Norm. Change" => :normchange,
        )[sel]
    end

    Label(baseline_row[1, 3], "  Window:"; fontsize = fsize, halign = :right)
    baseline_iv_slider = IntervalSlider(baseline_row[1, 4], range = -1.0:0.05:2.0, startvalues = (-1.0, 0.0))
    baseline_iv_lbl = Label(baseline_row[1, 5], "-1.0 – 0.0 s"; fontsize = fsize)
    on(baseline_iv_slider.interval) do iv
        baseline_interval_val[] = iv
        baseline_iv_lbl.text = "$(round(iv[1], digits=2)) – $(round(iv[2], digits=2)) s"
    end

    Label(baseline_row[1, 6], "  Freq scale:"; fontsize = fsize, halign = :right)
    log_toggle = Toggle(baseline_row[1, 7], active = false)
    Label(baseline_row[1, 8], "Log"; fontsize = fsize, halign = :left)
    on(log_toggle.active) do is_log
        ax_tf.yscale = is_log ? log10 : identity
        ylims!(ax_tf, 1.0, 50.0)
    end

    # Rows 3–6: component grid (shifted down from 2–5 to make room for baseline row)
    comp_grid = ctrl[3:6, 1:10] = GridLayout()

    # ── Column headers (row 1 of comp_grid) ────────────────────────────────
    Label(comp_grid[1, 1], ""; fontsize = fsize - 2, font = :bold, halign = :center)
    Label(comp_grid[1, 2], "Active"; fontsize = fsize - 2, font = :bold, halign = :center)
    Label(comp_grid[1, 3], "Freq (Hz)"; fontsize = fsize - 2, font = :bold, halign = :center)
    Label(comp_grid[1, 4], "Amp"; fontsize = fsize - 2, font = :bold, halign = :center)
    Label(comp_grid[1, 5:6], "Interval (s)"; fontsize = fsize - 2, font = :bold, halign = :center)

    # ── Helper: one component row inside comp_grid ──────────────────────────
    """Create a row of component controls (toggle, freq, amp, interval sliders) inside the grid."""
    function make_comp_row(grid_row, label, active_obs, freq_obs, amp_obs, interval_obs, init_active, init_freq, init_amp, init_interval)

        Label(comp_grid[grid_row, 1], label; fontsize = fsize, halign = :right)

        tb = Toggle(comp_grid[grid_row, 2])
        tb.active[] = init_active
        connect!(active_obs, tb.active)

        sl_f  = Slider(comp_grid[grid_row, 3], range = 1.0:1.0:40.0, startvalue = init_freq)
        lbl_f = Label(comp_grid[grid_row, 3, Bottom()], string(Int(init_freq)); fontsize = fsize - 2)
        on(sl_f.value) do v
            freq_obs[] = v
            lbl_f.text = string(Int(round(v)))
        end

        sl_a  = Slider(comp_grid[grid_row, 4], range = 0.0:0.5:5.0, startvalue = init_amp)
        lbl_a = Label(comp_grid[grid_row, 4, Bottom()], string(init_amp); fontsize = fsize - 2)
        on(sl_a.value) do v
            amp_obs[] = v
            lbl_a.text = string(round(v, digits = 1))
        end

        sl_iv  = IntervalSlider(comp_grid[grid_row, 5:6], range = -1.0:0.1:2.0, startvalues = init_interval)
        lbl_iv = Label(comp_grid[grid_row, 5:6, Bottom()], "$(init_interval[1]) – $(init_interval[2])"; fontsize = fsize - 2)
        on(sl_iv.interval) do iv
            interval_obs[] = iv
            lbl_iv.text = "$(round(iv[1], digits=1)) – $(round(iv[2], digits=1))"
        end
    end

    make_comp_row(2, "Comp 1", comp1_active, freq1_val, amp1_val, interval1_val, true, 10.0, 3.0, (0.0, 0.8))
    make_comp_row(3, "Comp 2", comp2_active, freq2_val, amp2_val, interval2_val, false, 25.0, 3.0, (0.5, 1.5))
    make_comp_row(4, "Comp 3", comp3_active, freq3_val, amp3_val, interval3_val, false, 40.0, 3.0, (1.0, 2.0))

    # Size columns now that content exists
    colsize!(comp_grid, 1, Fixed(70))      # component name
    colsize!(comp_grid, 2, Fixed(55))      # toggle
    colsize!(comp_grid, 3, Relative(0.27)) # Freq slider
    colsize!(comp_grid, 4, Relative(0.20)) # Amp slider
    colsize!(comp_grid, 5, Relative(0.25)) # Interval left half
    colsize!(comp_grid, 6, Relative(0.25)) # Interval right half

    rowsize!(fig.layout, 3, Relative(0.25))

    # ── Initial render ─────────────────────────────────────────────────────────
    run_simulation()

    display(fig)
    return fig, (time = ax_time, tf = ax_tf)
end
