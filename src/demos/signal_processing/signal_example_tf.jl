"""
    signal_example_tf()

Interactive time-frequency analysis demo comparing Morlet wavelets, STFT (Hanning), and
Multitaper (DPSS) decomposition across up to three synthetic oscillatory components, with
live resolution, noise, frequency-scale, and baseline-correction controls.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_tf()
```

# Returns
- `fig::Figure`: The Makie figure object
- `axes::NamedTuple`: Named tuple with `:time` and `:tf` axes
"""
function signal_example_tf()
    _check_interactive_demo_backend("Signal Time-Frequency demo") || return nothing

    const_samp_rate = 250
    const_n_trials = 20
    const_time_win = [-1.0, 2.0]
    const_freqs_out = collect(1.0:1.0:50.0)
    const_n_freqs = length(const_freqs_out)
    const_n_samples = round(Int, (const_time_win[2] - const_time_win[1]) * const_samp_rate)
    const_times_vec = collect(range(const_time_win[1], step = 1 / const_samp_rate, length = const_n_samples))

    const_step_out = 0.01
    const_sub_step = round(Int, const_step_out * const_samp_rate)
    const_t_sub = const_times_vec[1:const_sub_step:end]
    const_stft_ref_hz = 30.0

    method_val = Observable(:morlet)
    param_val = Observable(7.0)
    noise_val = Observable(0.5)
    baseline_method_val = Observable(:none)
    baseline_interval_val = Observable((-1.0, 0.0))

    comp1_active = Observable(true)
    freq1_val = Observable(10.0)
    amp1_val = Observable(3.0)
    interval1_val = Observable((0.0, 0.8))

    comp2_active = Observable(false)
    freq2_val = Observable(25.0)
    amp2_val = Observable(3.0)
    interval2_val = Observable((0.5, 1.5))

    comp3_active = Observable(false)
    freq3_val = Observable(40.0)
    amp3_val = Observable(3.0)
    interval3_val = Observable((1.0, 2.0))

    param_lbl_obs = Observable("Cycles: 7")

    fig = _eegfun_figure(size = (1350, 850))
    title_sz, lbl_sz, tick_sz, ctrl_sz, hdr_sz = _demo_adaptive_fonts(fig, 1350, (18, 13), (15, 11), (12, 10), (15, 11), (13, 10))

    ax_time = Axis(
        fig[1, 1],
        title = "Time-Domain Signal (trial average)",
        xlabel = "Time (s)",
        ylabel = "Amplitude",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible = true,
        ygridvisible = true,
    )
    ax_tf = Axis(
        fig[2, 1],
        title = "Time-Frequency Map",
        xlabel = "Time (s)",
        ylabel = "Frequency (Hz)",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        xticklabelsize = tick_sz,
        yticklabelsize = tick_sz,
        xgridvisible = false,
        ygridvisible = false,
    )

    vlines!(ax_time, [0.0], color = :black, linestyle = :dash, linewidth = 1.5)
    hlines!(ax_time, [0.0], color = (:black, 0.25), linewidth = 1.0)
    vlines!(ax_tf, [0.0], color = :white, linestyle = :dash, linewidth = 1.5)

    rowsize!(fig.layout, 1, Relative(0.20))
    rowsize!(fig.layout, 2, Relative(0.55))

    time_signal = Observable(zeros(const_n_samples))
    hm_times = Observable(const_t_sub)
    hm_freqs = Observable(const_freqs_out)
    tf_power_mat = Observable(zeros(length(hm_times[]), const_n_freqs))

    lines!(ax_time, const_times_vec, time_signal, color = :blue, linewidth = 2)
    hm = heatmap!(ax_tf, hm_times, hm_freqs, tf_power_mat, colormap = :jet, interpolate = true)
    colorbar_label = Observable("Power (a.u.)")
    Colorbar(fig[2, 2], hm, label = colorbar_label, labelsize = lbl_sz, ticklabelsize = tick_sz)

    function _update_tf_heatmap!(t_vec_new::Vector{Float64}, f_vec_new::Vector{Float64}, pwr_mat_new::Matrix{Float64})
        hm_times.val = t_vec_new
        hm_freqs.val = f_vec_new
        tf_power_mat.val = pwr_mat_new
        notify(hm_times)
        notify(hm_freqs)
        notify(tf_power_mat)
    end

    function _colorbar_label(m::Symbol)
        m === :db && return "Power (dB re baseline)"
        m === :percent && return "Power (% change)"
        m === :zscore && return "Power (z-score)"
        m === :absolute && return "Power (absolute change)"
        m === :relative && return "Power (relative power)"
        m === :relchange && return "Power (relative change)"
        (m === :normchange || m === :vssum) && return "Power (norm. change)"
        return "Power (a.u.)"
    end

    function _maybe_baseline(tf_dat::TimeFreqData)
        m = baseline_method_val[]
        m === :none && return tf_dat
        return Logging.with_logger(Logging.NullLogger()) do
            tf_baseline(tf_dat, baseline_interval_val[]; method = m)
        end
    end

    function _transpose_power_col(col::AbstractVector, n_f::Int, step::Int = 1)
        n_t_full = length(col) ÷ n_f
        t_indices = 1:step:n_t_full
        n_t = length(t_indices)
        out = Matrix{Float64}(undef, n_t, n_f)
        @inbounds for (ti_out, ti_src) in enumerate(t_indices)
            offset = (ti_src - 1) * n_f
            for fi = 1:n_f
                out[ti_out, fi] = Float64(col[offset+fi])
            end
        end
        return out
    end

    function _extract_power(tf_dat::TimeFreqData, n_f::Int)
        df = tf_dat.data_power
        times_unique = unique(df.time)
        return Float64.(times_unique), _transpose_power_col(df[!, :Ch1], n_f)
    end

    function run_simulation()
        freqs_list = Float64[]
        amps_list = Float64[]
        times_list = Vector{Float64}[]

        for (act, fr, am, iv_obs) in (
            (comp1_active, freq1_val, amp1_val, interval1_val),
            (comp2_active, freq2_val, amp2_val, interval2_val),
            (comp3_active, freq3_val, amp3_val, interval3_val),
        )
            act[] || continue
            push!(freqs_list, fr[])
            push!(amps_list, am[])
            iv = iv_obs[]
            push!(times_list, [iv[1], iv[2]])
        end

        if isempty(freqs_list)
            time_signal[] = zeros(const_n_samples)
            _update_tf_heatmap!(const_t_sub, const_freqs_out, zeros(length(const_t_sub), const_n_freqs))
            return
        end

        times_gen, signal_mat =
            generate_signal(const_n_trials, const_time_win, Float64(const_samp_rate), freqs_list, amps_list, times_list, noise_val[])

        time_signal[] = vec(mean(signal_mat, dims = 2))

        epoch_dat = signal_to_data(times_gen, signal_mat, :Ch1, const_samp_rate)
        m = method_val[]
        p = param_val[]

        if m === :morlet
            tf_dat = _maybe_baseline(tf_morlet(epoch_dat; frequencies = const_freqs_out, cycles = p))
            pwr_mat = _transpose_power_col(tf_dat.data_power[!, :Ch1], const_n_freqs, const_sub_step)
            _update_tf_heatmap!(const_t_sub, const_freqs_out, pwr_mat)
        elseif m === :stft
            win_s = p / const_stft_ref_hz
            tf_dat = _maybe_baseline(tf_stft(epoch_dat; frequencies = const_freqs_out, window_length = win_s, time_steps = const_step_out))
            t_out, pwr_mat = _extract_power(tf_dat, const_n_freqs)
            _update_tf_heatmap!(t_out, const_freqs_out, pwr_mat)
        elseif m === :multitaper
            tf_dat = _maybe_baseline(tf_multitaper(epoch_dat; frequencies = const_freqs_out, cycles = p, time_steps = const_step_out))
            t_out, pwr_mat = _extract_power(tf_dat, const_n_freqs)
            _update_tf_heatmap!(t_out, const_freqs_out, pwr_mat)
        end

        colorbar_label[] = _colorbar_label(baseline_method_val[])
        autolimits!(ax_time)
        xlims!(ax_tf, const_time_win[1], const_time_win[2])
        ylims!(ax_tf, 1.0, 50.0)
    end

    is_batch_updating = Ref(false)
    for obs in (
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
    )
        on(obs) do _
            is_batch_updating[] || run_simulation()
        end
    end

    ctrl = fig[3, 1:2] = GridLayout(tellheight = false)
    top_ctrl = ctrl[1, 1] = GridLayout()

    Label(top_ctrl[1, 1], "Method:"; fontsize = ctrl_sz, halign = :right)
    menu = Menu(top_ctrl[1, 2], options = ["Morlet Wavelet", "STFT (Hanning)", "Multitaper (DPSS)"], default = "Morlet Wavelet")
    on(menu.selection) do sel
        method_val[] = sel == "Morlet Wavelet" ? :morlet : sel == "STFT (Hanning)" ? :stft : :multitaper
        p = param_val[]
        param_lbl_obs[] = method_val[] === :stft ? "Window: $(round(Int, p / const_stft_ref_hz * 1000)) ms" : "Cycles: $(Int(p))"
    end

    Label(top_ctrl[1, 3], "  Param:"; fontsize = ctrl_sz, halign = :right)
    sl_param = Slider(top_ctrl[1, 4], range = 3.0:1.0:12.0, startvalue = 7.0)
    Label(top_ctrl[1, 5], param_lbl_obs; fontsize = ctrl_sz, halign = :left)
    on(sl_param.value) do v
        param_val[] = v
        param_lbl_obs[] = method_val[] === :stft ? "Window: $(round(Int, v / const_stft_ref_hz * 1000)) ms" : "Cycles: $(Int(v))"
    end

    Label(top_ctrl[1, 6], "  Noise:"; fontsize = ctrl_sz, halign = :right)
    sl_noise = Slider(top_ctrl[1, 7], range = 0.5:0.5:5.0, startvalue = 0.5)
    lbl_noise = Label(top_ctrl[1, 8], "0.5"; fontsize = ctrl_sz, halign = :left)
    on(sl_noise.value) do v
        noise_val[] = v
        lbl_noise.text = string(round(v, digits = 1))
    end

    Label(top_ctrl[2, 1], "Baseline:"; fontsize = ctrl_sz, halign = :right)
    baseline_menu = Menu(
        top_ctrl[2, 2],
        options = ["None", "dB", "% Change", "Z-score", "Absolute", "Relative", "Rel. Change", "Norm. Change"],
        default = "None",
    )
    on(baseline_menu.selection) do sel
        baseline_method_val[] =
            sel == "dB" ? :db :
            sel == "% Change" ? :percent :
            sel == "Z-score" ? :zscore :
            sel == "Absolute" ? :absolute :
            sel == "Relative" ? :relative : sel == "Rel. Change" ? :relchange : sel == "Norm. Change" ? :normchange : :none
    end

    Label(top_ctrl[2, 3], "  Window:"; fontsize = ctrl_sz, halign = :right)
    baseline_iv_slider = IntervalSlider(top_ctrl[2, 4], range = -1.0:0.05:2.0, startvalues = (-1.0, 0.0))
    baseline_iv_lbl = Label(top_ctrl[2, 5], "-1.0 – 0.0 s"; fontsize = ctrl_sz, halign = :left)
    on(baseline_iv_slider.interval) do iv
        baseline_interval_val[] = iv
        baseline_iv_lbl.text = "$(round(iv[1], digits=2)) – $(round(iv[2], digits=2)) s"
    end

    Label(top_ctrl[2, 6], "  Freq scale:"; fontsize = ctrl_sz, halign = :right)
    log_toggle = Toggle(top_ctrl[2, 7], active = false)
    Label(top_ctrl[2, 8], "Log"; fontsize = ctrl_sz, halign = :left)
    on(log_toggle.active) do is_log
        ax_tf.yscale = is_log ? log10 : identity
        ylims!(ax_tf, 1.0, 50.0)
    end

    colsize!(top_ctrl, 1, Fixed(80))
    colsize!(top_ctrl, 2, Auto(1))
    colsize!(top_ctrl, 3, Fixed(80))
    colsize!(top_ctrl, 4, Auto(1))
    colsize!(top_ctrl, 5, Fixed(130))
    colsize!(top_ctrl, 6, Fixed(100))
    colsize!(top_ctrl, 7, Auto(1))
    colsize!(top_ctrl, 8, Fixed(50))

    comp_grid = ctrl[2, 1] = GridLayout()
    Label(comp_grid[1, 1], ""; fontsize = hdr_sz, font = :bold, halign = :center)
    Label(comp_grid[1, 2], "Active"; fontsize = hdr_sz, font = :bold, halign = :center)
    Label(comp_grid[1, 3:4], "Freq (Hz)"; fontsize = hdr_sz, font = :bold, halign = :center)
    Label(comp_grid[1, 5:6], "Amp"; fontsize = hdr_sz, font = :bold, halign = :center)
    Label(comp_grid[1, 7:8], "Interval (s)"; fontsize = hdr_sz, font = :bold, halign = :center)

    function make_comp_row(grid_row, label, active_obs, freq_obs, amp_obs, interval_obs, init_active, init_freq, init_amp, init_interval)
        Label(comp_grid[grid_row, 1], label; fontsize = ctrl_sz, halign = :right)

        tb = Toggle(comp_grid[grid_row, 2], halign = :center)
        tb.active[] = init_active
        connect!(active_obs, tb.active)

        sl_f = Slider(comp_grid[grid_row, 3], range = 1.0:1.0:40.0, startvalue = init_freq)
        lbl_f = Label(comp_grid[grid_row, 4], string(Int(init_freq)); fontsize = ctrl_sz, halign = :left)
        on(sl_f.value) do v
            freq_obs[] = v
            lbl_f.text = string(Int(round(v)))
        end

        sl_a = Slider(comp_grid[grid_row, 5], range = 0.0:0.5:5.0, startvalue = init_amp)
        lbl_a = Label(comp_grid[grid_row, 6], string(init_amp); fontsize = ctrl_sz, halign = :left)
        on(sl_a.value) do v
            amp_obs[] = v
            lbl_a.text = string(round(v, digits = 1))
        end

        sl_iv = IntervalSlider(comp_grid[grid_row, 7], range = -1.0:0.1:2.0, startvalues = init_interval)
        lbl_iv = Label(comp_grid[grid_row, 8], "$(init_interval[1]) – $(init_interval[2])"; fontsize = ctrl_sz, halign = :left)
        on(sl_iv.interval) do iv
            interval_obs[] = iv
            lbl_iv.text = "$(round(iv[1], digits=1)) – $(round(iv[2], digits=1))"
        end
    end

    is_batch_updating[] = true
    try
        make_comp_row(2, "Comp 1", comp1_active, freq1_val, amp1_val, interval1_val, true, 10.0, 3.0, (0.0, 0.8))
        make_comp_row(3, "Comp 2", comp2_active, freq2_val, amp2_val, interval2_val, false, 25.0, 3.0, (0.5, 1.5))
        make_comp_row(4, "Comp 3", comp3_active, freq3_val, amp3_val, interval3_val, false, 40.0, 3.0, (1.0, 2.0))
    finally
        is_batch_updating[] = false
    end

    colsize!(comp_grid, 1, Fixed(70))
    colsize!(comp_grid, 2, Fixed(55))
    colsize!(comp_grid, 3, Auto(1))
    colsize!(comp_grid, 4, Fixed(35))
    colsize!(comp_grid, 5, Auto(1))
    colsize!(comp_grid, 6, Fixed(40))
    colsize!(comp_grid, 7, Auto(1.5))
    colsize!(comp_grid, 8, Fixed(100))

    rowsize!(fig.layout, 3, Relative(0.25))

    run_simulation()
    _display_demo_figure(fig)
    return fig, (time = ax_time, tf = ax_tf)
end
