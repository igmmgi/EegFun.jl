const _ERP_MEANPOWER = Float64[
    0.0015,
    0.0008,
    0.0006,
    0.0005,
    fill(0.0004, 3)...,
    fill(0.0003, 2)...,
    fill(0.0004, 2)...,
    fill(0.0003, 3)...,
    fill(0.0002, 42)...,
    fill(0.0001, 88)...,
]

function _add_peak_vec!(
    signal::Matrix{Float64},
    trials::Int,
    samples::Int,
    sample_rate::Int,
    peak_freq::Float64,
    peak_amp::Float64,
    peak_latency::Float64,
    jitter_amp::Float64,
    jitter_latency::Float64,
)
    trials == 0 && return signal
    amp = Vector{Float64}(undef, trials)
    @inbounds for tr = 1:trials
        amp[tr] = jitter_amp == 0.0 ? peak_amp : muladd(randn(), jitter_amp, peak_amp)
    end

    phase_scale = 2π * peak_freq / sample_rate
    half_pi = π / 2
    half_width = half_pi / phase_scale

    if jitter_latency == 0.0
        s_lo = max(1, floor(Int, peak_latency - half_width))
        s_hi = min(samples, ceil(Int, peak_latency + half_width))
        @inbounds for s = s_lo:s_hi
            ph = (Float64(s) - peak_latency) * phase_scale
            if abs(ph) <= half_pi
                c = cos(ph)
                for tr = 1:trials
                    signal[tr, s] = muladd(c, amp[tr], signal[tr, s])
                end
            end
        end
    else
        pos = Vector{Float64}(undef, trials)
        min_pos = Inf
        max_pos = -Inf
        @inbounds for tr = 1:trials
            p = muladd(randn(), jitter_latency, peak_latency)
            pos[tr] = p
            min_pos = min(min_pos, p)
            max_pos = max(max_pos, p)
        end
        s_lo = max(1, floor(Int, min_pos - half_width))
        s_hi = min(samples, ceil(Int, max_pos + half_width))
        @inbounds for s = s_lo:s_hi
            sf = Float64(s)
            for tr = 1:trials
                ph = (sf - pos[tr]) * phase_scale
                if abs(ph) <= half_pi
                    signal[tr, s] = muladd(cos(ph), amp[tr], signal[tr, s])
                end
            end
        end
    end
    return signal
end

function _peak_vec(
    trials::Int,
    samples::Int,
    sample_rate::Int,
    peak_freq::Float64,
    peak_amp::Float64,
    peak_latency::Float64,
    jitter_amp::Float64,
    jitter_latency::Float64,
)
    out = zeros(Float64, trials, samples)
    return _add_peak_vec!(out, trials, samples, sample_rate, peak_freq, peak_amp, peak_latency, jitter_amp, jitter_latency)
end

function _noise_vec(trials::Int, samples::Int, sample_rate::Int, noise_amp::Float64)
    noise_amp == 0.0 && return zeros(Float64, trials, samples)

    n_mp = length(_ERP_MEANPOWER)
    inv_mp1 = 1.0 / _ERP_MEANPOWER[1]
    dt_2pi = 2π / sample_rate
    sumsig = 50

    out = Matrix{Float64}(undef, trials, samples)
    col_buf = Vector{Float64}(undef, samples)

    @inbounds for trial = 1:trials
        fill!(col_buf, 0.0)
        f_acc = 0.0
        for _ = 1:sumsig
            f_acc += rand() * 4.0
            idx = min(max(1, ceil(Int, f_acc)), n_mp)
            a = _ERP_MEANPOWER[idx] * inv_mp1
            p = rand() * 2π
            ω = dt_2pi * f_acc
            sin_ω, cos_ω = sincos(ω)
            s_cur, c_cur = sincos(p)
            for s = 1:samples
                s_next = muladd(s_cur, cos_ω, c_cur * sin_ω)
                c_next = muladd(c_cur, cos_ω, -s_cur * sin_ω)
                col_buf[s] = muladd(a, s_next, col_buf[s])
                s_cur, c_cur = s_next, c_next
            end
        end
        for s = 1:samples
            out[trial, s] = col_buf[s] * noise_amp
        end
    end
    return out
end

function _run_erp_simulation(
    trials::Int,
    comps::Matrix{Float64};
    samp_freq::Int = 1000,
    sig_length::Float64 = 1.0,
    base_length::Float64 = 0.2,
    noise_base::Union{Matrix{Float64},Nothing} = nothing,
)
    n_samples = Int(samp_freq * (sig_length + base_length))
    noise_scale = comps[end, 6]
    if noise_base !== nothing
        my_signal = Matrix{Float64}(undef, trials, n_samples)
        if noise_scale == 0.0
            fill!(my_signal, 0.0)
        else
            @inbounds for s = 1:n_samples, tr = 1:trials
                my_signal[tr, s] = noise_base[tr, s] * noise_scale
            end
        end
    else
        my_signal = _noise_vec(trials, n_samples, samp_freq, noise_scale)
    end

    n_components = size(comps, 1)
    ms_to_samples = samp_freq / 1000.0
    base_offset = base_length * samp_freq
    for comp = 1:n_components
        _add_peak_vec!(
            my_signal,
            trials,
            n_samples,
            samp_freq,
            comps[comp, 1],
            comps[comp, 2],
            comps[comp, 3] * ms_to_samples + base_offset,
            comps[comp, 4],
            comps[comp, 5] * ms_to_samples,
        )
    end

    avg = Vector{Float64}(undef, n_samples)
    if trials > 0
        inv_trials = 1.0 / trials
        @inbounds for s = 1:n_samples
            acc = 0.0
            for tr = 1:trials
                acc += my_signal[tr, s]
            end
            avg[s] = acc * inv_trials
        end
    else
        fill!(avg, 0.0)
    end
    return avg, my_signal
end

mutable struct _SimulatedErpComponent
    active::Observable{Bool}
    freq::Observable{Float64}
    amp::Observable{Float64}
    latency::Observable{Float64}
    jitter_amp::Observable{Float64}
    jitter_lat::Observable{Float64}
end

"""
    simulate_erp()

Interactive ERP simulator for teaching signal averaging and component superposition.

Generates synthetic EEG trials from up to 5 cosine-shaped ERP components (with controllable
frequency, amplitude, latency, and trial-to-trial jitter) plus 1/f background noise, displaying
individual trials (or a 95% confidence interval) alongside the trial-averaged ERP.

# Examples
```julia
using GLMakie, EegFun
EegFun.simulate_erp()
```

# Returns
- `fig::Figure`: The Makie figure object
- `ax::Axis`: The main plot axis
"""
function simulate_erp()
    _check_interactive_demo_backend("ERP Simulator demo") || return nothing

    fig = _eegfun_figure(size = (1400, 1000))
    gl = fig[1:2, 1:3] = GridLayout()
    ax = Axis(gl[1, 1:2], xlabel = "Time (ms)", ylabel = "Amplitude (μV)", title = "ERP Simulation Demo")

    vlines!(ax, [0], color = :black, linestyle = :dash, linewidth = 1.5)
    hlines!(ax, [0], color = (:black, 0.3), linewidth = 1)

    comp_layout = gl[1:2, 3] = GridLayout()
    Label(comp_layout[1, 1:11], "Components", fontsize = 20)

    samp_freq = 1000
    sig_length = 1.0
    base_length = 0.2
    n_samples = Int(samp_freq * (sig_length + base_length))
    time = (collect(1:n_samples) .- (samp_freq * base_length)) .* (1000.0 / samp_freq)

    max_trials = 500
    noise_base = Observable(_noise_vec(max_trials, n_samples, samp_freq, 1.0))

    ctrl_layout = gl[2, 1:2] = GridLayout()

    trials = Observable(1)
    Label(ctrl_layout[1, 1], "Number of Trials", halign = :right)
    sl_trials = Slider(ctrl_layout[1, 2], range = 1:1:500, startvalue = 1)
    Label(ctrl_layout[1, 3], lift(x -> string(Int(x)), trials), halign = :left)
    connect!(trials, sl_trials.value)

    noise_amp = Observable(0.0)
    Label(ctrl_layout[1, 4], "Noise Amplitude", halign = :right)
    sl_noise = Slider(ctrl_layout[1, 5], range = 0.0:1:20.0, startvalue = 0.0)
    Label(ctrl_layout[1, 6], lift(x -> string(round(x, digits = 1)), noise_amp), halign = :left)
    connect!(noise_amp, sl_noise.value)

    btn_regen = Button(ctrl_layout[1, 7], label = "Regenerate", width = 120)

    tb_ci = Toggle(ctrl_layout[1, 8])
    tb_ci.active[] = false
    Label(ctrl_layout[1, 9], "Confidence Interval", halign = :left)
    show_ci = tb_ci.active

    colsize!(ctrl_layout, 2, Auto(1))
    colsize!(ctrl_layout, 5, Auto(1))

    components = _SimulatedErpComponent[]

    function create_component(parent, row)
        active = Observable(false)
        freq = Observable(3.0)
        amp = Observable(5.0)
        latency = Observable(100.0)
        jitter_amp = Observable(0.0)
        jitter_lat = Observable(0.0)

        tb = Toggle(parent[row, 1], halign = :center)
        tb.active[] = false
        connect!(active, tb.active)

        sl_freq = Slider(parent[row, 2], range = 0.1:0.1:20.0, startvalue = 3.0)
        Label(parent[row, 3], lift(x -> string(round(x, digits = 1)), freq), halign = :left)
        connect!(freq, sl_freq.value)

        sl_amp = Slider(parent[row, 4], range = -20.0:1.0:20.0, startvalue = 5.0)
        Label(parent[row, 5], lift(x -> string(round(x, digits = 1)), amp), halign = :left)
        connect!(amp, sl_amp.value)

        sl_latency = Slider(parent[row, 6], range = 0.0:10.0:1000.0, startvalue = 100.0)
        Label(parent[row, 7], lift(x -> string(Int(round(x))), latency), halign = :left)
        connect!(latency, sl_latency.value)

        sl_jitter_amp = Slider(parent[row, 8], range = 0.0:1.0:20.0, startvalue = 0.0)
        Label(parent[row, 9], lift(x -> string(round(x, digits = 1)), jitter_amp), halign = :left)
        connect!(jitter_amp, sl_jitter_amp.value)

        sl_jitter_lat = Slider(parent[row, 10], range = 0.0:1.0:150.0, startvalue = 0.0)
        Label(parent[row, 11], lift(x -> string(round(x, digits = 1)), jitter_lat), halign = :left)
        connect!(jitter_lat, sl_jitter_lat.value)

        return _SimulatedErpComponent(active, freq, amp, latency, jitter_amp, jitter_lat)
    end

    for (i, label) in enumerate(("Active", "Freq (Hz)", "Amp (μV)", "Latency (ms)", "Amp Jitter", "Lat Jitter"))
        Label(comp_layout[2, i == 1 ? 1 : (i - 1) * 2], label, font = :bold)
    end

    for i = 1:5
        push!(components, create_component(comp_layout, i + 2))
    end

    colsize!(comp_layout, 1, Fixed(60))
    for c in (2, 4, 6, 8, 10)
        colsize!(comp_layout, c, Auto(1))
        colsize!(comp_layout, c + 1, Fixed(45))
    end

    trial_plot_data = Observable(Point2f[])
    erp_plot_data = Observable(Point2f[])
    time_obs = Observable(Float64[])
    lower_obs = Observable(Float64[])
    upper_obs = Observable(Float64[])
    show_trials = @lift(!$show_ci)

    band!(ax, time_obs, lower_obs, upper_obs; color = (:blue, 0.2), visible = show_ci)
    lines!(ax, trial_plot_data, color = (:gray, 0.15), visible = show_trials)
    lines!(ax, erp_plot_data, color = :blue, linewidth = 3, label = "Empirical Average")

    function update_plot(new_signals, new_erp)
        n_trials, n_smpl = size(new_signals)
        if n_trials == 0
            trial_plot_data[] = Point2f[]
            erp_plot_data[] = Point2f[]
            time_obs[] = Float64[]
            lower_obs[] = Float64[]
            upper_obs[] = Float64[]
            autolimits!(ax)
            return
        end

        if n_trials > 1
            stds = vec(std(new_signals, dims = 1))
            ci_margin = 1.96 .* stds ./ sqrt(n_trials)
            lower_obs[] = new_erp .- ci_margin
            upper_obs[] = new_erp .+ ci_margin
        else
            lower_obs[] = new_erp
            upper_obs[] = new_erp
        end
        time_obs[] = time

        max_drawn_trials = min(n_trials, 50)
        points = Vector{Point2f}(undef, max_drawn_trials * (n_smpl + 1))
        @inbounds for trial = 1:max_drawn_trials
            offset = (trial - 1) * (n_smpl + 1)
            for s = 1:n_smpl
                points[offset+s] = Point2f(time[s], new_signals[trial, s])
            end
            points[offset+n_smpl+1] = Point2f(NaN32, NaN32)
        end
        trial_plot_data[] = points
        erp_plot_data[] = [Point2f(time[s], new_erp[s]) for s = 1:n_samples]
        autolimits!(ax)
    end

    function update_simulation()
        n_active = count(c -> c.active[], components)
        if n_active == 0
            update_plot(zeros(0, n_samples), zeros(n_samples))
            return
        end

        active_comps = Matrix{Float64}(undef, n_active, 6)
        row = 0
        for c in components
            c.active[] || continue
            row += 1
            active_comps[row, 1] = c.freq[]
            active_comps[row, 2] = c.amp[]
            active_comps[row, 3] = c.latency[]
            active_comps[row, 4] = c.jitter_amp[]
            active_comps[row, 5] = c.jitter_lat[]
            active_comps[row, 6] = 0.0
        end
        active_comps[end, 6] = noise_amp[]

        new_erp, new_signals = _run_erp_simulation(
            trials[],
            active_comps;
            samp_freq = samp_freq,
            sig_length = sig_length,
            base_length = base_length,
            noise_base = noise_base[],
        )
        update_plot(new_signals, new_erp)
    end

    update_trigger = Observable(0)
    throttled_trigger = Makie.Observables.async_latest(update_trigger)
    on(throttled_trigger) do _
        update_simulation()
    end

    for comp in components
        on(comp.active) do _
            update_trigger[] += 1
        end
        for param in (comp.freq, comp.amp, comp.latency, comp.jitter_amp, comp.jitter_lat)
            on(param) do _
                comp.active[] && (update_trigger[] += 1)
            end
        end
    end

    on(trials) do _
        update_trigger[] += 1
    end
    on(noise_amp) do _
        any(c -> c.active[], components) && (update_trigger[] += 1)
    end
    on(btn_regen.clicks) do _
        noise_base[] = _noise_vec(max_trials, n_samples, samp_freq, 1.0)
        any(c -> c.active[], components) && (update_trigger[] += 1)
    end

    update_plot(zeros(0, n_samples), zeros(n_samples))
    _display_demo_figure(fig)
    return fig, ax
end
