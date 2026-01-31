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
    # Vectorized trials handling
    rand_lat = randn(trials) .* jitter_latency
    pos = peak_latency .+ rand_lat

    # Create time matrix (trials x samples)
    t = (1:samples)'
    # Calculate phase for all trials and samples at once
    phase = (t .- pos) ./ sample_rate .* 2π .* peak_freq

    # Apply cosine and mask
    signal = cos.(phase)
    # Mask values outside [-pi/2, pi/2] phase range
    signal[abs.(phase) .> π/2] .= 0

    # Apply amplitude with jitter
    rand_amp = randn(trials) .* jitter_amp
    signal .*= (peak_amp .+ rand_amp)

    return signal
end

function _noise_vec(trials::Int, samples::Int, sample_rate::Int, noise_amp::Float64)

    noise_amp == 0 && return zeros(trials, samples)

    # Mean power spectrum of human EEG (hardcoded as in original MATLAB)
    meanpower = [
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

    sumsig = 50  # number of sinusoids
    signal = zeros(trials, samples)
    t_vec = (1:samples) ./ sample_rate .* 2π

    for trial = 1:trials
        freqs = cumsum(rand(sumsig) .* 4)
        indices = min.(ceil.(Int, freqs), length(meanpower))
        amps = meanpower[indices] ./ meanpower[1]
        phases = rand(sumsig) .* 2π
        for (f, a, p) in zip(freqs, amps, phases)
            signal[trial, :] .+= sin.(t_vec .* f .+ p) .* a
        end
    end

    return signal .* noise_amp
end

function _simulate_my_erp(trials::Int, comps::Matrix{Float64}; samp_freq::Int = 1000, sig_length::Float64 = 1.0, base_length::Float64 = 0.2)
    n_samples = Int(samp_freq * (sig_length + base_length))
    my_signal = _noise_vec(trials, n_samples, samp_freq, comps[end, 6])

    n_components = size(comps, 1)
    for comp = 1:n_components
        my_signal .+= _peak_vec(
            trials,
            n_samples,
            samp_freq,
            comps[comp, 1],  # peakFreq
            comps[comp, 2],  # peakAmp
            comps[comp, 3] + (base_length * samp_freq),  # peakLatency
            comps[comp, 4],  # jitterAmp
            comps[comp, 5],  # jitterLatency
        )
    end

    return vec(mean(my_signal, dims = 1)), my_signal
end

# Component parameters structure
mutable struct Component
    active::Observable{Bool}
    freq::Observable{Float64}
    amp::Observable{Float64}
    latency::Observable{Float64}
    jitter_amp::Observable{Float64}
    jitter_lat::Observable{Float64}
    controls::Vector{Any}
end

function simulate_erp()
    # Create the figure with more space for component controls
    fig = Figure(size = (1400, 1000))

    # Create layout
    gl = fig[1:2, 1:3] = GridLayout()

    # Create the main plot axis
    ax = Axis(gl[1, 1:2], xlabel = "Time (ms)", ylabel = "Amplitude (μV)", title = "ERP Simulation Demo")

    # Component controls layout (right side)
    comp_layout = gl[1:2, 3] = GridLayout()
    Label(comp_layout[1, 1], "Components", fontsize = 20)

    # Calculate time points based on parameters
    samp_freq = 1000
    sig_length = 1.0
    base_length = 0.2
    n_samples = Int(samp_freq * (sig_length + base_length))

    # Create observables for the plot data
    time = Observable(collect(1:n_samples) .- (samp_freq * base_length))

    # Trials control at the top
    trials = Observable(1)
    sl_trials = Slider(gl[2, 1], range = 1:1:500, startvalue = 1)
    Label(gl[2, 1, Top()], "Number of Trials")
    Label(gl[2, 1, Bottom()], lift(x -> string(Int(x)), trials))
    connect!(trials, sl_trials.value)

    # Noise control
    noise_amp = Observable(0.0)
    sl_noise = Slider(gl[2, 2], range = 0.0:1:20.0, startvalue = 0.0)
    Label(gl[2, 2, Top()], "Noise Amplitude")
    Label(gl[2, 2, Bottom()], lift(x -> string(round(x, digits = 1)), noise_amp))
    connect!(noise_amp, sl_noise.value)

    # List to store components
    components = Component[]

    function create_component(parent, row)
        # Initialize all observables - start INACTIVE
        active = Observable(false)
        freq = Observable(3.0)
        amp = Observable(5.0)
        latency = Observable(100.0)
        jitter_amp = Observable(0.0)
        jitter_lat = Observable(0.0)

        # Component controls
        controls = []

        # Active toggle - start INACTIVE
        tb = Toggle(parent[row, 1])
        tb.active[] = false
        active[] = false
        connect!(active, tb.active)
        push!(controls, tb)

        # Parameter sliders with value labels
        sl_freq = Slider(parent[row, 2], range = 0.1:0.1:5.0, startvalue = 3.0)
        lbl_freq = Label(parent[row, 2, Bottom()], lift(x -> string(round(x, digits = 1)), freq))
        connect!(freq, sl_freq.value)
        push!(controls, sl_freq)
        push!(controls, lbl_freq)

        sl_amp = Slider(parent[row, 3], range = -10.0:1.0:10.0, startvalue = 5.0)
        lbl_amp = Label(parent[row, 3, Bottom()], lift(x -> string(round(x, digits = 1)), amp))
        connect!(amp, sl_amp.value)
        push!(controls, sl_amp)
        push!(controls, lbl_amp)

        sl_latency = Slider(parent[row, 4], range = 0.0:10.0:1000.0, startvalue = 100.0)
        lbl_latency = Label(parent[row, 4, Bottom()], lift(x -> string(round(x, digits = 0)), latency))
        connect!(latency, sl_latency.value)
        push!(controls, sl_latency)
        push!(controls, lbl_latency)

        sl_jitter_amp = Slider(parent[row, 5], range = 0.0:1.0:20.0, startvalue = 0.0)
        lbl_jitter_amp = Label(parent[row, 5, Bottom()], lift(x -> string(round(x, digits = 1)), jitter_amp))
        connect!(jitter_amp, sl_jitter_amp.value)
        push!(controls, sl_jitter_amp)
        push!(controls, lbl_jitter_amp)

        sl_jitter_lat = Slider(parent[row, 6], range = 0.0:1.0:50.0, startvalue = 0.0)
        lbl_jitter_lat = Label(parent[row, 6, Bottom()], lift(x -> string(round(x, digits = 1)), jitter_lat))
        connect!(jitter_lat, sl_jitter_lat.value)
        push!(controls, sl_jitter_lat)
        push!(controls, lbl_jitter_lat)

        Component(active, freq, amp, latency, jitter_amp, jitter_lat, controls)
    end

    # Component parameter labels
    labels = ["Active", "Freq (Hz)", "Amp (μV)", "Latency (ms)", "Amp Jitter", "Lat Jitter"]
    for (i, label) in enumerate(labels)
        Label(comp_layout[2, i], label)
    end

    # Add initial components - all start INACTIVE
    for i = 1:5
        comp = create_component(comp_layout, i + 2)
        push!(components, comp)
    end

    # Update function
    function update_simulation()
        # Count active components
        n_active = count(c -> c.active[], components)

        if n_active == 0 # No active components - set everything to zero
            update_plot(zeros(0, n_samples), zeros(n_samples))
        else # Collect parameters from active components only
            active_comps_list = [[c.freq[] c.amp[] c.latency[] c.jitter_amp[] c.jitter_lat[] 0.0] for c in components if c.active[]]

            if isempty(active_comps_list)
                update_plot(zeros(trials[], n_samples), zeros(n_samples))
                return
            end

            active_comps = vcat(active_comps_list...)
            active_comps[end, 6] = noise_amp[]

            new_erp, new_signals =
                _simulate_my_erp(trials[], active_comps, samp_freq = samp_freq, sig_length = sig_length, base_length = base_length)

            # Update plot directly without clearing axis
            update_plot(new_signals, new_erp)
        end
    end

    # Create plot objects
    trial_plot_data = Observable(Point2f[])
    erp_plot_data = Observable(Point2f[])

    lines!(ax, trial_plot_data, color = (:gray, 0.2))
    lines!(ax, erp_plot_data, color = :blue, linewidth = 3)

    # Plot update function
    function update_plot(new_signals, new_erp)
        t = time[]
        n_trials, n_samples = size(new_signals)

        if n_trials == 0
            trial_plot_data[] = Point2f[]
            erp_plot_data[] = Point2f[]
            return
        end

        # Efficiently prepare trial data with NaN separators
        # We transform the matrix into a single vector of points
        points = Vector{Point2f}(undef, n_trials * (n_samples + 1))
        for trial = 1:n_trials
            offset = (trial - 1) * (n_samples + 1)
            for s = 1:n_samples
                points[offset+s] = Point2f(t[s], new_signals[trial, s])
            end
            points[offset+n_samples+1] = Point2f(NaN32, NaN32)
        end
        trial_plot_data[] = points

        # Prepare ERP data
        erp_plot_data[] = [Point2f(t[s], new_erp[s]) for s = 1:n_samples]

        autolimits!(ax)
    end

    # Connect parameter changes to simulation update
    for comp in components
        # Active state always triggers update
        on(comp.active) do _
            update_simulation()
        end

        # Other parameters only trigger update if component is active
        for param in [comp.freq, comp.amp, comp.latency, comp.jitter_amp, comp.jitter_lat]
            on(param) do _
                if comp.active[]
                    update_simulation()
                end
            end
        end
    end

    on(trials) do _
        update_simulation()
    end

    on(noise_amp) do _
        if any(c -> c.active[], components)
            update_simulation()
        end
    end

    # Initialize with empty plot
    update_plot(zeros(0, n_samples), zeros(n_samples))

    display(fig)
    return fig, ax

end
