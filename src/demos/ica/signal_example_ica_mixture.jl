"""
    signal_example_ica_mixture()

Interactive Demo — Part 2: What is a Mixture?

A simple visualization proving why component analysis is necessary. It shows
that a single EEG electrode physically records a summated mixture of multiple
underlying brain/artifact sources based on "mixing strengths" (conductance).

## Educational Value
Shows that an electrode does not record a "pure" signal. It records
`1.0 * Brain + W * Blink`. This establishes the forward mixing problem.

## See Also
- [`signal_example_ica_separation`](@ref) — next step: 2 electrodes, mixing matrices, and the ICA algorithm
- [`signal_example_ica_geometry`](@ref) — advanced: 3 sources, rotation geometry, scatter plots
- [`signal_example_ica_clt`](@ref) — background: why mixtures become Gaussian
- [`signal_example_ica_sphering`](@ref) — sphering: morphing data into perfect spheres
- [`signal_example_ica_optimization`](@ref) — optimization: gradient ascent on non-Gaussianity

# Examples
```julia
using EegFun
signal_example_ica_mixture()
```
"""
function signal_example_ica_mixture()
    n = 1000
    Random.seed!(42)

    # ── Styling ───────────────────────────────────────────────────────────────
    col_brain = :blue
    col_blink = :red
    col_mix   = :black

    fig = Figure(size = (1000, 800), figure_padding = (15, 15, 15, 15))

    lbl_sz   = Observable(14)
    tick_sz  = Observable(12)
    title_sz = Observable(16)
    ctrl_sz  = Observable(16)

    ax_args = (
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        leftspinevisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        titlesize = title_sz,
    )

    t_arr = range(0, 4π, length = n)

    # Generate signals
    brain_base = sin.(10.0 .* t_arr)
    # Give the brain signal some spindle-like envelopes
    env = exp.(-((t_arr .- π) .^ 2) ./ 0.5) .+ exp.(-((t_arr .- 3π) .^ 2) ./ 0.5)
    brain_raw = brain_base .* env .+ 0.1 .* randn(n)

    # Isolated blink shape
    blink_raw = 5.0 .* exp.(-((t_arr .- 2π) .^ 2) ./ 0.2)

    # ── Observables ───────────────────────────────────────────────────────────
    mix_obs = Observable(1.0)

    eeg_obs = Observable(zeros(n))

    function update!(mix_w)
        eeg_obs[] = ((1.0 - mix_w) .* brain_raw) .+ (mix_w .* blink_raw)
    end

    # ── UI Layout ─────────────────────────────────────────────────────────────
    Label(fig[0, 1], "What Is a Mixture?  —  Why EEG needs ICA", fontsize = @lift($title_sz + 4), font = :bold, halign = :center)

    # Source 1
    ax_s1 = Axis(fig[1, 1]; title = "Source 1: Brain Oscillation", ax_args...)
    lines!(ax_s1, t_arr, brain_raw; color = col_brain, linewidth = 2.0)
    ylims!(ax_s1, -6, 6)

    # Source 2
    ax_s2 = Axis(fig[2, 1]; title = "Source 2: Artifact (Blink)", ax_args...)
    lines!(ax_s2, t_arr, blink_raw; color = col_blink, linewidth = 2.0)
    ylims!(ax_s2, -6, 6)

    # Mixture (Electrode)
    ax_mix = Axis(fig[3, 1]; title = "The Scalp Electrode Mixture (Brain + Blink)", ax_args...)
    lines!(ax_mix, t_arr, eeg_obs; color = col_mix, linewidth = 2.5)
    ylims!(ax_mix, -10, 10)

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[4, 1]; colgap = 40)
    Label(ctrl[1, 1], "Mix Balance (100% Brain ↔ 100% Artifact):", fontsize = ctrl_sz, font = :bold)
    sl_mix = Slider(ctrl[1, 2], range = 0.0:0.01:1.0, startvalue = 0.5, width = 500)
    Label(ctrl[1, 3], @lift(string(round($(sl_mix.value); digits = 2))), fontsize = ctrl_sz, font = :bold, color = col_mix)

    on(sl_mix.value) do mix_w
        update!(mix_w)
    end

    update!(sl_mix.value[])

    display(fig)
    return fig
end
