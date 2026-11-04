"""
    signal_example_ica_mixture()

Interactive demonstration of linear signal mixing in EEG, showing how a single scalp electrode
records a weighted sum of an underlying brain oscillation and a blink artifact.

## See Also
- [`signal_example_ica_separation`](@ref), [`signal_example_ica_geometry`](@ref),
  [`signal_example_ica_clt`](@ref), [`signal_example_ica_sphering`](@ref),
  [`signal_example_ica_optimization`](@ref)

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_mixture()
```
"""
function signal_example_ica_mixture()
    _check_interactive_demo_backend("ICA Mixture demo") || return nothing

    n = 1000
    rng = Random.Xoshiro(42)

    fig = _eegfun_figure(size = (1000, 800), figure_padding = (15, 15, 15, 15))
    title_sz = Observable(16)
    ctrl_sz = Observable(16)

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

    t_arr = collect(range(0.0, 4π, length = n))
    noise = randn(rng, n)
    brain_raw = @. sin(10.0 * t_arr) * (exp(-((t_arr - π)^2) / 0.5) + exp(-((t_arr - 3π)^2) / 0.5)) + 0.1 * noise
    blink_raw = @. 5.0 * exp(-((t_arr - 2π)^2) / 0.2)

    eeg_obs = Observable(zeros(Float64, n))
    function update!(mix_w)
        eeg_obs[] = @. (1.0 - mix_w) * brain_raw + mix_w * blink_raw
    end

    Label(fig[0, 1], "What Is a Mixture?  —  Why EEG needs ICA", fontsize = @lift($title_sz + 4), font = :bold, halign = :center)

    ax_s1 = Axis(fig[1, 1]; title = "Source 1: Brain Oscillation", ax_args...)
    lines!(ax_s1, t_arr, brain_raw; color = :blue, linewidth = 2.0)
    ylims!(ax_s1, -6, 6)

    ax_s2 = Axis(fig[2, 1]; title = "Source 2: Artifact (Blink)", ax_args...)
    lines!(ax_s2, t_arr, blink_raw; color = :red, linewidth = 2.0)
    ylims!(ax_s2, -6, 6)

    ax_mix = Axis(fig[3, 1]; title = "The Scalp Electrode Mixture (Brain + Blink)", ax_args...)
    lines!(ax_mix, t_arr, eeg_obs; color = :black, linewidth = 2.5)
    ylims!(ax_mix, -10, 10)

    ctrl = GridLayout(fig[4, 1]; colgap = 40)
    Label(ctrl[1, 1], "Mix Balance (100% Brain ↔ 100% Artifact):", fontsize = ctrl_sz, font = :bold)
    sl_mix = Slider(ctrl[1, 2], range = 0.0:0.01:1.0, startvalue = 0.5, width = 500)
    Label(ctrl[1, 3], @lift(string(round($(sl_mix.value); digits = 2))), fontsize = ctrl_sz, font = :bold, color = :black)

    on(update!, sl_mix.value)
    update!(sl_mix.value[])

    _display_demo_figure(fig)
    return fig
end
