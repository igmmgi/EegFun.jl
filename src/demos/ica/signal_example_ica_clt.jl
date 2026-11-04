"""
    signal_example_ica_clt()

Interactive Central Limit Theorem (CLT) demonstration for Independent Component Analysis (ICA).

Mixes up to 50 independent super-Gaussian (Laplace) sources to show how their sum converges
toward a Gaussian distribution (excess kurtosis → 0), illustrating why ICA recovers sources by
maximizing non-Gaussianity.

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_clt()
```
"""
function signal_example_ica_clt()
    _check_interactive_demo_backend("ICA Central Limit Theorem demo") || return nothing

    n_samples = 100_000
    plot_samples = 150
    max_sources = 50
    rng = Random.Xoshiro(42)

    fig = _eegfun_figure(size = (1600, 800), figure_padding = (15, 15, 15, 15))
    title_sz = Observable(20)
    ctrl_sz = Observable(18)

    function _generate_sources()
        srcs = Vector{Vector{Float64}}(undef, max_sources)
        for k = 1:max_sources
            v = Vector{Float64}(undef, n_samples)
            @inbounds for j = 1:n_samples
                mag = -log(1.0 - rand(rng))
                v[j] = ifelse(rand(rng, Bool), mag, -mag)
            end
            srcs[k] = v
        end
        return srcs
    end

    sources_obs = Observable(_generate_sources())
    n_mix_obs = Observable(1)
    mix_data = Observable(zeros(Float64, n_samples))
    kurt_obs = Observable(0.0)
    mix_buf = Vector{Float64}(undef, n_samples)

    function update!()
        n = n_mix_obs[]
        curr_sources = sources_obs[]

        copyto!(mix_buf, curr_sources[1])
        for i = 2:n
            src_i = curr_sources[i]
            @inbounds @simd for j = 1:n_samples
                mix_buf[j] += src_i[j]
            end
        end

        μ = mean(mix_buf)
        inv_σ = 1.0 / std(mix_buf; mean = μ)
        z_mix = Vector{Float64}(undef, n_samples)
        m4_acc = 0.0
        @inbounds for j = 1:n_samples
            z = (mix_buf[j] - μ) * inv_σ
            z_mix[j] = z
            z2 = z * z
            m4_acc += z2 * z2
        end

        mix_data[] = z_mix
        kurt_obs[] = (m4_acc / n_samples) - 3.0
    end

    Label(fig[1, 1:3], "The Central Limit Theorem", fontsize = @lift($title_sz + 8), font = :bold)

    ctrl_grid = GridLayout(fig[2, 1:3])
    Label(ctrl_grid[1, 1], "Number of Sources Mixed Together (N):", fontsize = ctrl_sz, font = :bold)
    sl_mix = Slider(ctrl_grid[1, 2], range = 1:1:max_sources, startvalue = 1, width = 300)
    Label(ctrl_grid[1, 3], @lift(string($n_mix_obs)), fontsize = @lift($ctrl_sz + 4), font = :bold, color = :gray)
    btn_regen = Button(ctrl_grid[1, 4], label = "Regenerate Sources", width = 200)

    on(sl_mix.value) do val
        n_mix_obs[] = val
        update!()
    end
    on(btn_regen.clicks) do _
        sources_obs[] = _generate_sources()
        update!()
    end

    ax_src = Axis(fig[3, 1], title = "The 50 Pure Sources", titlesize = title_sz)
    x_plot = 1:plot_samples
    for i_val = 1:max_sources
        let i = i_val
            line_col = @lift($n_mix_obs >= i ? :blue : :lightgray)
            offset_seg = @lift(@view($(sources_obs)[i][1:plot_samples]) .+ (i * 15.0))
            lines!(ax_src, x_plot, offset_seg, color = line_col, linewidth = 1.5)
        end
    end
    ylims!(ax_src, 0.0, (max_sources + 1) * 15.0)
    hidedecorations!(ax_src)

    ax_ts = Axis(
        fig[3, 2],
        title = @lift("Sum of " * ($n_mix_obs == 1 ? "1 Source" : "$($n_mix_obs) Sources") * "\n(The Mixture)"),
        titlesize = title_sz,
    )
    lines!(ax_ts, x_plot, @lift($mix_data[1:plot_samples]), color = :gray, linewidth = 2)
    ylims!(ax_ts, -5, 5)
    hidedecorations!(ax_ts)

    ax_hist = Axis(fig[3, 3], title = @lift("Histogram of Mixture\n(Kurtosis: $(round($kurt_obs; digits = 2)))"), titlesize = title_sz)
    hist!(ax_hist, mix_data, bins = range(-4, 4, length = 60), color = (:gray, 0.6), normalization = :pdf)

    x_gauss = range(-4, 4, length = 200)
    y_gauss = @. (1.0 / sqrt(2π)) * exp(-0.5 * x_gauss^2)
    lines!(ax_hist, x_gauss, y_gauss, color = :black, linestyle = :dash, linewidth = 4, label = "Perfect Gaussian")
    xlims!(ax_hist, -4, 4)
    ylims!(ax_hist, 0, 0.7)
    axislegend(ax_hist, position = :rt)

    Label(fig[4, 3], "Gaussian Kurtosis = 0.0\nSub-Gaussian (Flat) < 0.0\nSuper-Gaussian (Pointy) > 0.0", color = :gray)

    colsize!(fig.layout, 1, Relative(0.20))
    colsize!(fig.layout, 2, Relative(0.40))
    colsize!(fig.layout, 3, Relative(0.40))
    rowsize!(fig.layout, 3, Relative(0.65))

    update!()
    _display_demo_figure(fig)
    return fig
end
