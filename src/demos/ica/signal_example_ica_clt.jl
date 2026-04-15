"""
    signal_example_ica_clt()

Interactive Demo — The Central Limit Theorem (CLT) in ICA.

Visually proves the core statistical engine of Independent Component Analysis.
The CLT states that the sum of multiple independent random variables approaches 
a Gaussian (normal bell-curve) distribution, no matter what their original shapes were.

Since a scalp electrode records a mixture of hundreds of independent brain, 
muscle, and eye sources, it effectively sums them up, creating a highly 
Gaussian recorded signal.

Because of this, ICA can recover the pure original sources simply by searching 
for the unmixing weights that make the signal look the *least* Gaussian 
(e.g., maximizing kurtosis).

# Examples
```julia
using EegFun
EegFun.signal_example_ica_clt()
```
"""
function signal_example_ica_clt()
    N_samples = 100_000
    Plot_samples = 150 # Drastically reduced to prevent visual spaghetti of white noise
    Max_Sources = 50
    Random.seed!(42)

    fig = Figure(size = (1600, 800), figure_padding = (15, 15, 15, 15))

    lbl_sz   = Observable(16)
    title_sz = Observable(20)
    ctrl_sz  = Observable(18)

    # ── Source Generation ─────────────────────────────────────────────────────
    # Generate 50 independent, highly Non-Gaussian sources.
    # We use a sparse, "pointy" Laplace distribution (typical of brain artifacts)
    # because it takes many more sums to converge to a Gaussian than a flat box.
    # A single source is a massive sharp peak with long tails.
    sources_obs = Observable([(-log.(rand(N_samples))) .* sign.(randn(N_samples)) for _ = 1:Max_Sources])

    function calc_kurtosis(x)
        v = var(x; corrected = false)
        E_x4 = mean((x .- mean(x)) .^ 4)
        return (E_x4 / (v^2)) - 3.0
    end

    # ── Observables ───────────────────────────────────────────────────────────
    n_mix_obs = Observable(1)
    mix_data = Observable(zeros(N_samples))
    kurt_obs = Observable(0.0)

    function update!()
        n = n_mix_obs[]
        curr_sources = sources_obs[]

        # Sum the first `n` independent sources together
        mix_sum = zeros(N_samples)
        for i = 1:n
            mix_sum .+= curr_sources[i]
        end

        # Standardize so the histogram stays visually scaled correctly
        z_mix = (mix_sum .- mean(mix_sum)) ./ std(mix_sum)

        mix_data[] = z_mix
        kurt_obs[] = calc_kurtosis(z_mix)
        notify(n_mix_obs) # trigger color updates
    end

    # ── Layout: UI Elements ───────────────────────────────────────────────────
    Label(fig[1, 1:3], "The Central Limit Theorem", fontsize = @lift($title_sz + 8), font = :bold)

    # Controls
    ctrl_grid = GridLayout(fig[2, 1:3])
    Label(ctrl_grid[1, 1], "Number of Sources Mixed Together (N):", fontsize = ctrl_sz, font = :bold)
    sl_mix = Slider(ctrl_grid[1, 2], range = 1:1:Max_Sources, startvalue = 1, width = 300)
    Label(ctrl_grid[1, 3], @lift(string($(n_mix_obs))), fontsize = @lift($ctrl_sz + 4), font = :bold, color = :gray)
    btn_regen = Button(ctrl_grid[1, 4], label = "Regenerate Sources", width = 200)

    on(sl_mix.value) do val
        n_mix_obs[] = val
        update!()
    end
    on(btn_regen.clicks) do _
        sources_obs[] = [(-log.(rand(N_samples))) .* sign.(randn(N_samples)) for _ = 1:Max_Sources]
        update!() # Redraw the UI
    end

    # ── Layout: Data Panels ───────────────────────────────────────────────────

    # 1. Individual Independent Sources Stack (Left)
    ax_src = Axis(fig[3, 1], title = "The 50 Pure Sources", titlesize = title_sz)
    for i_val = 1:Max_Sources
        let i = i_val # Safe scope capture for the observable macro
            # Color them blue if they are included in the mix, gray if not
            line_col = @lift($n_mix_obs >= i ? :blue : :lightgray)
            # Offset them vertically heavily since Laplace tails reach +/- 7
            offset_src = @lift($(sources_obs)[i] .+ (i * 15.0))
            lines!(ax_src, 1:Plot_samples, @lift($(offset_src)[1:Plot_samples]), color = line_col, linewidth = 1.5)
        end
    end
    ylims!(ax_src, 0.0, (Max_Sources + 1) * 15.0)
    hidedecorations!(ax_src)


    # 2. Mixed Time Series (Middle)
    ax_ts = Axis(
        fig[3, 2],
        title = @lift("Sum of " * ($(n_mix_obs) == 1 ? "1 Source" : "$($n_mix_obs) Sources") * "\n(The Mixture)"),
        titlesize = title_sz,
    )
    lines!(ax_ts, 1:Plot_samples, @lift($(mix_data)[1:Plot_samples]), color = :gray, linewidth = 2)
    ylims!(ax_ts, -5, 5)
    hidedecorations!(ax_ts)


    # 3. Histogram (The Proof) (Right)
    ax_hist = Axis(
        fig[3, 3],
        title = @lift("Histogram of Mixture\n(Kurtosis: " * string(round($kurt_obs; digits = 2)) * ")"),
        titlesize = title_sz,
    )

    bins = range(-4, 4, length = 60)

    # Live histogram of the mixture
    hist!(ax_hist, mix_data, bins = bins, color = (:gray, 0.6), normalization = :pdf)

    # Overlay Perfect Gaussian Bell Curve to show the target
    x_gauss = range(-4, 4, length = 200)
    y_gauss = (1.0 / sqrt(2π)) .* exp.(-0.5 .* x_gauss .^ 2)
    lines!(ax_hist, x_gauss, y_gauss, color = :black, linestyle = :dash, linewidth = 4, label = "Perfect Gaussian")

    xlims!(ax_hist, -4, 4)
    ylims!(ax_hist, 0, 0.7)
    axislegend(ax_hist, position = :rt)

    # Note below the histogram
    Label(fig[4, 3], "Gaussian Kurtosis = 0.0\nSub-Gaussian (Flat) < 0.0\nSuper-Gaussian (Pointy) > 0.0", color = :gray)

    # Sizing
    colsize!(fig.layout, 1, Relative(0.20))
    colsize!(fig.layout, 2, Relative(0.40))
    colsize!(fig.layout, 3, Relative(0.40))

    rowsize!(fig.layout, 3, Relative(0.65))

    update!() # bootstrap
    display(fig)
    return fig
end
