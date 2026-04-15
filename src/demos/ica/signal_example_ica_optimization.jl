"""
    signal_example_ica_optimization()

Interactive Demo — Inside the ICA Black Box (Optimization).

Demystifies the ICA algorithm by showing how it searches for independent sources
as an optimization problem. Instead of looking at time-series outputs, this demo
visualizes the **search space** and the **learning rule**.

## Layout

| Panel | Content |
|-------|---------|
| Top-left | 2D scatter plot of sphered mixture. Red line is the projection vector `w`. |
| Top-right | Histogram of the projected data `u = wᵀX`. |
| Bottom | Optimization landscape (Kurtosis vs. Angle) with a dot showing current state. |

## Educational Value

Demonstrates that ICA is just "hill climbing" on an objective function.
By rotating the vector, the projected 1D histogram changes shape. The independent
sources correspond to angles where the histogram is most "peaky" (non-Gaussian).
The algorithm simply calculates the local gradient slope and steps uphill.

## See Also
- [`signal_example_ica_mixture`](@ref) — Part 1: What Is a Mixture?
- [`signal_example_ica_separation`](@ref) — Part 2: Mixing & Unmixing
- [`signal_example_ica_geometry`](@ref) — Part 3: 3 Sources & Scatter Geometry
- [`signal_example_ica_sphering`](@ref) — Part 4: Sphering (Whitening)

# Examples
```julia
using EegFun
    signal_example_ica_optimization()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica_optimization()
    # ── Constants ─────────────────────────────────────────────────────────────
    N = 4000
    Random.seed!(42) # For stable shapes across resets

    fig = Figure(size = (1200, 900), figure_padding = (15, 15, 15, 15))

    # Responsive sizing
    lbl_sz   = Observable(14)
    tick_sz  = Observable(12)
    title_sz = Observable(16)
    ctrl_sz  = Observable(14)

    on(fig.scene.viewport) do area
        sf         = area.widths[1] / 1200
        lbl_sz[]   = max(10, round(Int, 14 * sf))
        tick_sz[]  = max(8, round(Int, 12 * sf))
        title_sz[] = max(11, round(Int, 16 * sf))
        ctrl_sz[]  = max(10, round(Int, 14 * sf))
    end

    # ── Source Generation ─────────────────────────────────────────────────────
    # Make two highly super-Gaussian sources
    s1 = Float64[]
    s2 = Float64[]
    # Source 1: sharp bursts (like blinks)
    for _ = 1:(N÷400)
        append!(s1, zeros(300))
        append!(s1, [exp(-abs(i)/10.0) for i = -50:49])
    end
    # Source 2: oscillatory bursts (like alpha)
    for _ = 1:(N÷500)
        append!(s2, zeros(300))
        append!(s2, [1.5 * sin(2π * i / 20) * exp(-(i/40)^2) for i = -100:99])
    end

    # Pad to N if necessary
    s1 = vcat(s1, zeros(N - length(s1)))[1:N]
    s2 = vcat(s2, zeros(N - length(s2)))[1:N]

    # Add minor noise
    s1 .+= 0.05 .* randn(N)
    s2 .+= 0.05 .* randn(N)

    # ── Mixture and Sphering ──────────────────────────────────────────────────
    S = vcat(s1', s2')
    M = [1.0 0.8; 0.6 1.0]
    X_mixed = M * S

    # Center
    X_mixed .-= mean(X_mixed, dims = 2)

    # Sphere the data (PCA/Whitening) so the mixing is purely a rotation
    cov_X = (X_mixed * X_mixed') / N
    sphere_mat = inv(sqrt(cov_X))
    Z = sphere_mat * X_mixed   # Z is now sphered (identity covariance)

    # ── Objective Function Calculation ────────────────────────────────────────
    angles_deg = 0:1:180
    angles_rad = deg2rad.(angles_deg)

    kurtosis_curve = zeros(length(angles_rad))
    for (i, theta) in enumerate(angles_rad)
        w = [cos(theta), sin(theta)]
        u = w' * Z
        # Kurtosis ≈ E[u^4] / E[u^2]^2 - 3
        # Since data is sphered, E[u^2] ≈ 1
        kurt_val = mean((u .- mean(u)) .^ 4) / (mean((u .- mean(u)) .^ 2))^2 - 3.0
        kurtosis_curve[i] = kurt_val
    end

    # ── State Observables ─────────────────────────────────────────────────────
    angle_obs = Observable(15.0)   # Current angle in degrees

    proj_x_obs = Observable(zeros(10)) # Line segment x
    proj_y_obs = Observable(zeros(10)) # Line segment y

    hist_u_obs = Observable(zeros(N))  # 1D projected data
    kurt_text  = Observable("Kurtosis: 0.0")

    # ── Update Logic ──────────────────────────────────────────────────────────
    function update_state()
        theta = deg2rad(angle_obs[])
        w = [cos(theta), sin(theta)]

        # Update line geometry (from origin extending outwards through the blob)
        max_r = max(maximum(abs.(Z[1, :])), maximum(abs.(Z[2, :]))) * 1.2
        proj_x_obs[] = [-max_r * w[1], max_r * w[1]]
        proj_y_obs[] = [-max_r * w[2], max_r * w[2]]

        # Project data
        u = vec(w' * Z)
        hist_u_obs[] = u

        # Calc kurtosis
        kurt_val = mean((u .- mean(u)) .^ 4) / (mean((u .- mean(u)) .^ 2))^2 - 3.0
        kurt_text[] = "Kurtosis (Peakiness): $(round(kurt_val, digits=2))"
    end

    # ── Step ICA Learning Rule ────────────────────────────────────────────────
    # FastICA single step (gradient ascent on Negentropy approx G(u) = log(cosh(u)))
    # w+ = E[Z * g(w'Z)] - E[g'(w'Z)] * w
    # where g(u) = tanh(u)
    function execute_ica_step()
        theta = deg2rad(angle_obs[])
        w_old = [cos(theta), sin(theta)]

        u = vec(w_old' * Z)
        g_u = tanh.(u)
        g_prime_u = 1.0 .- g_u .^ 2

        # Gradient term
        term1 = vec(mean(Z .* g_u', dims = 2))
        term2 = mean(g_prime_u) .* w_old
        w_new = term1 .- term2

        # Normalize
        w_new ./= norm(w_new)

        # Convert back to angle
        new_theta_rad = atan(w_new[2], w_new[1])
        # Wrap to [0, 180) to match visual plot domain
        if new_theta_rad < 0
            new_theta_rad += π
        end
        new_theta_deg = rad2deg(new_theta_rad)

        angle_obs[] = new_theta_deg
        update_state()
    end

    # ── Build Layout ──────────────────────────────────────────────────────────
    Label(fig[0, 1:2], "Inside the Black Box: How ICA Finds Sources", fontsize = @lift($title_sz + 4), font = :bold, halign = :center)

    # Panel 1: Scatter and Projection Vector
    ax_scatter = Axis(
        fig[1, 1];
        title = "The Search Space\n(Data is 'Sphered' to be mathematically round, so ICA only needs to search rotations)",
        xlabel = "Component 1",
        ylabel = "Component 2",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
    )
    scatter!(ax_scatter, Z[1, :], Z[2, :]; color = (:gray50, 0.3), markersize = 3)
    lines!(ax_scatter, proj_x_obs, proj_y_obs; color = :red, linewidth = 4)
    # Hide ticks
    hidedecorations!(ax_scatter, label = false)
    ax_scatter.aspect = DataAspect()

    # Panel 2: Histogram
    ax_hist = Axis(
        fig[1, 2];
        title = "Projected 1D Data Distribution",
        xlabel = "Amplitude",
        ylabel = "Density",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
    )

    density!(ax_hist, hist_u_obs; color = (:blue, 0.5), strokecolor = :blue, strokewidth = 2)
    # Target standard normal reference
    x_range = -5:0.1:5
    lines!(
        ax_hist,
        x_range,
        x -> exp(-x^2/2) / sqrt(2π);
        color = (:red, 0.5),
        linestyle = :dash,
        linewidth = 2,
        label = "Gaussian (Target=0)",
    )
    axislegend(ax_hist; position = :lt, labelsize = lbl_sz)
    ylims!(ax_hist, 0.0, 1.2)
    xlims!(ax_hist, -4, 4)

    Label(fig[2, 2], kurt_text; fontsize = @lift($title_sz + 2), font = :bold, color = :blue, tellwidth = false)

    # Panel 3: The Landscape
    ax_land = Axis(
        fig[3, 1:2];
        title = "Objective Function Landscape (Kurtosis vs. Angle)",
        xlabel = "Vector Angle (Degrees)",
        ylabel = "Non-Gaussianity (Kurtosis)",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
    )

    lines!(ax_land, angles_deg, kurtosis_curve; color = :black, linewidth = 3)

    # Red dot current position
    dot_y = @lift(kurtosis_curve[clamp(round(Int, $angle_obs) + 1, 1, 181)])
    scatter!(ax_land, angle_obs, dot_y; color = :red, markersize = 20)

    xlims!(ax_land, 0, 180)
    ylims!(ax_land, 0, maximum(kurtosis_curve) * 1.2)

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[4, 1:2])

    sl_angle = Slider(ctrl[1, 1]; range = 0.0:0.5:180.0, startvalue = 15.0)
    Label(ctrl[1, 2], @lift("Angle: $(round($angle_obs, digits=1))°"); fontsize = ctrl_sz, width = 120)

    btn_step = Button(ctrl[1, 3]; label = "Step Algorithm", fontsize = ctrl_sz)
    btn_run  = Button(ctrl[1, 4]; label = "Run to Convergence", fontsize = ctrl_sz)

    # ── Wiring ────────────────────────────────────────────────────────────────
    on(sl_angle.value) do v
        angle_obs[] = v
        update_state()
    end

    # Sync manual slider if angle_obs changed programmatically (via Step)
    on(angle_obs) do v
        set_close_to!(sl_angle, v)
    end

    on(btn_step.clicks) do _
        execute_ica_step()
    end

    on(btn_run.clicks) do _
        @async begin
            for _ = 1:20
                execute_ica_step()
                sleep(0.05)
            end
        end
    end

    # ── Final Layout Tweaks ───────────────────────────────────────────────────
    rowsize!(fig.layout, 0, Relative(0.05))
    rowsize!(fig.layout, 1, Relative(0.40))
    rowsize!(fig.layout, 2, Relative(0.05))
    rowsize!(fig.layout, 3, Relative(0.35))
    rowsize!(fig.layout, 4, Relative(0.10))

    update_state()
    display(fig)
    return fig
end
