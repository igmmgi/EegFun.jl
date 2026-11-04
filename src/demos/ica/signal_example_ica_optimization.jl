"""
    signal_example_ica_optimization()

Interactive ICA optimization demo visualizing the sphered 2D search space, projected 1D density,
kurtosis objective landscape, and single-step / convergent gradient ascent.

## See Also
- [`signal_example_ica_mixture`](@ref), [`signal_example_ica_separation`](@ref),
  [`signal_example_ica_geometry`](@ref), [`signal_example_ica_sphering`](@ref)

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_optimization()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica_optimization()
    _check_interactive_demo_backend("ICA Optimization demo") || return nothing

    N = 4000
    rng = Random.Xoshiro(42)

    fig = _eegfun_figure(size = (1200, 900), figure_padding = (15, 15, 15, 15))
    lbl_sz, _, title_sz, ctrl_sz = _demo_adaptive_fonts(fig, 1200, (14, 10), (12, 8), (16, 11), (14, 10))

    s1 = zeros(Float64, N)
    spike = [exp(-abs(i) / 10.0) for i = -50:49]
    for blk = 1:(N÷400)
        offset = (blk - 1) * 400 + 300
        @views s1[(offset+1):(offset+100)] .= spike
    end

    s2 = zeros(Float64, N)
    gabor = [1.5 * sin(2π * i / 20) * exp(-(i / 40)^2) for i = -100:99]
    for blk = 1:(N÷500)
        offset = (blk - 1) * 500 + 300
        @views s2[(offset+1):(offset+200)] .= gabor
    end

    s1 .+= 0.05 .* randn(rng, N)
    s2 .+= 0.05 .* randn(rng, N)

    S = vcat(s1', s2')
    X_mixed = [1.0 0.8; 0.6 1.0] * S
    X_mixed .-= mean(X_mixed, dims = 2)

    Z = inv(sqrt((X_mixed * X_mixed') / N)) * X_mixed
    z1 = Z[1, :]
    z2 = Z[2, :]
    max_r = max(maximum(abs, z1), maximum(abs, z2)) * 1.2

    function _projected_kurtosis(w1::Float64, w2::Float64)
        m1 = m2 = m4 = 0.0
        @inbounds for i = 1:N
            u = w1 * z1[i] + w2 * z2[i]
            m1 += u
        end
        μ = m1 / N
        @inbounds for i = 1:N
            d = (w1 * z1[i] + w2 * z2[i]) - μ
            d2 = d * d
            m2 += d2
            m4 += d2 * d2
        end
        v = m2 / N
        return (m4 / N) / (v * v) - 3.0
    end

    angles_deg = 0:1:180
    kurtosis_curve = Vector{Float64}(undef, length(angles_deg))
    for (i, deg) in enumerate(angles_deg)
        sθ, cθ = sincos(deg2rad(deg))
        kurtosis_curve[i] = _projected_kurtosis(cθ, sθ)
    end

    angle_obs = Observable(15.0)
    proj_x_obs = Observable(zeros(2))
    proj_y_obs = Observable(zeros(2))
    hist_u_obs = Observable(zeros(N))
    kurt_text = Observable("Kurtosis: 0.0")

    function update_state()
        sθ, cθ = sincos(deg2rad(angle_obs[]))
        proj_x_obs[] = [-max_r * cθ, max_r * cθ]
        proj_y_obs[] = [-max_r * sθ, max_r * sθ]

        u = Vector{Float64}(undef, N)
        @inbounds for i = 1:N
            u[i] = cθ * z1[i] + sθ * z2[i]
        end
        hist_u_obs[] = u
        kurt_text[] = "Kurtosis (Peakiness): $(round(_projected_kurtosis(cθ, sθ), digits=2))"
    end

    function execute_ica_step()
        sθ, cθ = sincos(deg2rad(angle_obs[]))
        t1_x = 0.0
        t1_y = 0.0
        gp_sum = 0.0
        @inbounds for i = 1:N
            gu = tanh(cθ * z1[i] + sθ * z2[i])
            t1_x += z1[i] * gu
            t1_y += z2[i] * gu
            gp_sum += 1.0 - gu * gu
        end
        inv_n = 1.0 / N
        gp_mean = gp_sum * inv_n
        wx = t1_x * inv_n - gp_mean * cθ
        wy = t1_y * inv_n - gp_mean * sθ
        inv_norm = 1.0 / hypot(wx, wy)
        wx *= inv_norm
        wy *= inv_norm

        new_theta_rad = atan(wy, wx)
        new_theta_rad < 0 && (new_theta_rad += π)
        angle_obs[] = rad2deg(new_theta_rad)
        update_state()
    end

    Label(fig[0, 1:2], "Inside the Black Box: How ICA Finds Sources", fontsize = @lift($title_sz + 4), font = :bold, halign = :center)

    ax_scatter = Axis(
        fig[1, 1];
        title = "The Search Space\n(Data is 'Sphered' to be mathematically round, so ICA only needs to search rotations)",
        xlabel = "Component 1",
        ylabel = "Component 2",
        titlesize = title_sz,
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        aspect = DataAspect(),
    )
    scatter!(ax_scatter, z1, z2; color = (:gray50, 0.3), markersize = 3)
    lines!(ax_scatter, proj_x_obs, proj_y_obs; color = :red, linewidth = 4)
    hidedecorations!(ax_scatter, label = false)

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
    x_range = -5:0.1:5
    lines!(
        ax_hist,
        x_range,
        [exp(-x^2 / 2) / sqrt(2π) for x in x_range];
        color = (:red, 0.5),
        linestyle = :dash,
        linewidth = 2,
        label = "Gaussian (Target=0)",
    )
    axislegend(ax_hist; position = :lt, labelsize = lbl_sz)
    ylims!(ax_hist, 0.0, 1.2)
    xlims!(ax_hist, -4, 4)

    Label(fig[2, 2], kurt_text; fontsize = @lift($title_sz + 2), font = :bold, color = :blue, tellwidth = false)

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
    dot_y = @lift(kurtosis_curve[clamp(round(Int, $angle_obs) + 1, 1, 181)])
    scatter!(ax_land, angle_obs, dot_y; color = :red, markersize = 20)
    xlims!(ax_land, 0, 180)
    ylims!(ax_land, 0, maximum(kurtosis_curve) * 1.2)

    ctrl = GridLayout(fig[4, 1:2])
    sl_angle = Slider(ctrl[1, 1]; range = 0.0:0.5:180.0, startvalue = 15.0)
    Label(ctrl[1, 2], @lift("Angle: $(round($angle_obs, digits=1))°"); fontsize = ctrl_sz, width = 120)
    btn_step = Button(ctrl[1, 3]; label = "Step Algorithm", fontsize = ctrl_sz)
    btn_run = Button(ctrl[1, 4]; label = "Run to Convergence", fontsize = ctrl_sz)

    is_syncing_angle = Ref(false)
    on(sl_angle.value) do v
        is_syncing_angle[] && return
        is_syncing_angle[] = true
        try
            angle_obs[] = v
            update_state()
        finally
            is_syncing_angle[] = false
        end
    end
    on(angle_obs) do v
        is_syncing_angle[] && return
        is_syncing_angle[] = true
        try
            set_close_to!(sl_angle, v)
        finally
            is_syncing_angle[] = false
        end
    end

    on(_ -> execute_ica_step(), btn_step.clicks)

    is_running = Ref(false)
    on(btn_run.clicks) do _
        is_running[] && return
        is_running[] = true
        @async begin
            try
                for _ = 1:20
                    (!is_running[] || !events(fig).window_open[]) && break
                    execute_ica_step()
                    sleep(0.05)
                end
            finally
                is_running[] = false
            end
        end
    end

    rowsize!(fig.layout, 0, Relative(0.05))
    rowsize!(fig.layout, 1, Relative(0.40))
    rowsize!(fig.layout, 2, Relative(0.05))
    rowsize!(fig.layout, 3, Relative(0.35))
    rowsize!(fig.layout, 4, Relative(0.10))

    update_state()
    _display_demo_figure(fig)
    return fig
end
