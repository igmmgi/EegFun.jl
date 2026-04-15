"""
    signal_example_ica_infomax()

Interactive Demo — Information Maximization (Infomax) ICA.

Demonstrates the exact mathematical mechanism defined by Information Theory.
Independent Component Analysis assumes that true sources are statistically 
independent. In Information Theory, two variables are independent if their 
Mutual Information is exactly zero.

Because the data is sphered (Joint Entropy is constant), minimizing the 
Mutual Information between components is mathematically equivalent to minimizing 
their individual (Marginal) Entropies. Since a Gaussian distribution has the 
maximum possible entropy, finding the lowest-entropy sources is the same as 
finding the least-Gaussian (highest kurtosis) sources!

## Layout
- Top Left: Histogram of Component 1 (with live Entropy).
- Top Center: 2D Scatter Plot of Component 1 vs Component 2.
- Top Right: Histogram of Component 2 (with live Entropy).
- Bottom: The Information Theory landscape (Mutual Information vs. Angle).

# Examples
```julia
using EegFun
using StatsBase

EegFun.signal_example_ica_infomax()
```
"""
function signal_example_ica_infomax()
    N = 50000
    Random.seed!(42)

    fig = Figure(size = (1600, 900), figure_padding = (15, 15, 15, 15))

    lbl_sz   = Observable(14)
    tick_sz  = Observable(12)
    title_sz = Observable(16)
    ctrl_sz  = Observable(16)

    # ── Source Generation ─────────────────────────────────────────────────────
    # Generate Highly Sparse Super-Gaussian Sources (Laplacian)
    s1 = (-log.(rand(N))) .* sign.(randn(N))
    s2 = (-log.(rand(N))) .* sign.(randn(N))

    s1 = (s1 .- mean(s1)) ./ std(s1)
    s2 = (s2 .- mean(s2)) ./ std(s2)

    # ── Mix and Sphere ────────────────────────────────────────────────────────
    S = vcat(s1', s2')
    M = [1.0 0.8; 0.6 1.0]
    X_mixed = M * S

    X_mixed .-= mean(X_mixed, dims=2)
    cov_X = (X_mixed * X_mixed') / N
    sphere_mat = inv(sqrt(cov_X))
    Z = sphere_mat * X_mixed   # Z is sphered

    # ── The Infomax Math ──────────────────────────────────────────────────────
    function calc_mutual_info(y1, y2)
        # We manually bin the space to calculate continuous Mutual Information
        edges = range(-5, 5, length=50)
        h2d = fit(Histogram, (y1, y2), (edges, edges)).weights
        p_xy = h2d ./ sum(h2d)
        
        p_x = sum(p_xy, dims=2)
        p_y = sum(p_xy, dims=1)
        
        H_xy = sum(-p_xy[p_xy .> 0] .* log2.(p_xy[p_xy .> 0]))
        H_x  = sum(-p_x[p_x .> 0] .* log2.(p_x[p_x .> 0]))
        H_y  = sum(-p_y[p_y .> 0] .* log2.(p_y[p_y .> 0]))
        
        # Mutual Info = H(X) + H(Y) - H(X,Y)
        return max(0.0, sum(H_x) + sum(H_y) - sum(H_xy)), sum(H_x), sum(H_y)
    end

    # Pre-compute Landscape
    angles_deg = 0:1:180
    angles_rad = deg2rad.(angles_deg)
    
    mi_curve = zeros(length(angles_rad))
    for (i, theta) in enumerate(angles_rad)
        # Rotation matrix
        W = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        Y = W * Z
        mi_val, _, _ = calc_mutual_info(vec(Y[1, :]), vec(Y[2, :]))
        mi_curve[i] = mi_val
    end

    # ── Observables ───────────────────────────────────────────────────────────
    angle_obs = Observable(30.0) 
    
    y1_obs = Observable(zeros(N))
    y2_obs = Observable(zeros(N))
    
    y1_subset = Observable(zeros(500)) # scatter performance
    y2_subset = Observable(zeros(500))
    
    mi_text = Observable("Mutual Information: 0.00 bits")
    hy1_text = Observable("H(Y1): 0.00 bits")
    hy2_text = Observable("H(Y2): 0.00 bits")

    # ── Update Logic ──────────────────────────────────────────────────────────
    function update_state()
        theta = deg2rad(angle_obs[])
        W = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        Y = W * Z
        
        y1 = vec(Y[1, :])
        y2 = vec(Y[2, :])
        
        y1_obs[] = y1
        y2_obs[] = y2
        
        y1_subset[] = y1[1:800]
        y2_subset[] = y2[1:800]
        
        mi_val, h1_val, h2_val = calc_mutual_info(y1, y2)
        
        mi_text[] = "Mutual Information: $(round(mi_val, digits=3)) bits"
        hy1_text[] = "Entropy H(y₁): $(round(h1_val, digits=2))"
        hy2_text[] = "Entropy H(y₂): $(round(h2_val, digits=2))"
    end

    # ── Build Layout ──────────────────────────────────────────────────────────
    Label(fig[0, 1:3], "Infomax ICA: Information Maximization (Bell & Sejnowski 1995)",
          fontsize = @lift($title_sz + 6), font = :bold)

    # Panel 1: Component 1 Histogram
    ax_y1 = Axis(fig[1, 1]; title = "Component 1 (Output)", 
                xlabel = "Amplitude", ylabel = "Probability Density",
                titlesize = title_sz, xlabelsize = lbl_sz, ylabelsize = lbl_sz)
    
    bins_u = range(-5, 5, length=80)
    hist!(ax_y1, y1_obs, bins=bins_u, color=(:blue, 0.5), normalization=:pdf)
    xlims!(ax_y1, -5, 5)
    ylims!(ax_y1, 0, 0.8)
    
    # Panel 2: The Scatter Plot (Joint Distribution)
    ax_scatter = Axis(fig[1, 2]; title = "Joint Distribution (Unmixed Space)", 
                  xlabel = "Component 1", ylabel = "Component 2",
                  titlesize = title_sz, xlabelsize = lbl_sz, ylabelsize = lbl_sz)
    
    scatter!(ax_scatter, y1_subset, y2_subset; color = (:black, 0.4), markersize = 5)
    xlims!(ax_scatter, -5, 5)
    ylims!(ax_scatter, -5, 5)
    ax_scatter.aspect = DataAspect()

    # Panel 3: Component 2 Histogram
    ax_y2 = Axis(fig[1, 3]; title = "Component 2 (Output)", 
                xlabel = "Amplitude", ylabel = "Probability Density",
                titlesize = title_sz, xlabelsize = lbl_sz, ylabelsize = lbl_sz)
    
    hist!(ax_y2, y2_obs, bins=bins_u, color=(:red, 0.5), normalization=:pdf)
    xlims!(ax_y2, -5, 5)
    ylims!(ax_y2, 0, 0.8)

    # Entropies below their respective blocks
    Label(fig[2, 1], hy1_text; fontsize = @lift($title_sz + 2), font = :bold, color = :blue)
    Label(fig[2, 2], mi_text; fontsize = @lift($title_sz + 6), font = :bold, color = :purple)
    Label(fig[2, 3], hy2_text; fontsize = @lift($title_sz + 2), font = :bold, color = :red)

    # ── Panel 4: The Landscape ────────────────────────────────────────────────
    ax_land = Axis(fig[3, 1:3]; title = "Objective Function: Mutual Information I(y₁, y₂)", 
                   subtitle = "Goal: Rotate unmixing matrix until signals share Zero Mutual Information.",
                   xlabel = "Search Angle (Degrees)", ylabel = "Mutual Information (Bits)",
                   titlesize = title_sz, xlabelsize = lbl_sz, ylabelsize = lbl_sz)
                   
    lines!(ax_land, angles_deg, mi_curve; color = :black, linewidth = 3)
    
    # Red dot current position
    dot_y = @lift(mi_curve[clamp(round(Int, $angle_obs) + 1, 1, 181)])
    scatter!(ax_land, angle_obs, dot_y; color = :purple, markersize = 20)
    
    xlims!(ax_land, 0, 180)
    
    min_ent = minimum(mi_curve)
    max_ent = maximum(mi_curve)
    ylims!(ax_land, max(0.0, min_ent - 0.05), max_ent + 0.05)

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[4, 1:3])
    
    Label(ctrl[1, 1], "Rotate Angle:"; fontsize = ctrl_sz, font=:bold, halign=:right)
    sl_angle = Slider(ctrl[1, 2]; range = 0.0:0.5:180.0, startvalue = 30.0)
    Label(ctrl[1, 3], @lift(" $(round($angle_obs, digits=1))°"); fontsize = ctrl_sz, halign=:left)
    
    # ── Wiring ────────────────────────────────────────────────────────────────
    on(sl_angle.value) do v
        angle_obs[] = v
        update_state()
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
