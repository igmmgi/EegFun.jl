"""
    signal_example_ica_sphering()

Interactive Demo — Preprocessing: Sphering and Whitening.

Bridges the gap between the chaotic mixture problem and the clean rotational search 
space of ICA optimization. Demonstrates how preprocessing the data drastically
simplifies the algorithm's job.

## Educational Value
Watching a messy, correlated, off-center EEG mixture mathematically morph into 
a perfectly round, uncorrelated circle proves *why* ICA can use simple rotational
algorithms to find the sources.

## See Also
- [`signal_example_ica_geometry`](@ref) — Part 3: 3 Sources & Scatter Geometry
- [`signal_example_ica_optimization`](@ref) — Part 5: Inside the Black Box (Optimization)

# Examples
```julia
using EegFun
    signal_example_ica_sphering()
```
"""
function signal_example_ica_sphering()
    N = 2500
    Random.seed!(101) # Stable random numbers

    fig = Figure(size = (1000, 700), figure_padding = (15, 15, 15, 15))

    lbl_sz   = Observable(14)
    tick_sz  = Observable(12)
    title_sz = Observable(16)
    ctrl_sz  = Observable(14)

    # ── Source Generation ─────────────────────────────────────────────────────
    s1 = [sin(2π * i / 20) * exp(-(i/40)^2) for i = -100:99]
    s1 = vcat([zeros(400); s1; zeros(400); s1]...)
    s1 = vcat(s1, zeros(max(0, N - length(s1))))[1:N]

    s2 = randn(N) .* 0.3
    for _ = 1:6
        idx = rand(200:(N-200))
        s2[idx:(idx+99)] .+= [exp(-abs(i)/15.0) for i = -50:49]
    end

    S = vcat(s1', s2')

    # Mix them and shift them specifically so Centering has something to do
    M = [1.2 0.8; 0.5 1.4]
    X_raw = M * S .+ [1.5, 0.8]

    # Calculate transformation anchors
    X_mean = mean(X_raw, dims = 2)
    X_centered = X_raw .- X_mean

    # PCA 
    C = (X_centered * X_centered') / N
    F = eigen(C)
    # Sort eigenvalues descending to ensure consistent rotation
    order = sortperm(F.values, rev = true)
    eigvals = F.values[order]
    V = F.vectors[:, order]

    # Ensure determinant is 1 (rotation not reflection)
    if det(V) < 0
        V[:, 2] .*= -1
    end

    X_rotated = V' * X_centered

    # Sphering
    D_inv_sqrt = diagm(1.0 ./ sqrt.(eigvals))
    X_sphered = D_inv_sqrt * X_rotated

    # ── Observables ───────────────────────────────────────────────────────────
    step_obs = Observable(1) # 1=Raw, 2=Centered, 3=Rotated, 4=Sphered

    # Actually continuous slider!
    morph_val = Observable(0.0) # 0.0 to 3.0

    current_X = Observable(copy(X_raw))
    title_text = Observable("Step 1: The Raw EEG Mixture")
    desc_text = Observable("The data arrives from the electrodes correlated (diagonal), stretched, and off-center.")

    # ── Morph Logic ───────────────────────────────────────────────────────────
    function update_morph(t)
        if t <= 1.0
            # 0.0 to 1.0 -> Morph Raw to Centered
            current_X[] = (1.0 - t) .* X_raw .+ t .* X_centered
            if t == 0.0
                title_text[] = "Step 1: The Raw EEG Mixture"
                desc_text[] = "The data arrives from the electrodes correlated (diagonal), stretched, and off-center."
            elseif t == 1.0
                title_text[] = "Step 2: Centered Data"
                desc_text[] = "We subtract the mean so the data blob is anchored perfectly at (0,0)."
            else
                title_text[] = "Morphing: Centering..."
                desc_text[] = "Subtracting the mean from all channels."
            end
        elseif t <= 2.0
            # 1.0 to 2.0 -> Morph Centered to Rotated
            t2 = t - 1.0
            current_X[] = (1.0 - t2) .* X_centered .+ t2 .* X_rotated
            if t == 2.0
                title_text[] = "Step 3: PCA Rotated (Decorrelated)"
                desc_text[] = "We rotate the data to perfectly align its longest arms with the X and Y axes. Now the two sensors are uncorrelated."
            else
                title_text[] = "Morphing: PCA Rotation..."
                desc_text[] = "Applying Principal Component Analysis to decorrelate the sensors."
            end
        else
            # 2.0 to 3.0 -> Morph Rotated to Sphered
            t3 = t - 2.0
            current_X[] = (1.0 - t3) .* X_rotated .+ t3 .* X_sphered
            if t == 3.0
                title_text[] = "Step 4: Sphered (Whitened)"
                desc_text[] = "We squish the axes by their variance until the blob is a perfect mathematical sphere. Now ICA only has to search for angles!"
            else
                title_text[] = "Morphing: Variance Scaling..."
                desc_text[] = "Dividing each axis by its standard deviation to equalize variance."
            end
        end
    end

    # ── Layout ────────────────────────────────────────────────────────────────
    Label(fig[0, 1], "Preprocessing: The 'Sphering' Mathematical Shortcut", fontsize = @lift($title_sz + 4), font = :bold)

    Label(fig[1, 1], desc_text, fontsize = title_sz, tellwidth = false)

    ax = Axis(
        fig[2, 1],
        title = title_text,
        titlesize = @lift($title_sz+2),
        xlabel = "Virtual Sensor 1",
        ylabel = "Virtual Sensor 2",
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
    )

    scatter!(ax, @lift($current_X[1, :]), @lift($current_X[2, :]); color = (:purple, 0.4), markersize = 4)

    # Add origin reference lines
    vlines!(ax, 0.0, color = (:black, 0.2), linewidth = 2, linestyle = :dash)
    hlines!(ax, 0.0, color = (:black, 0.2), linewidth = 2, linestyle = :dash)

    # Keep limits fixed so the morphing is obvious
    xlims!(ax, -4, 4)
    ylims!(ax, -4, 4)
    ax.aspect = DataAspect()

    # ── Controls ──────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[3, 1])

    Label(ctrl[1, 1], "Drag to Transformation Step:", font = :bold)
    sl_morph = Slider(ctrl[1, 2], range = 0.0:0.05:3.0, startvalue = 0.0)

    btn_box = GridLayout(ctrl[2, 1:2])
    btn1 = Button(btn_box[1, 1], label = "1. Raw")
    btn2 = Button(btn_box[1, 2], label = "2. Center")
    btn3 = Button(btn_box[1, 3], label = "3. PCA")
    btn4 = Button(btn_box[1, 4], label = "4. Sphere")

    # ── Wiring ────────────────────────────────────────────────────────────────
    on(sl_morph.value) do val
        update_morph(val)
    end

    on(btn1.clicks) do _
        ;
        set_close_to!(sl_morph, 0.0);
    end
    on(btn2.clicks) do _
        ;
        set_close_to!(sl_morph, 1.0);
    end
    on(btn3.clicks) do _
        ;
        set_close_to!(sl_morph, 2.0);
    end
    on(btn4.clicks) do _
        ;
        set_close_to!(sl_morph, 3.0);
    end

    rowsize!(fig.layout, 0, Relative(0.05))
    rowsize!(fig.layout, 1, Relative(0.05))
    rowsize!(fig.layout, 2, Relative(0.70))
    rowsize!(fig.layout, 3, Relative(0.15))

    update_morph(0.0)
    display(fig)
    return fig
end
