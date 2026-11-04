"""
    signal_example_ica_sphering()

Interactive ICA preprocessing demo showing how centering, PCA rotation, and variance scaling
(whitening/sphering) transform a correlated mixture into a rotationally symmetric search space.

## See Also
- [`signal_example_ica_geometry`](@ref), [`signal_example_ica_optimization`](@ref)

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_sphering()
```
"""
function signal_example_ica_sphering()
    _check_interactive_demo_backend("ICA Sphering demo") || return nothing

    N = 2500
    rng = Random.Xoshiro(101)

    fig = _eegfun_figure(size = (1000, 700), figure_padding = (15, 15, 15, 15))
    lbl_sz = Observable(14)
    title_sz = Observable(16)

    burst = [sin(2π * i / 20) * exp(-(i / 40)^2) for i = -100:99]
    s1 = zeros(Float64, N)
    s1[401:600] .= burst
    s1[1001:1200] .= burst

    s2 = randn(rng, N) .* 0.3
    spike = [exp(-abs(i) / 15.0) for i = -50:49]
    for _ = 1:6
        idx = rand(rng, 200:(N-200))
        @views s2[idx:(idx+99)] .+= spike
    end

    S = Matrix{Float64}(undef, 2, N)
    @inbounds for j = 1:N
        S[1, j] = s1[j]
        S[2, j] = s2[j]
    end

    M = [1.2 0.8; 0.5 1.4]
    X_raw = M * S .+ [1.5, 0.8]
    X_centered = X_raw .- mean(X_raw, dims = 2)

    F = eigen((X_centered * X_centered') / N)
    order = sortperm(F.values, rev = true)
    eigvals = F.values[order]
    V = F.vectors[:, order]
    det(V) < 0 && (V[:, 2] .*= -1)

    X_rotated = V' * X_centered
    X_sphered = diagm(1.0 ./ sqrt.(eigvals)) * X_rotated

    x_buf = copy(X_raw[1, :])
    y_buf = copy(X_raw[2, :])
    x_obs = Observable(x_buf)
    y_obs = Observable(y_buf)
    title_text = Observable("Step 1: The Raw EEG Mixture")
    desc_text = Observable("The data arrives from the electrodes correlated (diagonal), stretched, and off-center.")

    function _lerp_rows!(A::Matrix{Float64}, B::Matrix{Float64}, α::Float64)
        wA = 1.0 - α
        @inbounds for j = 1:N
            x_buf[j] = muladd(wA, A[1, j], α * B[1, j])
            y_buf[j] = muladd(wA, A[2, j], α * B[2, j])
        end
        notify(x_obs)
        notify(y_obs)
    end

    function update_morph(t)
        if t <= 1.0
            _lerp_rows!(X_raw, X_centered, t)
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
            _lerp_rows!(X_centered, X_rotated, t - 1.0)
            if t == 2.0
                title_text[] = "Step 3: PCA Rotated (Decorrelated)"
                desc_text[] = "We rotate the data to align its principal axes with X and Y so the two sensors are uncorrelated."
            else
                title_text[] = "Morphing: PCA Rotation..."
                desc_text[] = "Applying Principal Component Analysis to decorrelate the sensors."
            end
        else
            _lerp_rows!(X_rotated, X_sphered, t - 2.0)
            if t == 3.0
                title_text[] = "Step 4: Sphered (Whitened)"
                desc_text[] = "We scale each axis by its standard deviation until the blob is spherical — now ICA only searches angles!"
            else
                title_text[] = "Morphing: Variance Scaling..."
                desc_text[] = "Dividing each axis by its standard deviation to equalize variance."
            end
        end
    end

    Label(fig[0, 1], "Preprocessing: The 'Sphering' Mathematical Shortcut", fontsize = @lift($title_sz + 4), font = :bold)
    Label(fig[1, 1], desc_text, fontsize = title_sz, tellwidth = false)

    ax = Axis(
        fig[2, 1],
        title = title_text,
        titlesize = @lift($title_sz + 2),
        xlabel = "Virtual Sensor 1",
        ylabel = "Virtual Sensor 2",
        xlabelsize = lbl_sz,
        ylabelsize = lbl_sz,
        aspect = DataAspect(),
    )
    scatter!(ax, x_obs, y_obs; color = (:purple, 0.4), markersize = 4)
    vlines!(ax, 0.0, color = (:black, 0.2), linewidth = 2, linestyle = :dash)
    hlines!(ax, 0.0, color = (:black, 0.2), linewidth = 2, linestyle = :dash)
    xlims!(ax, -4, 4)
    ylims!(ax, -4, 4)

    ctrl = GridLayout(fig[3, 1])
    Label(ctrl[1, 1], "Drag to Transformation Step:", font = :bold)
    sl_morph = Slider(ctrl[1, 2], range = 0.0:0.05:3.0, startvalue = 0.0)

    btn_box = GridLayout(ctrl[2, 1:2])
    for (col, lbl, target) in ((1, "1. Raw", 0.0), (2, "2. Center", 1.0), (3, "3. PCA", 2.0), (4, "4. Sphere", 3.0))
        btn = Button(btn_box[1, col], label = lbl)
        let tgt = target
            on(_ -> set_close_to!(sl_morph, tgt), btn.clicks)
        end
    end

    on(update_morph, sl_morph.value)

    rowsize!(fig.layout, 0, Relative(0.05))
    rowsize!(fig.layout, 1, Relative(0.05))
    rowsize!(fig.layout, 2, Relative(0.70))
    rowsize!(fig.layout, 3, Relative(0.15))

    update_morph(0.0)
    _display_demo_figure(fig)
    return fig
end
