"""
    signal_example_ica()

Interactive ICA Demo — Blind Source Separation.

Demonstrates how Independent Component Analysis (ICA) recovers independent
source signals from their linear mixtures. Uses a 3-source / 3-sensor
scenario where mixing is a composition of three 2D rotations (one per
source pair), making the geometry tractable while realistic.

## Layout

| Column | Row 1             | Row 2                  | Row 3                   |
|--------|-------------------|------------------------|-------------------------|
| Left   | True sources      | EEG recordings (mixed) | Recovered components    |
| Middle | Sources S1 vs S2  | Mixed M1 vs M2         | Recovered C1 vs C2      |
| Right  | Sources S1 vs S3  | Mixed M1 vs M3         | Recovered C1 vs C3      |

The scatter plots show pairwise joint distributions. Independent sources show
a distinctive **non-blob** shape. Mixing rotates them into correlated blobs.
Correct unmixing restores the shape.

## Sources

- **S1**: Neural oscillation — sine wave at user-controlled frequency
- **S2**: Blink artifact — sparse gamma-shaped pulses (high kurtosis)
- **S3**: Muscle burst — Gabor pulse (sinusoid × Gaussian envelope, high kurtosis)

## Controls

| Control         | Range    | Description                                                  |
|:----------------|:---------|:-------------------------------------------------------------|
| Mix α (S1-S2)   | 0 – 90°  | Rotation blending S1 and S2 into the recordings              |
| Mix β (S1-S3)   | 0 – 90°  | Rotation blending S1 and S3 into the recordings              |
| Mix γ (S2-S3)   | 0 – 90°  | Rotation blending S2 and S3 into the recordings              |
| Unmix φ (S1-S2) | 0 – 90°  | Inverse rotation used to recover S1/S2 (manual mode)        |
| Unmix ψ (S1-S3) | 0 – 90°  | Inverse rotation used to recover S1/S3 (manual mode)        |
| Unmix χ (S2-S3) | 0 – 90°  | Inverse rotation used to recover S2/S3 (manual mode)        |
| Signal Freq     | 1 – 20 Hz| Frequency of the oscillatory source (S1)                     |
| Noise           | 0 – 1    | Additive Gaussian noise on all sources                       |
| Auto-ICA button | —        | Runs Infomax ICA on the mixed signals to find W              |

## Teaching Notes

**The problem**: EEG electrodes record *mixtures* of all underlying sources.
Each electrode signal is a weighted sum — you never observe sources directly.

**The geometry**: With N sources and N sensors, mixing is an N×N rotation.
For 3 sources we compose three pairwise 2D rotations (α, β, γ).
Unmixing reverses the rotations in the opposite order with the same angles.

**ICA's trick**: It finds an unmixing matrix W such that components W·M are
maximally *non-Gaussian*. Independent signals have higher combined non-Gaussianity
than their mixtures (Central Limit Theorem in reverse).

**The scatter plot**: Scatter of two independent non-Gaussian signals shows a
distinctive cross / T-shape. Mixed signals produce a correlated blob. The Auto-ICA
button runs the real Infomax algorithm (Bell & Sejnowski, 1995) — the same one
used by `run_ica` in EegFun — to find the unmixing matrix W directly.

**What to try**:
1. Set all mixing angles > 0 and watch the scatters collapse into blobs.
2. Drag unmixing sliders manually to match the mixing angles — sources reappear.
3. Reset unmixing angles to 0, then click Auto-ICA — Infomax recovers the sources.
4. Add Noise — Infomax still works because S2 and S3 are highly non-Gaussian.
5. Compare manual rotation (angle sliders) with Infomax (Auto-ICA) — Infomax
   finds a general unmixing matrix, not just a rotation.

## See Also

- Bell, A. J., & Sejnowski, T. J. (1995). An information-maximization approach
  to blind separation. *Neural Computation*, *7*(6), 1129–1159.
- Hyvärinen, A., & Oja, E. (2000). Independent component analysis: algorithms and
  applications. *Neural Networks*, *13*(4-5), 411–430.

# Examples
```julia
using EegFun
EegFun.signal_example_ica()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica()

    FS = 512.0
    T  = 3.0

    fig = Figure(size = (2100, 900), figure_padding = (8, 8, 8, 8))

    # ── Responsive font sizing ─────────────────────────────────────────────────
    lbl_sz   = Observable(13)
    tick_sz  = Observable(11)
    title_sz = Observable(14)
    ctrl_sz  = Observable(13)

    on(fig.scene.viewport) do area
        sf         = area.widths[1] / 3600
        lbl_sz[]   = max(10, round(Int, 13 * sf))
        tick_sz[]  = max(8, round(Int, 11 * sf))
        title_sz[] = max(11, round(Int, 14 * sf))
        ctrl_sz[]  = max(10, round(Int, 13 * sf))
    end

    t_arr = collect(range(0.0, T - 1.0 / FS, step = 1.0 / FS))
    n     = length(t_arr)

    # ── Observables ────────────────────────────────────────────────────────────
    mix_α    = Observable(25.0)   # S1-S2 plane
    mix_β    = Observable(15.0)   # S1-S3 plane
    mix_γ    = Observable(20.0)   # S2-S3 plane
    unmix_α  = Observable(25.0)
    unmix_β  = Observable(15.0)
    unmix_γ  = Observable(20.0)
    sig_freq  = Observable(8.0)
    noise_lvl = Observable(0.1)

    # Infomax state: when true rec signals come from the ICA result,
    # not the manual rotation sliders. Moving any unmix slider reverts to false.
    use_infomax   = Observable(false)
    ica_rec       = Observable(zeros(3, n))   # 3 × n recovered matrix
    in_ica_update = Ref(false)                # true while ICA moves sliders programmatically

    # ── Source generation ──────────────────────────────────────────────────────
    #
    # All three sources are sparse / episodic so that pairwise scatter plots
    # show clean cross shapes (they are rarely large at the same time):
    #
    # S1: oscillatory bursts (e.g. alpha events) at t = 0.25, 1.05, 2.55 s
    # S2: blink artifact gamma pulses        at t = 0.65, 1.50, 2.10 s
    # S3: muscle burst (Gabor)               at t = 1.80 s
    #
    # sig_freq controls the oscillation frequency inside the S1 Gabor envelopes.
    #
    """Generate three synthetic source signals: oscillatory bursts, blink artifacts, and muscle burst."""
    function make_sources(freq, noise)
        # S1: short oscillatory bursts (Gabor: sinusoid × Gaussian envelope)
        s1 = zeros(n)
        for center in [0.25, 1.05, 2.55]
            s1 .+= 2.0 .* sin.(2π .* freq .* (t_arr .- center)) .* exp.(.-((t_arr .- center) .^ 2) ./ (2 .* 0.09^2))
        end

        # S2: blink artifact
        s2 = zeros(n)
        for onset in round.(Int, [0.65, 1.50, 2.10] .* FS)
            width = round(Int, 0.20 * FS)
            for k = 1:width
                idx = onset + k
                if 1 <= idx <= n
                    s2[idx] += 3.5 * (k / 25.0) * exp(-(k - 25.0) / 45.0)
                end
            end
        end

        # S3: muscle burst
        burst_center = 1.80
        burst_σ     = 0.11
        burst_freq   = 22.0
        s3           = 2.5 .* sin.(2π .* burst_freq .* (t_arr .- burst_center)) .* exp.(.-((t_arr .- burst_center) .^ 2) ./ (2 .* burst_σ^2))

        if noise > 0.0
            s1 .+= noise .* randn(n)
            s2 .+= noise .* 0.35 .* randn(n)
            s3 .+= noise .* 0.5 .* randn(n)
        end

        return s1, s2, s3
    end

    # ── 3D rotation ────────────────────────────────────────────────────────────
    #
    # Mixing:   R₂₃(γ) · R₁₃(β) · R₁₂(α) · S
    # Unmixing: R₁₂ᵀ(φ) · R₁₃ᵀ(ψ) · R₂₃ᵀ(χ) · M
    # When φ=α, ψ=β, χ=γ → perfect recovery.
    #
    """Apply 3D mixing rotation R₂₃(γ) · R₁₃(β) · R₁₂(α) to three source signals."""
    function apply_rotation_3d(x, y, z, α_deg, β_deg, γ_deg)
        α, β, γ = deg2rad(α_deg), deg2rad(β_deg), deg2rad(γ_deg)
        x1 = cos(α) .* x .- sin(α) .* y
        y1 = sin(α) .* x .+ cos(α) .* y
        z1 = z
        x2 = cos(β) .* x1 .- sin(β) .* z1
        y2 = y1
        z2 = sin(β) .* x1 .+ cos(β) .* z1
        x3 = x2
        y3 = cos(γ) .* y2 .- sin(γ) .* z2
        z3 = sin(γ) .* y2 .+ cos(γ) .* z2
        return x3, y3, z3
    end

    """Apply inverse 3D rotation to recover source signals from mixtures."""
    function apply_unmixing_3d(x, y, z, φ_deg, ψ_deg, χ_deg)
        φ, ψ, χ = deg2rad(φ_deg), deg2rad(ψ_deg), deg2rad(χ_deg)
        x1 = x
        y1 = cos(χ) .* y .+ sin(χ) .* z
        z1 = -sin(χ) .* y .+ cos(χ) .* z
        x2 = cos(ψ) .* x1 .+ sin(ψ) .* z1
        y2 = y1
        z2 = -sin(ψ) .* x1 .+ cos(ψ) .* z1
        x3 = cos(φ) .* x2 .+ sin(φ) .* y2
        y3 = -sin(φ) .* x2 .+ cos(φ) .* y2
        z3 = z2
        return x3, y3, z3
    end

    # ── Excess kurtosis ────────────────────────────────────────────────────────
    """Compute excess kurtosis (0 for Gaussian, >0 for super-Gaussian signals)."""
    function kurtosis_excess(x)
        μ = mean(x)
        σ = std(x)
        σ < 1e-10 && return 0.0
        return mean(((x .- μ) ./ σ) .^ 4) - 3.0
    end

    # ── Data observables ───────────────────────────────────────────────────────
    src1, src2, src3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))
    mix1, mix2, mix3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))
    rec1, rec2, rec3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))

    kurt_src_text = Observable("Excess kurtosis:  S1=–  S2=–  S3=–")
    kurt_rec_text = Observable("Excess kurtosis:  C1=–  C2=–  C3=–")
    mix_title     = Observable("EEG Recordings (mixtures)  ·  α=25°  β=15°  γ=20°")
    rec_title     = Observable("Recovered Components  ·  φ=25°  ψ=15°  χ=20°")

    # ── Main update ────────────────────────────────────────────────────────────
    """Recompute sources, mixtures, and recovered signals from current parameter values."""
    function update()
        s1, s2, s3 = make_sources(sig_freq[], noise_lvl[])
        m1, m2, m3 = apply_rotation_3d(s1, s2, s3, mix_α[], mix_β[], mix_γ[])

        src1[], src2[], src3[] = s1, s2, s3
        mix1[], mix2[], mix3[] = m1, m2, m3

        if use_infomax[]
            R = ica_rec[]
            c1, c2, c3 = R[1, :], R[2, :], R[3, :]
        else
            c1, c2, c3 = apply_unmixing_3d(m1, m2, m3, unmix_α[], unmix_β[], unmix_γ[])
        end

        rec1[], rec2[], rec3[] = c1, c2, c3

        ks1, ks2, ks3 = round.(kurtosis_excess.([s1, s2, s3]), digits = 1)
        kc1, kc2, kc3 = round.(kurtosis_excess.([c1, c2, c3]), digits = 1)

        kurt_src_text[] = "Excess kurtosis:  S1=$(ks1)  S2=$(ks2)  S3=$(ks3)"
        kurt_rec_text[] = "Excess kurtosis:  C1=$(kc1)  C2=$(kc2)  C3=$(kc3)"
        mix_title[] = "EEG Recordings (mixtures)  ·  α=$(round(Int, mix_α[]))°  β=$(round(Int, mix_β[]))°  γ=$(round(Int, mix_γ[]))°"
        rec_title[] =
            use_infomax[] ? "Recovered by Infomax ICA  (move an unmix slider to return to manual mode)" :
            "Recovered Components  ·  φ=$(round(Int, unmix_α[]))°  ψ=$(round(Int, unmix_β[]))°  χ=$(round(Int, unmix_γ[]))°"
    end

    # Mixing / signal changes invalidate any stored Infomax result
    for obs in (mix_α, mix_β, mix_γ, sig_freq, noise_lvl)
        on(obs) do _
            use_infomax[] = false
            update()
        end
    end
    for obs in (unmix_α, unmix_β, unmix_γ)
        on(obs) do _
            update()
        end
    end

    # ── Colours ────────────────────────────────────────────────────────────────
    col_a = RGBf(0.20, 0.44, 0.69)   # blue  — S1/M1/C1
    col_b = RGBf(0.75, 0.22, 0.17)   # red   — S2/M2/C2
    col_c = RGBf(0.18, 0.63, 0.34)   # green — S3/M3/C3

    # ── Axis helpers ───────────────────────────────────────────────────────────
    """Create a time-series Axis in the specified row."""
    function signal_axis(row, title_obs)
        Axis(
            fig[row, 1];
            title = title_obs,
            xlabel = "Time (s)",
            ylabel = "Amplitude",
            titlesize = title_sz,
            xlabelsize = lbl_sz,
            ylabelsize = lbl_sz,
            xticklabelsize = tick_sz,
            yticklabelsize = tick_sz,
        )
    end

    """Create a scatter-plot Axis for pairwise joint distributions."""
    function scatter_axis(row, col, title_str, xl, yl)
        Axis(
            fig[row, col];
            title = title_str,
            xlabel = xl,
            ylabel = yl,
            titlesize = title_sz,
            xlabelsize = lbl_sz,
            ylabelsize = lbl_sz,
            xticklabelsize = tick_sz,
            yticklabelsize = tick_sz,
        )
    end

    ax_src = signal_axis(1, "True Sources  (S1 = alpha bursts,  S2 = blink,  S3 = muscle burst)")
    ax_mix = signal_axis(2, mix_title)
    ax_rec = signal_axis(3, rec_title)
    linkxaxes!(ax_src, ax_mix, ax_rec)

    ax_sc_src_12 = scatter_axis(1, 2, "Sources S1 vs S2  →  non-blob", "S1", "S2")
    ax_sc_mix_12 = scatter_axis(2, 2, "Mixed M1 vs M2  →  blob", "M1", "M2")
    ax_sc_rec_12 = scatter_axis(3, 2, "Recovered C1 vs C2", "C1", "C2")
    ax_sc_src_13 = scatter_axis(1, 3, "Sources S1 vs S3  →  non-blob", "S1", "S3")
    ax_sc_mix_13 = scatter_axis(2, 3, "Mixed M1 vs M3  →  blob", "M1", "M3")
    ax_sc_rec_13 = scatter_axis(3, 3, "Recovered C1 vs C3", "C1", "C3")
    ax_sc_src_23 = scatter_axis(1, 4, "Sources S2 vs S3  →  non-blob", "S2", "S3")
    ax_sc_mix_23 = scatter_axis(2, 4, "Mixed M2 vs M3  →  blob", "M2", "M3")
    ax_sc_rec_23 = scatter_axis(3, 4, "Recovered C2 vs C3", "C2", "C3")

    # ── Signal plots ───────────────────────────────────────────────────────────
    lines!(ax_src, t_arr, src1; color = col_a, linewidth = 1.5, label = "S1  (alpha bursts)")
    lines!(ax_src, t_arr, src2; color = col_b, linewidth = 1.5, label = "S2  (blink)")
    lines!(ax_src, t_arr, src3; color = col_c, linewidth = 1.5, label = "S3  (muscle burst)")
    axislegend(ax_src; position = :rt, labelsize = lbl_sz)

    lines!(ax_mix, t_arr, mix1; color = col_a, linewidth = 1.5, label = "M1")
    lines!(ax_mix, t_arr, mix2; color = col_b, linewidth = 1.5, label = "M2")
    lines!(ax_mix, t_arr, mix3; color = col_c, linewidth = 1.5, label = "M3")
    axislegend(ax_mix; position = :rt, labelsize = lbl_sz)

    lines!(ax_rec, t_arr, rec1; color = col_a, linewidth = 1.5, label = "C1")
    lines!(ax_rec, t_arr, rec2; color = col_b, linewidth = 1.5, label = "C2")
    lines!(ax_rec, t_arr, rec3; color = col_c, linewidth = 1.5, label = "C3")
    axislegend(ax_rec; position = :rt, labelsize = lbl_sz)

    # ── Scatter plots ──────────────────────────────────────────────────────────
    scatter!(ax_sc_src_12, src1, src2; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_mix_12, mix1, mix2; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_rec_12, rec1, rec2; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_src_13, src1, src3; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_mix_13, mix1, mix3; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_rec_13, rec1, rec3; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_src_23, src2, src3; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_mix_23, mix2, mix3; color = (:black, 0.25), markersize = 3)
    scatter!(ax_sc_rec_23, rec2, rec3; color = (:black, 0.25), markersize = 3)

    text!(
        ax_sc_mix_12,
        0.03,
        0.97;
        text = kurt_src_text,
        space = :relative,
        fontsize = @lift($tick_sz - 1),
        color = :black,
        align = (:left, :top),
    )
    text!(
        ax_sc_rec_12,
        0.03,
        0.97;
        text = kurt_rec_text,
        space = :relative,
        fontsize = @lift($tick_sz - 1),
        color = :black,
        align = (:left, :top),
    )

    # ── Controls ───────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[4, 1:4]; colgap = 12)

    """Create a slider with a header label above and a value label below."""
    function labelled_slider(parent, col, header, range_vals, startval, fmt)
        Label(parent[1, col], header; fontsize = ctrl_sz, halign = :center)
        sl  = Slider(parent[2, col]; range = range_vals, startvalue = startval)
        lbl = Label(parent[3, col], fmt(startval); fontsize = ctrl_sz, halign = :center)
        return sl, lbl
    end

    sl_mix_α, lbl_mix_α     = labelled_slider(ctrl, 1, "Mix α\n(S1-S2)", 0.0:1.0:90.0, 25.0, v -> "$(round(Int, v))°")
    sl_mix_β, lbl_mix_β     = labelled_slider(ctrl, 2, "Mix β\n(S1-S3)", 0.0:1.0:90.0, 15.0, v -> "$(round(Int, v))°")
    sl_mix_γ, lbl_mix_γ     = labelled_slider(ctrl, 3, "Mix γ\n(S2-S3)", 0.0:1.0:90.0, 20.0, v -> "$(round(Int, v))°")
    sl_unmix_α, lbl_unmix_α = labelled_slider(ctrl, 4, "Unmix φ\n(S1-S2)", 0.0:1.0:90.0, 25.0, v -> "$(round(Int, v))°")
    sl_unmix_β, lbl_unmix_β = labelled_slider(ctrl, 5, "Unmix ψ\n(S1-S3)", 0.0:1.0:90.0, 15.0, v -> "$(round(Int, v))°")
    sl_unmix_γ, lbl_unmix_γ = labelled_slider(ctrl, 6, "Unmix χ\n(S2-S3)", 0.0:1.0:90.0, 20.0, v -> "$(round(Int, v))°")
    sl_freq, lbl_freq         = labelled_slider(ctrl, 7, "Burst Freq", 1.0:1.0:20.0, 8.0, v -> "$(round(Int, v)) Hz")
    sl_noise, lbl_noise       = labelled_slider(ctrl, 8, "Noise", 0.0:0.05:1.0, 0.1, v -> "$(round(v; digits = 2))")

    Label(ctrl[1, 9], "Auto-ICA"; fontsize = ctrl_sz, halign = :center)
    btn_ica     = Button(ctrl[2, 9]; label = "Infomax", fontsize = ctrl_sz)
    btn_ica_ext = Button(ctrl[3, 9]; label = "Extended", fontsize = ctrl_sz)

    # ── Wire sliders ───────────────────────────────────────────────────────────
    on(sl_mix_α.value) do v
        mix_α[] = v
        lbl_mix_α.text = "$(round(Int, v))°"
    end
    on(sl_mix_β.value) do v
        mix_β[] = v
        lbl_mix_β.text = "$(round(Int, v))°"
    end
    on(sl_mix_γ.value) do v
        mix_γ[] = v
        lbl_mix_γ.text = "$(round(Int, v))°"
    end
    on(sl_unmix_α.value) do v
        in_ica_update[] || (use_infomax[] = false)
        unmix_α[] = v
        lbl_unmix_α.text = "$(round(Int, v))°"
    end
    on(sl_unmix_β.value) do v
        in_ica_update[] || (use_infomax[] = false)
        unmix_β[] = v
        lbl_unmix_β.text = "$(round(Int, v))°"
    end
    on(sl_unmix_γ.value) do v
        in_ica_update[] || (use_infomax[] = false)
        unmix_γ[] = v
        lbl_unmix_γ.text = "$(round(Int, v))°"
    end
    on(sl_freq.value) do v
        sig_freq[] = v
        lbl_freq.text = "$(round(Int, v)) Hz"
    end
    on(sl_noise.value) do v
        noise_lvl[] = v
        lbl_noise.text = "$(round(v; digits = 2))"
    end

    # ── Auto-ICA ──────────────────────────────────────────────────────────────
    # After running Infomax/Extended, a kurtosis grid search finds the nearest
    # rotation angles so the sliders give visual feedback. The displayed signals
    # still use the Infomax W matrix (stored in ica_rec[]) not the rotation.
    """Run Infomax ICA, store the unmixing result, and find closest rotation angles for slider feedback."""
    function _run_ica_and_apply(ica_fn)
        dat_mat      = vcat(mix1[]', mix2[]', mix3[]')   # 3 × n
        dummy_df     = DataFrame(label = [:M1, :M2, :M3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0])
        dummy_layout = Layout(dummy_df, nothing, nothing, nothing)
        result       = ica_fn(dat_mat, dummy_layout, "demo"; n_components = 3, params = IcaPrms(max_iter = 300))
        dat_c        = (dat_mat .- result.mean) ./ result.scale
        recovered    = result.unmixing * dat_c
        for i = 1:3
            recovered[i, :] ./= max(std(recovered[i, :]), 1e-6)
        end
        ica_rec[]     = recovered
        use_infomax[] = true
        update()

        # Find the closest rotation angles to the Infomax unmixing.
        # These are shown on the sliders as visual feedback — the actual
        # recovery uses the Infomax W matrix stored in ica_rec[].
        m1v, m2v, m3v = mix1[], mix2[], mix3[]
        """Sum of absolute excess kurtosis across three unmixed signals for a given rotation."""
        function kurtosis_score(α, β, γ)
            c1t, c2t, c3t = apply_unmixing_3d(m1v, m2v, m3v, α, β, γ)
            return abs(kurtosis_excess(c1t)) + abs(kurtosis_excess(c2t)) + abs(kurtosis_excess(c3t))
        end
        best_score = -Inf
        best_α, best_β, best_γ = unmix_α[], unmix_β[], unmix_γ[]
        for α = 0.0:5.0:90.0, β = 0.0:5.0:90.0, γ = 0.0:5.0:90.0
            s = kurtosis_score(α, β, γ)
            if s > best_score
                best_score = s
                best_α, best_β, best_γ = α, β, γ
            end
        end
        fine_α, fine_β, fine_γ = best_α, best_β, best_γ
        best_score2 = -Inf
        for α = max(0.0, best_α-6.0):1.0:min(90.0, best_α+6.0),
            β = max(0.0, best_β-6.0):1.0:min(90.0, best_β+6.0),
            γ = max(0.0, best_γ-6.0):1.0:min(90.0, best_γ+6.0)

            s = kurtosis_score(α, β, γ)
            if s > best_score2
                best_score2 = s
                fine_α, fine_β, fine_γ = α, β, γ
            end
        end
        in_ica_update[] = true
        set_close_to!(sl_unmix_α, fine_α)
        set_close_to!(sl_unmix_β, fine_β)
        set_close_to!(sl_unmix_γ, fine_γ)
        in_ica_update[] = false
    end

    on(btn_ica.clicks) do _
        _run_ica_and_apply(infomax_ica)
    end
    on(btn_ica_ext.clicks) do _
        _run_ica_and_apply(infomax_extended_ica)
    end

    # ── Column / row sizing ────────────────────────────────────────────────────
    for col = 1:9
        colsize!(ctrl, col, Relative(1 / 9))
    end

    rowsize!(fig.layout, 1, Relative(0.26))
    rowsize!(fig.layout, 2, Relative(0.26))
    rowsize!(fig.layout, 3, Relative(0.26))
    rowsize!(fig.layout, 4, Relative(0.14))

    colsize!(fig.layout, 1, Relative(0.47))
    colsize!(fig.layout, 2, Relative(0.177))
    colsize!(fig.layout, 3, Relative(0.177))
    colsize!(fig.layout, 4, Relative(0.177))

    update()
    display(fig)
    return fig
end
