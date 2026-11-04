"""
    signal_example_ica_geometry()

Interactive 3-source / 3-sensor ICA blind source separation demo illustrating pairwise 3D
rotations, joint-distribution scatter geometry, excess kurtosis, and Infomax / Extended Infomax
unmixing.

## See Also
- [`signal_example_ica_mixture`](@ref), [`signal_example_ica_separation`](@ref),
  [`signal_example_ica_clt`](@ref), [`signal_example_ica_sphering`](@ref),
  [`signal_example_ica_optimization`](@ref)

# Examples
```julia
using GLMakie, EegFun
EegFun.signal_example_ica_geometry()
```

# Returns
- `fig::Figure`: The Makie figure object
"""
function signal_example_ica_geometry()
    _check_interactive_demo_backend("ICA Geometry demo") || return nothing

    FS = 512.0
    T = 3.0

    fig = _eegfun_figure(size = (1400, 850), figure_padding = (8, 8, 8, 8))
    lbl_sz, tick_sz, title_sz, ctrl_sz = _demo_adaptive_fonts(fig, 1400, (13, 10), (11, 8), (14, 11), (13, 10))

    t_arr = collect(range(0.0, T - 1.0 / FS, step = 1.0 / FS))
    n = length(t_arr)

    noise1_buf = randn(n)
    noise2_buf = randn(n)
    noise3_buf = randn(n)

    mix_α, mix_β, mix_γ = Observable(25.0), Observable(15.0), Observable(20.0)
    unmix_α, unmix_β, unmix_γ = Observable(25.0), Observable(15.0), Observable(20.0)
    sig_freq = Observable(8.0)
    noise_lvl = Observable(0.1)

    use_infomax = Observable(false)
    ica_rec = Observable(zeros(3, n))
    in_ica_update = Ref(false)

    # Precompute static S2 (blinks) and S3 (muscle Gabor) base templates
    s2_base = zeros(Float64, n)
    blink_width = round(Int, 0.20 * FS)
    for onset in round.(Int, (0.65, 1.50, 2.10) .* FS)
        @inbounds for k = 1:blink_width
            idx = onset + k
            if 1 <= idx <= n
                s2_base[idx] += 3.5 * (k / 25.0) * exp(-(k - 25.0) / 45.0)
            end
        end
    end
    s3_base = @. 2.5 * sin(2π * 22.0 * (t_arr - 1.80)) * exp(-((t_arr - 1.80)^2) / (2.0 * 0.11^2))
    s1_env = [@. 2.0 * exp(-((t_arr - c)^2) / (2.0 * 0.09^2)) for c in (0.25, 1.05, 2.55)]

    function make_sources(freq::Float64, noise::Float64)
        s1 = zeros(Float64, n)
        for (c, env) in zip((0.25, 1.05, 2.55), s1_env)
            @. s1 += sin(2π * freq * (t_arr - c)) * env
        end
        s2 = copy(s2_base)
        s3 = copy(s3_base)
        if noise > 0.0
            @. s1 += noise * noise1_buf
            @. s2 += noise * 0.35 * noise2_buf
            @. s3 += noise * 0.5 * noise3_buf
        end
        return s1, s2, s3
    end

    function _mixing_matrix(α_deg::Float64, β_deg::Float64, γ_deg::Float64)
        sα, cα = sincos(deg2rad(α_deg))
        sβ, cβ = sincos(deg2rad(β_deg))
        sγ, cγ = sincos(deg2rad(γ_deg))
        return (
            cβ * cα,
            -cβ * sα,
            -sβ,
            cγ * sα - sγ * sβ * cα,
            cγ * cα + sγ * sβ * sα,
            -sγ * cβ,
            sγ * sα + cγ * sβ * cα,
            sγ * cα - cγ * sβ * sα,
            cγ * cβ,
        )
    end

    function _apply_mat_3d(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, m::NTuple{9,Float64})
        m11, m12, m13, m21, m22, m23, m31, m32, m33 = m
        o1 = Vector{Float64}(undef, n)
        o2 = Vector{Float64}(undef, n)
        o3 = Vector{Float64}(undef, n)
        @inbounds for i = 1:n
            xi, yi, zi = x[i], y[i], z[i]
            o1[i] = m11 * xi + m12 * yi + m13 * zi
            o2[i] = m21 * xi + m22 * yi + m23 * zi
            o3[i] = m31 * xi + m32 * yi + m33 * zi
        end
        return o1, o2, o3
    end

    function apply_rotation_3d(x, y, z, α_deg, β_deg, γ_deg)
        return _apply_mat_3d(x, y, z, _mixing_matrix(α_deg, β_deg, γ_deg))
    end

    function apply_unmixing_3d(x, y, z, φ_deg, ψ_deg, χ_deg)
        m11, m12, m13, m21, m22, m23, m31, m32, m33 = _mixing_matrix(φ_deg, ψ_deg, χ_deg)
        return _apply_mat_3d(x, y, z, (m11, m21, m31, m12, m22, m32, m13, m23, m33))
    end

    function kurtosis_excess(x::AbstractVector{Float64})
        len = length(x)
        μ = mean(x)
        m2 = 0.0
        m4 = 0.0
        @inbounds for i = 1:len
            d = x[i] - μ
            d2 = d * d
            m2 += d2
            m4 += d2 * d2
        end
        σ2 = m2 / (len - 1)
        σ2 < 1e-20 && return 0.0
        return (m4 / len) / (σ2 * σ2) - 3.0
    end

    function _kurtosis_score_unmix(
        x::Vector{Float64},
        y::Vector{Float64},
        z::Vector{Float64},
        φ_deg::Float64,
        ψ_deg::Float64,
        χ_deg::Float64,
    )
        # x, y, z are mean-centered, so every rotated projection has zero mean
        m11, m12, m13, m21, m22, m23, m31, m32, m33 = _mixing_matrix(φ_deg, ψ_deg, χ_deg)
        s1_2 = s1_4 = 0.0
        s2_2 = s2_4 = 0.0
        s3_2 = s3_4 = 0.0
        @inbounds for i = 1:n
            xi, yi, zi = x[i], y[i], z[i]
            c1 = m11 * xi + m21 * yi + m31 * zi
            c2 = m12 * xi + m22 * yi + m32 * zi
            c3 = m13 * xi + m23 * yi + m33 * zi
            c1_2 = c1 * c1
            c2_2 = c2 * c2
            c3_2 = c3 * c3
            s1_2 += c1_2
            s1_4 += c1_2 * c1_2
            s2_2 += c2_2
            s2_4 += c2_2 * c2_2
            s3_2 += c3_2
            s3_4 += c3_2 * c3_2
        end
        inv_n = 1.0 / n
        inv_nm1 = 1.0 / (n - 1)
        score = 0.0
        for (m2, m4) in ((s1_2, s1_4), (s2_2, s2_4), (s3_2, s3_4))
            σ2 = m2 * inv_nm1
            if σ2 > 1e-20
                score += abs((m4 * inv_n) / (σ2 * σ2) - 3.0)
            end
        end
        return score
    end

    src1, src2, src3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))
    mix1, mix2, mix3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))
    rec1, rec2, rec3 = Observable(zeros(n)), Observable(zeros(n)), Observable(zeros(n))

    kurt_src_text = Observable("Excess kurtosis:  S1=–  S2=–  S3=–")
    kurt_rec_text = Observable("Excess kurtosis:  C1=–  C2=–  C3=–")
    mix_title = Observable("EEG Recordings (mixtures)  ·  α=25°  β=15°  γ=20°")
    rec_title = Observable("Recovered Components  ·  φ=25°  ψ=15°  χ=20°")

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

        ks1, ks2, ks3 =
            round(kurtosis_excess(s1); digits = 1), round(kurtosis_excess(s2); digits = 1), round(kurtosis_excess(s3); digits = 1)
        kc1, kc2, kc3 =
            round(kurtosis_excess(c1); digits = 1), round(kurtosis_excess(c2); digits = 1), round(kurtosis_excess(c3); digits = 1)

        kurt_src_text[] = "Excess kurtosis:  S1=$(ks1)  S2=$(ks2)  S3=$(ks3)"
        kurt_rec_text[] = "Excess kurtosis:  C1=$(kc1)  C2=$(kc2)  C3=$(kc3)"
        mix_title[] = "EEG Recordings (mixtures)  ·  α=$(round(Int, mix_α[]))°  β=$(round(Int, mix_β[]))°  γ=$(round(Int, mix_γ[]))°"
        rec_title[] =
            use_infomax[] ? "Recovered by Infomax ICA  (move an unmix slider to return to manual mode)" :
            "Recovered Components  ·  φ=$(round(Int, unmix_α[]))°  ψ=$(round(Int, unmix_β[]))°  χ=$(round(Int, unmix_γ[]))°"
    end

    for obs in (mix_α, mix_β, mix_γ, sig_freq, noise_lvl)
        on(obs) do _
            use_infomax[] = false
            update()
        end
    end
    for obs in (unmix_α, unmix_β, unmix_γ)
        on(_ -> update(), obs)
    end

    col_a = RGBf(0.20, 0.44, 0.69)
    col_b = RGBf(0.75, 0.22, 0.17)
    col_c = RGBf(0.18, 0.63, 0.34)

    ax_common = (titlesize = title_sz, xlabelsize = lbl_sz, ylabelsize = lbl_sz, xticklabelsize = tick_sz, yticklabelsize = tick_sz)
    signal_axis(row, title_obs) = Axis(fig[row, 1]; title = title_obs, xlabel = "Time (s)", ylabel = "Amplitude", ax_common...)
    scatter_axis(row, col, title_str, xl, yl) = Axis(fig[row, col]; title = title_str, xlabel = xl, ylabel = yl, ax_common...)

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

    for (ax, a, b, c, la, lb, lc) in (
        (ax_src, src1, src2, src3, "S1  (alpha bursts)", "S2  (blink)", "S3  (muscle burst)"),
        (ax_mix, mix1, mix2, mix3, "M1", "M2", "M3"),
        (ax_rec, rec1, rec2, rec3, "C1", "C2", "C3"),
    )
        lines!(ax, t_arr, a; color = col_a, linewidth = 1.5, label = la)
        lines!(ax, t_arr, b; color = col_b, linewidth = 1.5, label = lb)
        lines!(ax, t_arr, c; color = col_c, linewidth = 1.5, label = lc)
        axislegend(ax; position = :rt, labelsize = lbl_sz)
    end

    for (ax, x_obs, y_obs) in (
        (ax_sc_src_12, src1, src2),
        (ax_sc_mix_12, mix1, mix2),
        (ax_sc_rec_12, rec1, rec2),
        (ax_sc_src_13, src1, src3),
        (ax_sc_mix_13, mix1, mix3),
        (ax_sc_rec_13, rec1, rec3),
        (ax_sc_src_23, src2, src3),
        (ax_sc_mix_23, mix2, mix3),
        (ax_sc_rec_23, rec2, rec3),
    )
        scatter!(ax, x_obs, y_obs; color = (:black, 0.25), markersize = 3)
    end

    annot_sz = @lift($tick_sz - 1)
    text!(ax_sc_mix_12, 0.03, 0.97; text = kurt_src_text, space = :relative, fontsize = annot_sz, color = :black, align = (:left, :top))
    text!(ax_sc_rec_12, 0.03, 0.97; text = kurt_rec_text, space = :relative, fontsize = annot_sz, color = :black, align = (:left, :top))

    ctrl = GridLayout(fig[4, 1:4]; colgap = 12)
    deg_fmt = v -> "$(round(Int, v))°"
    sl_mix_α, lbl_mix_α = _demo_labelled_slider(ctrl, 1, "Mix α\n(S1-S2)", 0.0:1.0:90.0, 25.0, deg_fmt, ctrl_sz)
    sl_mix_β, lbl_mix_β = _demo_labelled_slider(ctrl, 2, "Mix β\n(S1-S3)", 0.0:1.0:90.0, 15.0, deg_fmt, ctrl_sz)
    sl_mix_γ, lbl_mix_γ = _demo_labelled_slider(ctrl, 3, "Mix γ\n(S2-S3)", 0.0:1.0:90.0, 20.0, deg_fmt, ctrl_sz)
    sl_unmix_α, lbl_unmix_α = _demo_labelled_slider(ctrl, 4, "Unmix φ\n(S1-S2)", 0.0:1.0:90.0, 25.0, deg_fmt, ctrl_sz)
    sl_unmix_β, lbl_unmix_β = _demo_labelled_slider(ctrl, 5, "Unmix ψ\n(S1-S3)", 0.0:1.0:90.0, 15.0, deg_fmt, ctrl_sz)
    sl_unmix_γ, lbl_unmix_γ = _demo_labelled_slider(ctrl, 6, "Unmix χ\n(S2-S3)", 0.0:1.0:90.0, 20.0, deg_fmt, ctrl_sz)
    sl_freq, lbl_freq = _demo_labelled_slider(ctrl, 7, "Burst Freq", 1.0:1.0:20.0, 8.0, v -> "$(round(Int, v)) Hz", ctrl_sz)
    sl_noise, lbl_noise = _demo_labelled_slider(ctrl, 8, "Noise", 0.0:0.05:1.0, 0.1, v -> "$(round(v; digits = 2))", ctrl_sz)

    Label(ctrl[1, 9], "Auto-ICA"; fontsize = ctrl_sz, halign = :center)
    btn_ica = Button(ctrl[2, 9]; label = "Infomax", fontsize = ctrl_sz)
    btn_ica_ext = Button(ctrl[3, 9]; label = "Extended", fontsize = ctrl_sz)

    for (sl, obs, lbl) in ((sl_mix_α, mix_α, lbl_mix_α), (sl_mix_β, mix_β, lbl_mix_β), (sl_mix_γ, mix_γ, lbl_mix_γ))
        let o = obs, l = lbl
            on(sl.value) do v
                o[] = v
                l.text = deg_fmt(v)
            end
        end
    end
    for (sl, obs, lbl) in ((sl_unmix_α, unmix_α, lbl_unmix_α), (sl_unmix_β, unmix_β, lbl_unmix_β), (sl_unmix_γ, unmix_γ, lbl_unmix_γ))
        let o = obs, l = lbl
            on(sl.value) do v
                in_ica_update[] || (use_infomax[] = false)
                o[] = v
                l.text = deg_fmt(v)
            end
        end
    end
    on(sl_freq.value) do v
        sig_freq[] = v
        lbl_freq.text = "$(round(Int, v)) Hz"
    end
    on(sl_noise.value) do v
        noise_lvl[] = v
        lbl_noise.text = "$(round(v; digits = 2))"
        randn!(noise1_buf)
        randn!(noise2_buf)
        randn!(noise3_buf)
    end

    dummy_layout = Layout(DataFrame(label = [:M1, :M2, :M3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]), nothing, nothing, nothing)

    function _run_ica_and_apply(ica_fn)
        dat_mat = vcat(mix1[]', mix2[]', mix3[]')
        result = ica_fn(dat_mat, dummy_layout, "demo"; n_components = 3, params = IcaPrms(max_iter = 300))
        recovered = result.unmixing * ((dat_mat .- result.mean) ./ result.scale)
        for i = 1:3
            @views recovered[i, :] ./= max(std(recovered[i, :]), 1e-6)
        end
        ica_rec[] = recovered
        use_infomax[] = true
        update()

        m1c = mix1[] .- mean(mix1[])
        m2c = mix2[] .- mean(mix2[])
        m3c = mix3[] .- mean(mix3[])
        best_score = -Inf
        best_α, best_β, best_γ = unmix_α[], unmix_β[], unmix_γ[]
        for α = 0.0:5.0:90.0, β = 0.0:5.0:90.0, γ = 0.0:5.0:90.0
            s = _kurtosis_score_unmix(m1c, m2c, m3c, α, β, γ)
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

            s = _kurtosis_score_unmix(m1c, m2c, m3c, α, β, γ)
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

    on(_ -> _run_ica_and_apply(infomax_ica), btn_ica.clicks)
    on(_ -> _run_ica_and_apply(infomax_extended_ica), btn_ica_ext.clicks)

    for col = 1:9
        colsize!(ctrl, col, Relative(1 / 9))
    end
    for r = 1:3
        rowsize!(fig.layout, r, Relative(0.26))
    end
    rowsize!(fig.layout, 4, Relative(0.14))
    colsize!(fig.layout, 1, Relative(0.47))
    for c = 2:4
        colsize!(fig.layout, c, Relative(0.177))
    end

    update()
    _display_demo_figure(fig)
    return fig
end
