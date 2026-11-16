# Validation script for Mass-Univariate Mixed Models

using Test
using EegFun
using MixedModels
using StatsModels
using DataFrames
using Distributions
using Random
using Logging
using LinearAlgebra
using SparseArrays
using CairoMakie

const LMM_EXT = Base.get_extension(EegFun, :EegFunMixedModelsExt)

# Simulated subjects × trials data: random intercepts (and optionally slopes), a within-subject
# `condition` (±0.5), a between-subject `group` (±0.5) and a continuous `rt`.
function _sim_lmm(rng; S = 10, n_s = 20, n_ch = 2, n_tp = 3, slope_sd = 0.0, cond_effect = 0.0)
    n = S * n_s
    meta = DataFrame(
        subject = repeat(string.("S", lpad.(1:S, 2, '0')), inner = n_s),
        condition = repeat([-0.5, 0.5], n ÷ 2),
        group = repeat([s <= S ÷ 2 ? -0.5 : 0.5 for s = 1:S], inner = n_s),
        rt = randn(rng, n),
    )
    sid = repeat(1:S, inner = n_s)
    eeg = zeros(n_ch, n_tp, n)
    for ch = 1:n_ch, tp = 1:n_tp
        b0, b1 = randn(rng, S), slope_sd .* randn(rng, S)
        eeg[ch, tp, :] .= b0[sid] .+ (b1[sid] .+ cond_effect) .* meta.condition .+ randn(rng, n)
    end
    return eeg, meta
end

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

function _mock_epochs(rng; S = 8, n_s = 6, n_channels = 4, n_timepoints = 10)
    dfs = DataFrame[]
    for s = 1:S, k = 1:n_s
        df = DataFrame(randn(rng, n_timepoints, n_channels), [Symbol("Ch$j") for j = 1:n_channels])
        df.time = collect(range(0.0, 0.9, length = n_timepoints))
        df.subject = fill(string("Sub", s), n_timepoints)
        df.condition = fill(isodd(k) ? "A" : "B", n_timepoints)
        df.rt = fill(rand(rng), n_timepoints)
        push!(dfs, df)
    end
    # Valid non-collinear 2D coordinates so topographic interpolation works end-to-end
    incs = [35.0, 35.0, 55.0, 55.0][1:n_channels]
    azis = [-45.0, 45.0, -135.0, 135.0][1:n_channels]
    layout_df = DataFrame(label = [Symbol("Ch$j") for j = 1:n_channels], inc = incs, azi = azis)
    layout = EegFun.Layout(layout_df, nothing, nothing, nothing)
    EegFun._ensure_coordinates_2d!(layout)
    EegFun._ensure_coordinates_3d!(layout)
    return EegFun.EpochData("mock", 1, "mock", dfs, layout, 100, EegFun.AnalysisInfo()), dfs
end

@testset "MixedModels Extension Validation" begin
    rng = MersenneTwister(1)
    epochs, dfs = _mock_epochs(rng)
    f = @formula(amplitude ~ condition + rt + (1 | subject))

    result = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 0))

    # The wrapper must reproduce a single MixedModels fit at one point exactly
    meta_df = DataFrame(
        subject = [df.subject[1] for df in dfs],
        condition = [df.condition[1] for df in dfs],
        rt = [df.rt[1] for df in dfs],
        amplitude = [df[1, :Ch1] for df in dfs],
    )
    m_single = quiet(() -> fit(MixedModel, f, meta_df; progress = false))
    @test isapprox(result.beta[1, 1, :], coef(m_single), atol = 1e-5)
    @test isapprox(result.z_values[1, 1, :], coeftable(m_single).cols[3], atol = 1e-5)
    @test isnothing(result.p_corrected)
    @test result.inference.correction === :none

    # Between-Within parametric degrees of freedom:
    # S = 8 subjects, N = 48 epochs, 1 between col (Intercept) -> df = 7; 2 within cols (condition: B, rt) -> df = 48 - 8 - 2 = 38
    @test result.inference.df_parametric == [7.0, 38.0, 38.0]
    for k = 1:3
        z11 = result.z_values[1, 1, k]
        @test result.p_parametric[1, 1, k] ≈ 2.0 * ccdf(TDist(result.inference.df_parametric[k]), abs(z11)) atol = 1e-12
    end

    # Accessors and show
    @test EegFun.channel_labels(result) == [Symbol("Ch$j") for j = 1:4]
    @test EegFun.n_channels(result) == 4
    @test EegFun.n_samples(result) == 10
    @test length(EegFun.time_vector(result)) == 10
    str = sprint(show, result)
    @test occursin("LmmStatsResult", str) && occursin("4 channels", str) && occursin("Wald z only", str)

    # With permutations
    res_p = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 99, rng = MersenneTwister(2)))
    plain = sprint(show, MIME"text/plain"(), res_p)
    @test occursin("sign flipping by :subject", plain)
    @test occursin("Significant points (corrected", plain)
    @test res_p.epochs === epochs
end

@testset "fit_mass_lmm: closed-form pivotal CR2 sign-flip null equals dense GLS + CR2" begin
    rng = MersenneTwister(3)
    eeg, meta = _sim_lmm(rng; S = 9, n_s = 16, n_ch = 1, n_tp = 1, slope_sd = 0.7)
    meta.amplitude = eeg[1, 1, :]
    for (form, _) in (
        (@formula(amplitude ~ 1 + condition + rt + (1 + condition | subject)), "slopes"),
        (@formula(amplitude ~ 1 + condition + group + (1 | subject)), "between"),
    )
        m = quiet(() -> fit(MixedModel, form, meta; progress = false))
        X = Matrix{Float64}(modelmatrix(m)[:, invperm(MixedModels.pivot(m))])
        p = size(X, 2)
        ids = LMM_EXT._block_ids(meta.subject)
        design = LMM_EXT.BlockDesign(m, X, ids)
        XtVX_s, XtVy_s, _ = LMM_EXT._block_gls_stats(design, [Matrix(A.λ) for A in m.reterms], m.y)
        S = maximum(ids)
        H = hcat(ones(S), rand(rng, [-1.0, 1.0], S, 30))

        # Dense N×N reference: build full V = I + G G' at θ̂, run dense Freedman–Lane + block CR2
        G_re = zeros(length(m.y), 0)
        for s = 1:S
            Gs = zeros(length(m.y), size(design.Z[s], 2))
            Λ = zeros(size(design.Z[s], 2), size(design.Z[s], 2))
            for (k, off) in design.lambda_blocks[s]
                q = size(m.reterms[k].λ, 1)
                Λ[(off+1):(off+q), (off+1):(off+q)] .= Matrix(m.reterms[k].λ)
            end
            Gs[design.rows[s], :] .= design.Z[s] * Λ
            G_re = hcat(G_re, Gs)
        end
        Vi = inv(I + G_re * G_re')
        XtVX_dense = X' * Vi * X
        chol_d = cholesky(Symmetric(XtVX_dense))
        M_d = inv(chol_d)
        R_d = Matrix(chol_d.U) \ Matrix{Float64}(I, p, p)

        for c = 1:p
            z_null = zeros(size(H, 2))
            z_obs = LMM_EXT._cr2_pivotal_z!(z_null, XtVX_s, XtVy_s, c, H)

            # Identity column H[:, 1] == ones(S) must match z_obs to machine precision
            @test isapprox(z_null[1], z_obs; atol = 1e-12)

            # Dense Freedman–Lane pseudo-observations y*(η) = X_{-c} γ̂ + diag(η) r_{-c}
            idx = setdiff(1:p, c)
            Xr = X[:, idx]
            γ = (Xr' * Vi * Xr) \ (Xr' * Vi * m.y)
            r = m.y - Xr * γ

            # Dense block CR2 weights and Satterthwaite df
            vs_dense = Vector{Vector{Float64}}(undef, S)
            Omega_dense = zeros(S, S)
            cs_list = Vector{Vector{Float64}}(undef, S)
            Hs_list = Vector{Matrix{Float64}}(undef, S)
            for s = 1:S
                rs = design.rows[s]
                Vis = inv(Matrix(I + G_re[rs, :] * G_re[rs, :]'))
                Hs_d = Matrix(Symmetric(R_d' * X[rs, :]' * Vis * X[rs, :] * R_d))
                eig = eigen(Symmetric(Hs_d))
                As_d = eig.vectors * Diagonal(1.0 ./ sqrt.(max.(1.0 .- eig.values, 1e-8))) * eig.vectors'
                cs_d = As_d * R_d[c, :]
                Hs_list[s] = Hs_d
                cs_list[s] = cs_d
                vs_dense[s] = R_d * cs_d
            end
            for s = 1:S, t = 1:S
                Omega_dense[s, t] =
                    (s == t ? dot(cs_list[s], Hs_list[s] * cs_list[s]) : 0.0) - dot(Hs_list[s] * cs_list[s], Hs_list[t] * cs_list[t])
            end
            df_d = clamp(tr(Omega_dense)^2 / tr(Omega_dense^2), 1.0, Float64(S - 1))

            z_dense = map(1:size(H, 2)) do col
                ys = Xr * γ + H[ids, col] .* r
                β = M_d * (X' * Vi * ys)
                e = ys - X * β
                ss = 0.0
                for s = 1:S
                    rs = design.rows[s]
                    Vis = inv(Matrix(I + G_re[rs, :] * G_re[rs, :]'))
                    gs = dot(vs_dense[s], X[rs, :]' * (Vis * e[rs]))
                    ss += gs^2
                end
                LMM_EXT._t_to_z(β[c] / sqrt(ss), df_d)
            end
            @test maximum(abs, z_null .- z_dense) < 1e-9
        end
    end
end

@testset "fit_mass_lmm: exchangeability blocks and sign-flip enumeration" begin
    rng = MersenneTwister(4)
    eeg, meta = _sim_lmm(rng; S = 6, n_s = 10)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))

    # 2^(6-1) = 32 ≤ n_perms → exact enumeration; the identity is in the null set
    res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 1000))
    @test res.inference.exhaustive
    @test res.inference.n_perms == 32
    @test res.inference.flip_by === :subject
    @test res.inference.n_blocks == 6
    ic = findfirst(==("condition"), res.coefficients)
    @test any(isapprox(maximum(abs, res.inference.z_perm_obs[:, :, ic]); rtol = 1e-8), res.null_max[:, ic])
    @test minimum(res.p_corrected[:, :, ic]) >= 1 / 32

    # --- Sign matrix tests ---
    H, exh = LMM_EXT._sign_flip_matrix(4, 100, rng)
    @test exh && size(H) == (4, 8) && all(H[1, :] .== 1) && length(unique(eachcol(H))) == 8

    all_vecs4 = [[(b >> i) & 1 == 1 ? 1.0 : -1.0 for i = 0:3] for b = 0:15]
    canonical4 = [v[1] > 0 ? v : -v for v in all_vecs4]
    classes4 = unique(canonical4)
    @test length(classes4) == 8
    @test count(==(ones(4)), canonical4) == 2

    non_id_vecs4 = filter(v -> !all(==(1.0), v) && !all(==(-1.0), v), all_vecs4)
    @test length(non_id_vecs4) == 14
    non_id_classes4 = unique([v[1] > 0 ? v : -v for v in non_id_vecs4])
    @test length(non_id_classes4) == 7
    @test !any(==(ones(4)), non_id_classes4)

    # Monte Carlo rejection sampler: deterministic exclusion of identity patterns
    H, exh = LMM_EXT._sign_flip_matrix(12, 100, rng)
    @test !exh && size(H) == (12, 100)
    @test !any(all(==(1.0), col) for col in eachcol(H))
    @test !any(all(==(-1.0), col) for col in eachcol(H))

    # Reproducible seeded test of the rejection sampler
    H_seed1, exh_seed1 = LMM_EXT._sign_flip_matrix(6, 10, MersenneTwister(42))
    H_seed2, _ = LMM_EXT._sign_flip_matrix(6, 10, MersenneTwister(42))
    @test !exh_seed1 && size(H_seed1) == (6, 10)
    @test H_seed1 == H_seed2
    @test !any(all(==(1.0), col) for col in eachcol(H_seed1))
    @test !any(all(==(-1.0), col) for col in eachcol(H_seed1))

    # Supplementary goodness-of-fit check across quotient classes
    cols_mc = [col[1] > 0 ? col : -col for _ = 1:500 for col in eachcol(LMM_EXT._sign_flip_matrix(6, 30, rng)[1])]
    unique_mc = unique(cols_mc)
    @test length(unique_mc) == 31
    counts_mc = [count(==(u), cols_mc) for u in unique_mc]
    expected_mc = length(cols_mc) / 31
    chisq_mc = sum((counts_mc .- expected_mc) .^ 2 ./ expected_mc)
    @test chisq_mc < 65.0

    # Too few blocks
    eeg4, meta4 = _sim_lmm(rng; S = 4, n_s = 10)
    @test_throws EegFun.EegFunError quiet(() -> EegFun.fit_mass_lmm(eeg4, meta4, f; n_perms = 100))
    @test quiet(() -> EegFun.fit_mass_lmm(eeg4, meta4, f; n_perms = 0)) isa EegFun.LmmStatsResult

    # Nested subject/session: flip_by is the outer factor
    meta.session = repeat(repeat(["a", "b"], inner = 5), 6)
    res_n = quiet(
        () ->
            EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + (1 | subject) + (1 | subject & session)); n_perms = 50),
    )
    @test res_n.inference.flip_by === :subject

    # Crossed subject × item: no exact scheme → error with n_perms > 0
    meta.item = string.(repeat(1:10, 6))
    f_cross = @formula(amplitude ~ 1 + condition + (1 | subject) + (1 | item))
    @test_throws EegFun.EegFunError quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_cross; n_perms = 50))
    @test_throws EegFun.EegFunError quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_cross; n_perms = 50, flip_by = :subject))
    @test quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_cross; n_perms = 0)) isa EegFun.LmmStatsResult

    # Unknown flip_by column
    @test_throws EegFun.EegFunError quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 50, flip_by = :nope))
end

@testset "fit_mass_lmm: between-subject predictor" begin
    rng = MersenneTwister(5)
    eeg, meta = _sim_lmm(rng; S = 12, n_s = 10)
    res = quiet(
        () -> EegFun.fit_mass_lmm(
            eeg,
            meta,
            @formula(amplitude ~ 1 + condition + group + (1 | subject));
            n_perms = 200,
            rng = MersenneTwister(6),
        ),
    )
    ig = findfirst(==("group"), res.coefficients)
    @test res.inference.df_parametric[ig] == 10.0  # S - 2 = 10 exact two-sample df
    @test all(0 .< res.p_corrected[:, :, ig] .<= 1)
    @test res.tested == [findfirst(==("condition"), res.coefficients), ig]
end

@testset "fit_mass_lmm: rank-deficient designs" begin
    rng = Random.MersenneTwister(11)
    n_sub, n_per = 10, 30
    n = n_sub * n_per
    meta = DataFrame(Subject = repeat(string.("S", 1:n_sub), inner = n_per))
    meta.a = randn(rng, n)
    meta.b = 2 .* meta.a              # exactly collinear with a -> rank deficient
    meta.c = randn(rng, n)
    sub_eff = Dict(s => randn(rng) for s in unique(meta.Subject))

    n_ch, n_tp = 2, 3
    eeg = zeros(n_ch, n_tp, n)
    for ch = 1:n_ch, tp = 1:n_tp
        eeg[ch, tp, :] .= 0.4 .* meta.b .- 0.3 .* meta.c .+ [sub_eff[s] for s in meta.Subject] .+ randn(rng, n)
    end

    f_rd = @formula(amplitude ~ 1 + a + b + c + (1 | Subject))
    f_fr = @formula(amplitude ~ 1 + b + c + (1 | Subject))

    probe = copy(meta)
    probe.amplitude = randn(rng, n)
    m_probe = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f_rd, probe)
    end
    @test m_probe.feterm.rank == 3
    @test MixedModels.pivot(m_probe) != 1:4

    res_rd =
        @test_logs (:warn, r"rank deficient") match_mode = :any EegFun.fit_mass_lmm(eeg, meta, f_rd; n_perms = 50, rng = MersenneTwister(5))
    res_fr = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_fr; n_perms = 50, rng = MersenneTwister(5)))

    ia = findfirst(==("a"), res_rd.coefficients)
    for fld in (:beta, :se, :z_values, :p_parametric, :p_corrected)
        @test all(isnan, getfield(res_rd, fld)[:, :, ia])
    end
    @test all(isnan, res_rd.null_max[:, ia])
    @test ia ∉ res_rd.tested

    for name in ("(Intercept)", "b", "c")
        i_rd = findfirst(==(name), res_rd.coefficients)
        i_fr = findfirst(==(name), res_fr.coefficients)
        @test res_rd.beta[:, :, i_rd] ≈ res_fr.beta[:, :, i_fr] atol = 1e-8
        @test res_rd.z_values[:, :, i_rd] ≈ res_fr.z_values[:, :, i_fr] atol = 1e-8
    end
    for name in ("b", "c")
        i_rd = findfirst(==(name), res_rd.coefficients)
        i_fr = findfirst(==(name), res_fr.coefficients)
        @test res_rd.null_max[:, i_rd] ≈ res_fr.null_max[:, i_fr] rtol = 1e-5
        @test res_rd.p_corrected[:, :, i_rd] == res_fr.p_corrected[:, :, i_fr]
    end
end

@testset "fit_mass_lmm: input validation and error handling" begin
    empty_epochs = EegFun.EpochData(
        "empty",
        1,
        "empty",
        DataFrame[],
        EegFun.Layout(DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]), nothing, nothing, nothing),
        100,
        EegFun.AnalysisInfo(),
    )
    f_dummy = @formula(amplitude ~ condition + (1 | subject))
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(empty_epochs, f_dummy)

    eeg_bad = zeros(2, 2, 5)
    meta_bad = DataFrame(a = 1:10)
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg_bad, meta_bad, @formula(amplitude ~ a))

    eeg, meta = _sim_lmm(MersenneTwister(7))
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + (1 | subject)); correction = :bogus)

    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(
        eeg,
        meta,
        @formula(amplitude ~ 1 + condition + (1 | subject));
        dim_order = (:channels, :channels, :epochs),
    )
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(
        eeg,
        meta,
        @formula(amplitude ~ 1 + condition + (1 | subject));
        dim_order = (:channels, :timepoints),
    )

    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(
        eeg,
        meta,
        @formula(amplitude ~ 1 + condition + (1 | subject));
        channel_names = [:OnlyOne],
    )
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + (1 | subject)); time_points = [1.0])
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(
        eeg,
        meta,
        @formula(amplitude ~ 1 + condition + (1 | subject));
        time_points = [1.0, NaN, 3.0],
    )

    eeg_flat = copy(eeg)
    eeg_flat[1, 1, :] .= 0.0
    res_flat = quiet(() -> EegFun.fit_mass_lmm(eeg_flat, meta, @formula(amplitude ~ 1 + condition + (1 | subject)); n_perms = 0))
    @test res_flat.singular_fits isa BitMatrix || res_flat.singular_fits isa Array{Bool}
    @test res_flat.singular_fits[1, 1] == false
    @test res_flat.fit_status[1, 1] == LMM_EXT.STATUS_FLAT
end

@testset "fit_mass_lmm: cluster and TFCE corrections" begin
    rng = MersenneTwister(42)
    eeg, meta = _sim_lmm(rng; S = 10, n_s = 8, n_ch = 3, n_tp = 4)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    conn = sparse([1, 2, 2, 3], [2, 1, 3, 2], trues(4), 3, 3)
    for corr in (:cluster, :tfce)
        res = quiet(
            () -> EegFun.fit_mass_lmm(
                eeg,
                meta,
                f;
                n_perms = 20,
                correction = corr,
                cluster_threshold = 1.0,
                spatial_connectivity = conn,
                rng = MersenneTwister(1),
            ),
        )
        ic = findfirst(==("condition"), res.coefficients)
        ii = findfirst(==("(Intercept)"), res.coefficients)
        @test res.inference.correction === corr
        @test size(res.null_max) == (20, 2)
        @test all(res.null_max[:, ic] .>= 0.0)
        @test all(isnan, res.null_max[:, ii])
        @test all(0.0 .< res.p_corrected[:, :, ic] .<= 1.0)
        @test all(isnan, res.p_corrected[:, :, ii])
    end
end

@testset "fit_mass_lmm: extract_predictor_stats, StatsAPI, and plotting integration" begin
    rng = MersenneTwister(8)
    epochs, _ = _mock_epochs(rng; S = 8, n_s = 6)
    f = @formula(amplitude ~ condition + (1 | subject))
    conn = sparse(Bool.(I(4)))
    for corr in (:maxt, :cluster, :tfce)
        res = quiet(
            () -> EegFun.fit_mass_lmm(
                epochs,
                f;
                n_perms = 30,
                correction = corr,
                cluster_threshold = 1.5,
                spatial_connectivity = conn,
                rng = MersenneTwister(3),
            ),
        )
        pr = quiet(() -> EegFun.extract_predictor_stats(res, "condition: B"))
        @test pr isa EegFun.PermutationResult
        @test pr.stat_matrix.t == res.inference.z_perm_obs[:, :, 2]
        # Extracted ERP pair encodes [β̂, 0], so StatsAPI.coef(pr) == β̂
        @test EegFun.StatsAPI.coef(pr) ≈ res.beta[:, :, 2] atol = 1e-12
        if corr !== :cluster
            @test pr.masks.positive == ((res.p_corrected[:, :, 2] .<= 0.05) .& (res.inference.z_perm_obs[:, :, 2] .> 0))
            @test EegFun.StatsAPI.pvalue(pr) == res.p_corrected[:, :, 2]
        end
        @test pr.critical_t == (1.5, -1.5)

        # End-to-end plotting smoke tests (no display)
        plt_erp = EegFun.plot_erp_stats(res; coef = "condition: B", display_plot = false)
        @test plt_erp.fig isa Figure
        plt_hm = EegFun.plot_stat_heatmap(res; coef = 2, display_plot = false)
        @test plt_hm.fig isa Figure
        plt_topo_t = EegFun.plot_topography_stats(res; coef = 2, topo_data = :tvalues, n_topos = 2, display_plot = false)
        @test plt_topo_t.fig isa Figure
        plt_topo_b = EegFun.plot_topography_stats(res; coef = "condition: B", topo_data = :beta, n_topos = 2, display_plot = false)
        @test plt_topo_b.fig isa Figure
    end

    # Flat/ineligible pixel consistency in extract_predictor_stats (:cluster)
    epochs_bad, _ = _mock_epochs(MersenneTwister(19); S = 8, n_s = 6)
    for df in epochs_bad.data
        df[2, :Ch1] = 0.0
    end
    res_bad = quiet(
        () -> EegFun.fit_mass_lmm(
            epochs_bad,
            f;
            n_perms = 30,
            correction = :cluster,
            cluster_threshold = 0.5,
            spatial_connectivity = conn,
            rng = MersenneTwister(4),
        ),
    )
    pr_bad = EegFun.extract_predictor_stats(res_bad, 2; alpha = 0.5)
    @test !pr_bad.masks.positive[1, 2] && !pr_bad.masks.negative[1, 2]

    res = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 30, rng = MersenneTwister(3)))
    @test_throws EegFun.EegFunError EegFun.extract_predictor_stats(res, "(Intercept)")
    @test_throws EegFun.EegFunError EegFun.extract_predictor_stats(res, "nope")

    @test EegFun.StatsAPI.pvalue(res) === res.p_corrected
    @test EegFun.StatsAPI.pvalue(res; corrected = false) === res.p_parametric
    @test EegFun.StatsAPI.nobs(res) == 48
    @test EegFun.StatsAPI.coef(res) === res.beta
    ci_lo, ci_hi = EegFun.StatsAPI.confint(res; level = 0.95)
    @test size(ci_lo) == size(res.beta) && all(ci_lo .< ci_hi)
    res0 = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 0))
    @test EegFun.StatsAPI.pvalue(res0) === res0.p_parametric
    @test_throws EegFun.EegFunError EegFun.StatsAPI.pvalue(res0; corrected = true)
end

@testset "fit_mass_lmm: multi-pixel FWER null calibration (random slopes, misspecified, between-subject)" begin
    rng = MersenneTwister(2026)
    n_sim = 80

    # 1. Misspecified random-intercept model when true data has random slopes (20 pixels per dataset!)
    n_rej_fwer_mis, n_rej_pix1_wald = 0, 0
    f_int = @formula(amplitude ~ 1 + condition + (1 | subject))
    for _ = 1:n_sim
        eeg, meta = _sim_lmm(rng; S = 12, n_s = 20, n_ch = 5, n_tp = 4, slope_sd = 0.7)
        res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_int; n_perms = 149, rng = rng, progress = false))
        n_rej_fwer_mis += any(<=(0.05), @view res.p_corrected[:, :, 2])
        n_rej_pix1_wald += (res.p_parametric[1, 1, 2] <= 0.05)
    end
    @test n_rej_fwer_mis / n_sim <= 0.125
    @test n_rej_pix1_wald > n_rej_fwer_mis

    # 2. Random-slopes model (1 + condition | subject) across 20 pixels (previously ~47% FWER under frozen-θ!)
    n_rej_fwer_slp = 0
    f_slp = @formula(amplitude ~ 1 + condition + (1 + condition | subject))
    for _ = 1:n_sim
        eeg, meta = _sim_lmm(rng; S = 12, n_s = 20, n_ch = 5, n_tp = 4, slope_sd = 0.7)
        res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_slp; n_perms = 149, rng = rng, progress = false))
        n_rej_fwer_slp += any(<=(0.05), @view res.p_corrected[:, :, 2])
    end
    @test n_rej_fwer_slp / n_sim <= 0.125

    # 3. Between-subject group predictor across 20 pixels (previously ~60% FWER under frozen-θ!)
    n_rej_fwer_grp = 0
    f_grp = @formula(amplitude ~ 1 + condition + group + (1 | subject))
    for _ = 1:n_sim
        eeg, meta = _sim_lmm(rng; S = 12, n_s = 16, n_ch = 5, n_tp = 4, slope_sd = 0.0)
        res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_grp; tested_coefs = :group, n_perms = 149, rng = rng, progress = false))
        ig = findfirst(==("group"), res.coefficients)
        n_rej_fwer_grp += any(<=(0.05), @view res.p_corrected[:, :, ig])
    end
    @test n_rej_fwer_grp / n_sim <= 0.125
end

@testset "fit_mass_lmm: tested_coefs name matching" begin
    names = ["(Intercept)", "condition: B", "condition: C", "rt", "condition: B & rt", "condition: C & rt"]
    @test LMM_EXT._resolve_tested_coefs(names, :condition, false) == [2, 3]
    @test LMM_EXT._resolve_tested_coefs(names, :rt, false) == [4]
    @test LMM_EXT._resolve_tested_coefs(names, "condition: B", false) == [2]
    @test LMM_EXT._resolve_tested_coefs(names, "condition: B & rt", false) == [5]
    @test_throws EegFun.EegFunError LMM_EXT._resolve_tested_coefs(names, :nonexistent, false)
end

@testset "fit_mass_lmm: flat channel and bad pixels are isolated (NaN, one warning)" begin
    rng = MersenneTwister(21)
    eeg, meta = _sim_lmm(rng; S = 10, n_s = 12, n_ch = 4, n_tp = 4)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    conn = sparse(Bool.(I(3)))
    eeg_flat = copy(eeg)
    eeg_flat[2, :, :] .= 0.0
    keep = [1, 3, 4]
    ic = 2

    for corr in (:maxt, :cluster, :tfce)
        kw = (n_perms = 30, correction = corr, cluster_threshold = 1.0)
        res = @test_logs (:warn, r"could not be fitted") match_mode = :any EegFun.fit_mass_lmm(
            eeg_flat,
            meta,
            f;
            kw...,
            spatial_connectivity = sparse(Bool.(I(4))),
            rng = MersenneTwister(5),
        )
        ref = quiet(() -> EegFun.fit_mass_lmm(eeg[keep, :, :], meta, f; kw..., spatial_connectivity = conn, rng = MersenneTwister(5)))
        @test all(==(LMM_EXT.STATUS_FLAT), res.fit_status[2, :])
        for fld in (:beta, :se, :z_values, :p_parametric, :p_corrected)
            @test all(isnan, getfield(res, fld)[2, :, ic])
        end
        @test res.beta[keep, :, :] ≈ ref.beta atol = 1e-10
        @test res.z_values[keep, :, :] ≈ ref.z_values atol = 1e-10
        @test res.null_max[:, ic] ≈ ref.null_max[:, ic] atol = 1e-6
        @test res.p_corrected[keep, :, ic] == ref.p_corrected[:, :, ic]
        @test !any(isnan, res.p_corrected[keep, :, ic])
    end

    eeg_bad = copy(eeg)
    eeg_bad[1, 2, :] .= 3.0
    eeg_bad[3, 4, 7] = NaN
    res = @test_logs (:warn, r"2 of 16 channel") match_mode = :any EegFun.fit_mass_lmm(
        eeg_bad,
        meta,
        f;
        n_perms = 30,
        rng = MersenneTwister(5),
    )
    flat = res.fit_status .== LMM_EXT.STATUS_FLAT
    @test findall(flat) == [CartesianIndex(1, 2), CartesianIndex(3, 4)]
    @test all(isnan, res.z_values[1, 2, :]) && all(isnan, res.z_values[3, 4, :])
    @test all(isnan, res.p_corrected[1, 2, ic]) && all(isnan, res.p_corrected[3, 4, ic])
    @test !any(isnan, res.z_values[.!flat, ic])
    @test !any(isnan, res.p_corrected[.!flat, ic])
    @test occursin("flat_response", sprint(show, res))
end

@testset "fit_mass_lmm: convergence classification and reporting" begin
    C = LMM_EXT
    @test C._classify_returnvalue(:FTOL_REACHED) == C.STATUS_CONVERGED
    @test C._classify_returnvalue(:XTOL_REACHED) == C.STATUS_CONVERGED
    @test C._classify_returnvalue(:SUCCESS) == C.STATUS_CONVERGED
    @test C._classify_returnvalue(:ROUNDOFF_LIMITED) == C.STATUS_ROUNDOFF
    @test C._classify_returnvalue(:MAXEVAL_REACHED) == C.STATUS_MAXEVAL
    @test C._classify_returnvalue(:FAILURE) == C.STATUS_OPTFAIL
    @test C._classify_returnvalue(:SOMETHING_NEW) == C.STATUS_OPTFAIL

    rng = MersenneTwister(22)
    eeg, meta = _sim_lmm(rng; S = 8, n_s = 10, n_ch = 2, n_tp = 3)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 0))
    @test all(<=(C.STATUS_ROUNDOFF), res.fit_status)
    @test size(res.fit_status) == (2, 3)

    meta_i = copy(meta)
    meta_i.amplitude = eeg[1, 1, :]
    scr = C._create_scratch(f, meta_i, [1, 2], nothing, 0, 0, Dict{Symbol,Any}())
    scr.m.optsum.maxfeval = 2
    scr.y .= eeg[1, 1, :]
    @test C._compute_pixel!(scr, (REML = false, n_flips = 0, est_cols = [1, 2], tested = Int[], tested_pos = Int[])) == C.STATUS_MAXEVAL

    status = fill(C.STATUS_CONVERGED, 3, 4)
    status[1, :] .= C.STATUS_FLAT
    status[2, 2] = C.STATUS_ERROR
    status[3, 1] = C.STATUS_MAXEVAL
    status[3, 3] = C.STATUS_OPTFAIL
    @test_logs (:warn, r"(?s)5 of 12.*Ch1 \(4/4 time points\).*Ch2 \(1/4.*boom.*2 of 12 fits did not converge") C._report_fit_problems(
        status,
        [:Ch1, :Ch2, :Ch3],
        "boom",
    )
    @test_logs C._report_fit_problems(fill(C.STATUS_CONVERGED, 2, 2), [:A, :B], nothing)
    @test_logs C._report_fit_problems(fill(C.STATUS_ROUNDOFF, 2, 2), [:A, :B], nothing)
end

@testset "fit_mass_lmm: power (injected effect is detected, null channel is not)" begin
    rng = MersenneTwister(23)
    eeg, meta = _sim_lmm(rng; S = 12, n_s = 20, n_ch = 3, n_tp = 4)
    eeg[1, 2, :] .+= 1.2 .* meta.condition
    eeg[1, 3, :] .+= 1.2 .* meta.condition
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    conn = sparse(Bool.(I(3)))

    for corr in (:maxt, :cluster, :tfce)
        res = quiet(
            () -> EegFun.fit_mass_lmm(
                eeg,
                meta,
                f;
                n_perms = 199,
                correction = corr,
                cluster_threshold = 1.5,
                spatial_connectivity = conn,
                rng = MersenneTwister(10),
            ),
        )
        ic = findfirst(==("condition"), res.coefficients)
        @test res.p_corrected[1, 2, ic] <= 0.05
        @test res.p_corrected[1, 3, ic] <= 0.05
        @test all(res.p_corrected[2:3, :, ic] .> 0.05)
    end
end

@testset "fit_mass_lmm: REML estimation" begin
    rng = MersenneTwister(24)
    eeg, meta = _sim_lmm(rng; S = 8, n_s = 10, n_ch = 2, n_tp = 2)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))

    res_reml = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; REML = true, n_perms = 30, rng = MersenneTwister(1)))
    @test res_reml.inference.REML == true

    meta_11 = copy(meta)
    meta_11.amplitude = eeg[1, 1, :]
    m_reml = quiet(() -> fit(MixedModel, f, meta_11; REML = true, progress = false))
    @test isapprox(res_reml.beta[1, 1, :], coef(m_reml), atol = 1e-5)
    @test isapprox(res_reml.se[1, 1, :], stderror(m_reml), atol = 1e-5)
    @test isapprox(res_reml.z_values[1, 1, :], coeftable(m_reml).cols[3], atol = 1e-5)
end

@testset "fit_mass_lmm: dim_order" begin
    rng = MersenneTwister(25)
    eeg, meta = _sim_lmm(rng; S = 8, n_s = 10, n_ch = 2, n_tp = 3)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))

    res_std = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 0))

    eeg_p1 = permutedims(eeg, (2, 1, 3))
    res_p1 = quiet(() -> EegFun.fit_mass_lmm(eeg_p1, meta, f; dim_order = (:timepoints, :channels, :epochs), n_perms = 0))
    @test res_std.beta == res_p1.beta
    @test res_std.z_values == res_p1.z_values

    eeg_p2 = permutedims(eeg, (3, 1, 2))
    res_p2 = quiet(() -> EegFun.fit_mass_lmm(eeg_p2, meta, f; dim_order = (:epochs, :channels, :timepoints), n_perms = 0))
    @test res_std.beta == res_p2.beta
    @test res_std.z_values == res_p2.z_values

    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg, meta, f; dim_order = (:channels, :timepoints, :bad))
end

@testset "fit_mass_lmm: contrasts" begin
    rng = MersenneTwister(26)
    S, n_s = 10, 15
    n = S * n_s
    meta = DataFrame(
        subject = repeat(string.("S", lpad.(1:S, 2, '0')), inner = n_s),
        cond3 = repeat(["a", "b", "c"], n ÷ 3),
        amplitude = randn(rng, n),
    )
    eeg = reshape(meta.amplitude, 1, 1, n)
    f = @formula(amplitude ~ 1 + cond3 + (1 | subject))

    contrasts = Dict(:cond3 => EffectsCoding(base = "a"))
    res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; contrasts = contrasts, n_perms = 30, rng = MersenneTwister(1)))
    @test res.inference.contrasts == contrasts

    m_direct = quiet(() -> fit(MixedModel, f, meta; contrasts = contrasts, progress = false))
    @test isapprox(res.beta[1, 1, :], coef(m_direct), atol = 1e-5)
    @test isapprox(res.se[1, 1, :], stderror(m_direct), atol = 1e-5)

    res_term =
        quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; contrasts = contrasts, tested_coefs = :cond3, n_perms = 30, rng = MersenneTwister(1)))
    @test length(res_term.tested) == 2
    @test all(0 .< res_term.p_corrected[1, 1, res_term.tested] .<= 1)
end

@testset "fit_mass_lmm: memory guard and small n_perms warning" begin
    C = LMM_EXT
    @test_logs (:warn, r"The cluster/TFCE null needs") C._check_zperm_memory(4, 10, 100, 1; free = 20_000)
    @test_throws EegFun.EegFunError C._check_zperm_memory(4, 10, 100, 1; free = 17_000)
    @test_logs C._check_zperm_memory(4, 10, 100, 1; free = 1_000_000)

    @test_logs (:warn, r"no point can reach p < 0.05") C._warn_small_n_perms(15, false)
    @test_logs (:warn, r"resolution of") C._warn_small_n_perms(50, false)
    @test_logs C._warn_small_n_perms(16, true)
    @test_logs C._warn_small_n_perms(100, false)

    rng = MersenneTwister(27)
    eeg, meta = _sim_lmm(rng; S = 12, n_s = 6, n_ch = 1, n_tp = 1)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    @test_logs (:warn, r"Only 15 random sign flips") match_mode = :any EegFun.fit_mass_lmm(
        eeg,
        meta,
        f;
        n_perms = 15,
        rng = MersenneTwister(1),
    )
end

@testset "fit_mass_lmm: multi-chunk thread execution" begin
    rng = MersenneTwister(28)
    eeg, meta = _sim_lmm(rng; S = 6, n_s = 8, n_ch = 3, n_tp = 4)
    eeg[2, 2, :] .= 0.0
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    conn = sparse(Bool.(I(3)))

    for corr in (:maxt, :cluster, :tfce)
        res = quiet(
            () -> EegFun.fit_mass_lmm(
                eeg,
                meta,
                f;
                n_perms = 32,
                correction = corr,
                spatial_connectivity = conn,
                rng = MersenneTwister(7),
                progress = false,
            ),
        )
        ic = 2
        @test res.fit_status[2, 2] == LMM_EXT.STATUS_FLAT
        @test isnan(res.z_values[2, 2, ic])
        @test isnan(res.p_corrected[2, 2, ic])
        @test all(isfinite, res.z_values[1, :, ic])
        @test all(0 .< res.p_corrected[1, :, ic] .<= 1)
    end
end

@testset "fit_mass_lmm: all-pixel constant data is safely handled" begin
    rng = MersenneTwister(99)
    eeg, meta = _sim_lmm(rng; S = 8, n_s = 6, n_ch = 2, n_tp = 2)
    eeg .= 0.0
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    res = @test_logs (:warn, r"4 of 4 channel") match_mode = :any EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 30, rng = MersenneTwister(1))
    @test all(==(LMM_EXT.STATUS_FLAT), res.fit_status)
    @test all(isnan, res.z_values)
    @test all(isnan, res.p_corrected)
end

@testset "_lmm_perm_pvalue: edge cases, exact ties, and zeros" begin
    null = [1.0, 2.0, 3.0, 4.0]
    @test EegFun._lmm_perm_pvalue(null, 4.0, true) == 1 / 4
    @test EegFun._lmm_perm_pvalue(null, 3.0, true) == 2 / 4
    @test EegFun._lmm_perm_pvalue(null, 2.0, true) == 3 / 4
    @test EegFun._lmm_perm_pvalue(null, 1.0, true) == 4 / 4
    @test EegFun._lmm_perm_pvalue(null, 0.5, true) == 1.0
    @test EegFun._lmm_perm_pvalue(null, 5.0, true) == 0.0

    @test EegFun._lmm_perm_pvalue(null, 4.0, false) == (1 + 1) / (4 + 1)
    @test EegFun._lmm_perm_pvalue(null, 5.0, false) == (0 + 1) / (4 + 1)

    null_near = [1.0, 2.0, 3.0 - 1e-11, 4.0]
    @test EegFun._lmm_perm_pvalue(null_near, 3.0, true) == 2 / 4

    null_zeros = [0.0, 0.0, 0.0, 1.5, 3.0]
    @test EegFun._lmm_perm_pvalue(null_zeros, 0.0, true) == 1.0
    @test EegFun._lmm_perm_pvalue(null_zeros, 0.0, false) == 1.0
end
