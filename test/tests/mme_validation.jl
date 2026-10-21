# Validation script for Mass-Univariate Mixed Models

using Test
using EegFun
using MixedModels
using StatsModels
using DataFrames
using Random
using Logging
using LinearAlgebra
using SparseArrays

const LMM_EXT = Base.get_extension(EegFun, :EegFunMixedModelsExt)

# Simulated subjects × trials data: random intercepts (and optionally slopes), a within-subject
# `condition` (±0.5), a between-subject `group` (±0.5) and a continuous `rt`.
function _sim_lmm(rng; S = 10, n_s = 20, n_ch = 2, n_tp = 3, slope_sd = 0.0, cond_effect = 0.0)
    n = S * n_s
    meta = DataFrame(
        subject = repeat(string.("S", lpad.(1:S, 2, '0')), inner = n_s),
        condition = repeat([-0.5, 0.5], n ÷ 2),
        group = repeat([s <= S ÷ 2 ? -0.5 : 0.5 for s in 1:S], inner = n_s),
        rt = randn(rng, n),
    )
    sid = repeat(1:S, inner = n_s)
    eeg = zeros(n_ch, n_tp, n)
    for ch in 1:n_ch, tp in 1:n_tp
        b0, b1 = randn(rng, S), slope_sd .* randn(rng, S)
        eeg[ch, tp, :] .= b0[sid] .+ (b1[sid] .+ cond_effect) .* meta.condition .+ randn(rng, n)
    end
    return eeg, meta
end

quiet(f) = Logging.with_logger(f, Logging.NullLogger())

function _mock_epochs(rng; S = 8, n_s = 6, n_channels = 4, n_timepoints = 10)
    dfs = DataFrame[]
    for s in 1:S, k in 1:n_s
        df = DataFrame(randn(rng, n_timepoints, n_channels), [Symbol("Ch$j") for j in 1:n_channels])
        df.time = range(0, 1, length = n_timepoints)
        df.subject = fill(string("Sub", s), n_timepoints)
        df.condition = fill(isodd(k) ? "A" : "B", n_timepoints)
        df.rt = fill(rand(rng), n_timepoints)
        push!(dfs, df)
    end
    layout = EegFun.Layout(
        DataFrame(label = [Symbol("Ch$j") for j in 1:n_channels], x = zeros(n_channels), y = zeros(n_channels),
            z = zeros(n_channels), theta = zeros(n_channels), radius = zeros(n_channels)),
        nothing, nothing, nothing)
    return EegFun.EpochData("mock", 1, "mock", dfs, layout, 100, EegFun.AnalysisInfo()), dfs
end

@testset "MixedModels Extension Validation" begin
    rng = MersenneTwister(1)
    epochs, dfs = _mock_epochs(rng)
    f = @formula(amplitude ~ condition + rt + (1 | subject))

    result = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 0))

    # The wrapper must reproduce a single MixedModels fit at one point exactly
    meta_df = DataFrame(subject = [df.subject[1] for df in dfs], condition = [df.condition[1] for df in dfs],
        rt = [df.rt[1] for df in dfs], amplitude = [df[1, :Ch1] for df in dfs])
    m_single = quiet(() -> fit(MixedModel, f, meta_df; progress = false))
    @test isapprox(result.beta[1, 1, :], coef(m_single), atol = 1e-5)
    @test isapprox(result.z_values[1, 1, :], coeftable(m_single).cols[3], atol = 1e-5)
    @test isnothing(result.p_corrected)
    @test result.inference.correction === :none

    # Accessors and show
    @test EegFun.channel_labels(result) == [Symbol("Ch$j") for j in 1:4]
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

@testset "fit_mass_lmm: closed-form sign-flip null equals dense GLS" begin
    rng = MersenneTwister(3)
    eeg, meta = _sim_lmm(rng; S = 9, n_s = 16, n_ch = 1, n_tp = 1, slope_sd = 0.7)
    meta.amplitude = eeg[1, 1, :]
    for (form, label) in ((@formula(amplitude ~ 1 + condition + rt + (1 + condition | subject)), "slopes"),
                          (@formula(amplitude ~ 1 + condition + group + (1 | subject)), "between"))
        m = quiet(() -> fit(MixedModel, form, meta; progress = false))
        X = Matrix{Float64}(modelmatrix(m)[:, invperm(MixedModels.pivot(m))])
        ids = LMM_EXT._block_ids(meta.subject)
        design = LMM_EXT.BlockDesign(m, X, ids)
        stats = LMM_EXT._block_gls_stats(design, [Matrix(A.λ) for A in m.reterms], m.y)
        S = maximum(ids)
        H = hcat(ones(S), rand(rng, [-1.0, 1.0], S, 30))
        ν = MixedModels.ssqdenom(m)
        # Dense reference: V = I + G G' at θ̂, Freedman–Lane with reduced-model GLS nuisance
        G = zeros(length(m.y), 0)
        for s in 1:S
            Gs = zeros(length(m.y), size(design.Z[s], 2))
            Λ = zeros(size(design.Z[s], 2), size(design.Z[s], 2))
            for (k, off) in design.lambda_blocks[s]
                q = size(m.reterms[k].λ, 1)
                Λ[(off + 1):(off + q), (off + 1):(off + q)] .= Matrix(m.reterms[k].λ)
            end
            Gs[design.rows[s], :] .= design.Z[s] * Λ
            G = hcat(G, Gs)
        end
        Vi = inv(I + G * G')
        for c in 1:size(X, 2)
            z = zeros(size(H, 2))
            LMM_EXT._null_z!(z, stats..., c, H, ν)
            idx = setdiff(1:size(X, 2), c)
            Xr = X[:, idx]
            γ = (Xr' * Vi * Xr) \ (Xr' * Vi * m.y)
            r = m.y - Xr * γ
            M = inv(X' * Vi * X)
            P = Vi - Vi * X * M * X' * Vi
            zd = map(1:size(H, 2)) do p
                ys = Xr * γ + H[ids, p] .* r
                β = M * X' * Vi * ys
                β[c] / sqrt(dot(ys, P * ys) / ν * M[c, c])
            end
            @test maximum(abs, z .- zd) < 1e-9
            @test z[1] ≈ coef(m)[c] / stderror(m)[c] atol = 1e-9   # identity flip = the LMM's own z
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
    @test any(isapprox(maximum(abs, res.z_values[:, :, ic]); rtol = 1e-8), res.null_max[:, ic])
    @test minimum(res.p_corrected[:, :, ic]) >= 1 / 32

    # Sign matrix
    H, exh = LMM_EXT._sign_flip_matrix(4, 100, rng)
    @test exh && size(H) == (4, 8) && all(H[1, :] .== 1) && length(unique(eachcol(H))) == 8
    H, exh = LMM_EXT._sign_flip_matrix(12, 100, rng)
    @test !exh && size(H) == (12, 100)

    # Too few blocks
    eeg4, meta4 = _sim_lmm(rng; S = 4, n_s = 10)
    @test_throws EegFun.EegFunError quiet(() -> EegFun.fit_mass_lmm(eeg4, meta4, f; n_perms = 100))
    # ... but Wald-only works
    @test quiet(() -> EegFun.fit_mass_lmm(eeg4, meta4, f; n_perms = 0)) isa EegFun.LmmStatsResult

    # Nested subject/session: flip_by is the outer factor
    meta.session = repeat(repeat(["a", "b"], inner = 5), 6)
    res_n = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + (1 | subject) + (1 | subject & session)); n_perms = 50))
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
    res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + group + (1 | subject));
        n_perms = 200, rng = MersenneTwister(6)))
    ig = findfirst(==("group"), res.coefficients)
    # Whole-block flips change the sign of the group effect: the null is centred, not collapsed onto z_obs
    @test all(0 .< res.p_corrected[:, :, ig] .<= 1)
    @test minimum(res.null_max[:, ig]) < 0.5 * maximum(abs, res.z_values[:, :, ig]) || maximum(abs, res.z_values[:, :, ig]) < 1
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
    for ch in 1:n_ch, tp in 1:n_tp
        eeg[ch, tp, :] .= 0.4 .* meta.b .- 0.3 .* meta.c .+ [sub_eff[s] for s in meta.Subject] .+ randn(rng, n)
    end

    f_rd = @formula(amplitude ~ 1 + a + b + c + (1 | Subject))
    f_fr = @formula(amplitude ~ 1 + b + c + (1 | Subject))

    # Guard: this design must produce a NON-identity pivot (the case the old code got wrong)
    probe = copy(meta); probe.amplitude = randn(rng, n)
    m_probe = Logging.with_logger(Logging.NullLogger()) do
        LinearMixedModel(f_rd, probe)
    end
    @test m_probe.feterm.rank == 3
    @test MixedModels.pivot(m_probe) != 1:4

    # Same RNG seed and number of blocks → identical sign-flip matrices for both fits
    res_rd = @test_logs (:warn, r"rank deficient") match_mode = :any EegFun.fit_mass_lmm(
        eeg, meta, f_rd; n_perms = 50, rng = MersenneTwister(5))
    res_fr = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f_fr; n_perms = 50, rng = MersenneTwister(5)))

    ia = findfirst(==("a"), res_rd.coefficients)
    # Aliased coefficient: every output is NaN, and it is not permuted
    for fld in (:beta, :se, :z_values, :p_parametric, :p_corrected)
        @test all(isnan, getfield(res_rd, fld)[:, :, ia])
    end
    @test all(isnan, res_rd.null_max[:, ia])
    @test ia ∉ res_rd.tested

    # Estimable coefficients: identical to the equivalent full-rank model
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
    # Empty EpochData
    empty_epochs = EegFun.EpochData(
        "empty", 1, "empty", DataFrame[],
        EegFun.Layout(DataFrame(label=[:Ch1]), nothing, nothing, nothing),
        100, EegFun.AnalysisInfo()
    )
    f_dummy = @formula(amplitude ~ condition + (1 | subject))
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(empty_epochs, f_dummy)

    # Mismatched meta_df and eeg tensor
    eeg_bad = zeros(2, 2, 5)
    meta_bad = DataFrame(a = 1:10)
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg_bad, meta_bad, @formula(amplitude ~ a))

    # Unknown correction
    eeg, meta = _sim_lmm(MersenneTwister(7))
    @test_throws EegFun.EegFunError EegFun.fit_mass_lmm(eeg, meta, @formula(amplitude ~ 1 + condition + (1 | subject)); correction = :bogus)
end

@testset "fit_mass_lmm: cluster and TFCE corrections" begin
    rng = MersenneTwister(42)
    eeg, meta = _sim_lmm(rng; S = 10, n_s = 8, n_ch = 3, n_tp = 4)
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    conn = sparse([1, 2, 2, 3], [2, 1, 3, 2], trues(4), 3, 3)
    for corr in (:cluster, :tfce)
        res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 20, correction = corr, cluster_threshold = 1.0,
            spatial_connectivity = conn, rng = MersenneTwister(1)))
        ic = findfirst(==("condition"), res.coefficients)
        ii = findfirst(==("(Intercept)"), res.coefficients)
        @test res.inference.correction === corr
        @test size(res.null_max) == (20, 2)
        @test all(res.null_max[:, ic] .>= 0.0)
        @test all(isnan, res.null_max[:, ii])
        @test all(0.0 .< res.p_corrected[:, :, ic] .<= 1.0)
        # Untested coefficient: no permutation p-value
        @test all(isnan, res.p_corrected[:, :, ii])
    end
end

@testset "fit_mass_lmm: extract_predictor_stats and StatsAPI" begin
    rng = MersenneTwister(8)
    epochs, _ = _mock_epochs(rng; S = 8, n_s = 6)
    f = @formula(amplitude ~ condition + (1 | subject))
    conn = sparse(Bool.(I(4)))
    for corr in (:maxt, :cluster, :tfce)
        res = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 30, correction = corr, cluster_threshold = 1.5,
            spatial_connectivity = conn, rng = MersenneTwister(3)))
        pr = quiet(() -> EegFun.extract_predictor_stats(res, "condition: B"))
        @test pr isa EegFun.PermutationResult
        @test pr.stat_matrix.t == res.z_values[:, :, 2]
        if corr !== :cluster
            @test pr.masks.positive == ((res.p_corrected[:, :, 2] .<= 0.05) .& (res.z_values[:, :, 2] .> 0))
        end
        # Settings come from the fit, not from function defaults
        @test pr.critical_t == (1.5, -1.5)
    end
    res = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 30, rng = MersenneTwister(3)))
    @test_throws EegFun.EegFunError EegFun.extract_predictor_stats(res, "(Intercept)")   # not tested
    @test_throws EegFun.EegFunError EegFun.extract_predictor_stats(res, "nope")

    @test EegFun.StatsAPI.pvalue(res) === res.p_corrected
    @test EegFun.StatsAPI.pvalue(res; corrected = false) === res.p_parametric
    @test EegFun.StatsAPI.nobs(res) == 48
    @test EegFun.StatsAPI.coef(res) === res.beta
    res0 = quiet(() -> EegFun.fit_mass_lmm(epochs, f; n_perms = 0))
    @test EegFun.StatsAPI.pvalue(res0) === res0.p_parametric
    @test_throws EegFun.EegFunError EegFun.StatsAPI.pvalue(res0; corrected = true)
end

@testset "fit_mass_lmm: null calibration (random slopes, misspecified model)" begin
    # Data with random slopes, model with random intercepts only: the Wald test is badly liberal here,
    # whole-subject sign flipping should not be. Loose bound: 120 datasets, α = 0.05.
    rng = MersenneTwister(2026)
    n_rej, n_rej_wald, n_sim = 0, 0, 120
    f = @formula(amplitude ~ 1 + condition + (1 | subject))
    for _ in 1:n_sim
        eeg, meta = _sim_lmm(rng; S = 12, n_s = 20, n_ch = 1, n_tp = 1, slope_sd = 0.6)
        res = quiet(() -> EegFun.fit_mass_lmm(eeg, meta, f; n_perms = 199, rng = rng))
        n_rej += res.p_corrected[1, 1, 2] <= 0.05
        n_rej_wald += res.p_parametric[1, 1, 2] <= 0.05
    end
    @test n_rej / n_sim <= 0.11
    @test n_rej_wald > n_rej
end

@testset "fit_mass_lmm: tested_coefs name matching" begin
    names = ["(Intercept)", "condition: B", "condition: C", "rt", "condition: B & rt", "condition: C & rt"]
    # Term name selects its contrast columns but not interactions
    @test LMM_EXT._resolve_tested_coefs(names, :condition, false) == [2, 3]
    # Exact name: only that column (previously substring matching also captured the interactions)
    @test LMM_EXT._resolve_tested_coefs(names, :rt, false) == [4]
    @test LMM_EXT._resolve_tested_coefs(names, "condition: B", false) == [2]
    # Interaction requested explicitly
    @test LMM_EXT._resolve_tested_coefs(names, "condition: B & rt", false) == [5]
    @test_throws EegFun.EegFunError LMM_EXT._resolve_tested_coefs(names, :nonexistent, false)
end
