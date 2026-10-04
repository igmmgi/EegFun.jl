# Validation script for Mass-Univariate Mixed Models

using Test
using EegFun
using MixedModels
using StatsModels
using DataFrames
using BenchmarkTools
using Random
using Logging

@testset "MixedModels Extension Validation" begin
    # 1. Create a dummy dataset
    n_epochs = 20
    n_channels = 4
    n_timepoints = 10
    
    # We create dummy EegFun DataFrames
    dfs = DataFrame[]
    for i in 1:n_epochs
        # Random data
        df = DataFrame(randn(n_timepoints, n_channels), :auto)
        rename!(df, [Symbol("Ch$j") for j in 1:n_channels])
        
        # Add metadata
        df.time = range(0, 1, length=n_timepoints)
        df.subject = fill(string("Sub", mod(i, 5)), n_timepoints) # 5 subjects
        df.condition = fill(i % 2 == 0 ? "A" : "B", n_timepoints)
        df.rt = fill(rand(), n_timepoints)
        
        push!(dfs, df)
    end
    
    # Create EpochData mock (assuming Layout and other metadata are minimal)
    # EegFun.EpochData typically takes more arguments, so we mock it.
    epochs = EegFun.EpochData(
        "mock", 
        1, 
        "mock", 
        dfs, 
        EegFun.Layout(
            DataFrame(
                label = [Symbol("Ch$j") for j in 1:4],
                x = zeros(4), y = zeros(4), z = zeros(4),
                theta = zeros(4), radius = zeros(4)
            ),
            nothing,
            nothing,
            nothing
        ), 
        100, 
        EegFun.AnalysisInfo()
    )
    
    f = @formula(amplitude ~ condition + rt + (1 | subject))
    
    # 2. Run the Mass-Univariate wrapper
    @info "Running fit_mass_lmm..."
    result = EegFun.fit_mass_lmm(epochs, f)
    
    # 3. Assertions (Internal Correctness)
    # We manually test one single point (Ch1, Timepoint 1) to ensure the 
    # multithreaded extraction and refit! logic works mathematically perfectly.
    
    y = [df[1, :Ch1] for df in dfs]
    meta_df = DataFrame(
        subject = [df.subject[1] for df in dfs],
        condition = [df.condition[1] for df in dfs],
        rt = [df.rt[1] for df in dfs]
    )
    meta_df.amplitude = y
    
    m_single = fit(MixedModel, f, meta_df)
    
    # Extract the coefficients and compare
    expected_coefs = coef(m_single)
    expected_t = coeftable(m_single).cols[3]
    
    # Get the wrapper's result for Ch1 (index 1), Timepoint 1
    actual_coefs = result.beta[1, 1, :]
    actual_t = result.t_values[1, 1, :]
    
    @test isapprox(actual_coefs, expected_coefs, atol=1e-5)
    @test isapprox(actual_t, expected_t, atol=1e-5)
    
    @info "Mathematical equivalence verified!"
    
    # 4. Benchmarking
    # @btime fit_mass_lmm($epochs, $f)
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

    perms = EegFun.generate_permutation_matrix(meta, :Subject; n_perms = 50, rng = Random.MersenneTwister(5))

    # :full (tight tolerance) so the comparison tests the pivot handling rather than optimizer stopping noise
    res_rd = @test_logs (:warn, r"rank deficient") match_mode = :any EegFun.fit_mass_lmm(
        eeg, meta, f_rd; n_perms = 50, perm_matrix = perms)
    res_fr = EegFun.fit_mass_lmm(eeg, meta, f_fr; n_perms = 50, perm_matrix = perms)

    ia = findfirst(==("a"), res_rd.coef_names)
    # Aliased coefficient: every output is NaN, and it is not permuted
    for fld in (:beta, :se, :t, :p, :p_uncorrected, :p_corrected)
        @test all(isnan, getfield(res_rd, fld)[:, :, ia])
    end
    @test all(iszero, res_rd.max_t_null[:, ia])

    # Estimable coefficients: identical to the equivalent full-rank model
    for name in ("(Intercept)", "b", "c")
        i_rd = findfirst(==(name), res_rd.coef_names)
        i_fr = findfirst(==(name), res_fr.coef_names)
        @test res_rd.beta[:, :, i_rd] ≈ res_fr.beta[:, :, i_fr] atol = 1e-8
        @test res_rd.t[:, :, i_rd] ≈ res_fr.t[:, :, i_fr] atol = 1e-8
    end
    for name in ("b", "c")
        i_rd = findfirst(==(name), res_rd.coef_names)
        i_fr = findfirst(==(name), res_fr.coef_names)
        @test res_rd.max_t_null[:, i_rd] ≈ res_fr.max_t_null[:, i_fr] rtol = 1e-5
        @test res_rd.p_corrected[:, :, i_rd] == res_fr.p_corrected[:, :, i_fr]
    end
end

