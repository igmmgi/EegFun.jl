# Validation script for Mass-Univariate Mixed Models

using Test
using EegFun
using MixedModels
using StatsModels
using DataFrames
using BenchmarkTools

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
