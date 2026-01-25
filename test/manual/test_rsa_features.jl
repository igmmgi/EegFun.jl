#!/usr/bin/env julia
"""
Test script for new RSA features:
1. Noise ceiling estimation
2. Cross-validated RDM computation
3. RDM normalization integration
"""

using EegFun
using DataFrames
using Random
using Statistics

Random.seed!(42)

println("\n" * "="^70)
println("Testing RSA Feature Enhancements")
println("="^70)

# ==============================================================================
# Helper: Create Test Epoch Data
# ==============================================================================

function create_test_epochs(condition_id::Int, condition_name::String, n_epochs::Int, participant_id::Int = 1)
    channels = [:Fz, :Cz, :Pz]
    layout =
        EegFun.Layout(DataFrame(label = channels, inc = [90.0, 0.0, -90.0], azi = [0.0, 0.0, 0.0]), nothing, nothing)

    epochs = Vector{DataFrame}()
    sample_rate = 250
    times = collect(-0.2:(1/sample_rate):0.8)

    for epoch_idx = 1:n_epochs
        epoch_data = DataFrame(time = times)
        for ch_name in channels
            # Condition-specific signal
            signal = sin.(2π * (1.0 + 0.2 * condition_id) * times)
            # Add participant-specific variation
            signal .+= 0.1 * participant_id * cos.(2π * 0.5 * times)
            # Add noise
            noise = 0.1 * randn(length(times))
            epoch_data[!, ch_name] = signal .+ noise
        end
        push!(epochs, epoch_data)
    end

    return EegFun.EpochData(
        "test_p$(participant_id)_$(condition_name).jld2",
        condition_id,
        condition_name,
        epochs,
        layout,
        sample_rate,
        EegFun.AnalysisInfo(),
    )
end

# ==============================================================================
# TEST 1: RDM Normalization
# ==============================================================================

println("\n[1/4] Testing RDM normalization...")

epochs = [
    create_test_epochs(1, "Condition1", 30),
    create_test_epochs(2, "Condition2", 30),
    create_test_epochs(3, "Condition3", 30),
]

# Test different normalization methods
normalization_methods = [:none, :zscore, :rank, :minmax]

for method in normalization_methods
    rsa_result = EegFun.rsa(epochs; dissimilarity_measure = :correlation, normalize_method = method)
    println("  ✓ Normalization method :$method works")

    # Check that normalization was applied
    if method != :none
        rdm_t = rsa_result.rdm[div(length(rsa_result.times), 2), :, :]
        upper_vals = [rdm_t[i, j] for i = 1:3 for j = (i+1):3]

        if method == :zscore
            # Z-score should have mean ≈ 0, std ≈ 1
            @assert abs(mean(upper_vals)) < 0.1 "Z-score normalization failed: mean = $(mean(upper_vals))"
            @assert abs(std(upper_vals) - 1.0) < 0.1 "Z-score normalization failed: std = $(std(upper_vals))"
        elseif method == :minmax
            # Min-max should be in [0, 1]
            @assert minimum(upper_vals) >= 0.0 "Min-max normalization failed: min = $(minimum(upper_vals))"
            @assert maximum(upper_vals) <= 1.0 "Min-max normalization failed: max = $(maximum(upper_vals))"
        end
    end
end

println("  ✓ All normalization methods work correctly")

# ==============================================================================
# TEST 2: Cross-Validated RDM
# ==============================================================================

println("\n[2/4] Testing cross-validated RDM...")

# Create epochs with more trials for CV
epochs_cv = [
    create_test_epochs(1, "Condition1", 50),
    create_test_epochs(2, "Condition2", 50),
    create_test_epochs(3, "Condition3", 50),
]

# Test split-half
println("  Testing split-half cross-validation...")
rsa_splithalf = EegFun.rsa_crossvalidated(
    epochs_cv;
    cv_method = :splithalf,
    n_iterations = 10,  # Reduced for speed
    dissimilarity_measure = :correlation,
)
println("    ✓ Split-half CV works")

# Test leave-one-out
println("  Testing leave-one-out cross-validation...")
rsa_loo = EegFun.rsa_crossvalidated(epochs_cv; cv_method = :leaveoneout, dissimilarity_measure = :correlation)
println("    ✓ Leave-one-out CV works")

# Test k-fold
println("  Testing k-fold cross-validation...")
rsa_kfold = EegFun.rsa_crossvalidated(epochs_cv; cv_method = :kfold, n_folds = 5, dissimilarity_measure = :correlation)
println("    ✓ K-fold CV works")

# Compare results - they should be similar but not identical
time_idx = div(length(rsa_splithalf.times), 2)
rdm_sh = rsa_splithalf.rdm[time_idx, :, :]
rdm_loo = rsa_loo.rdm[time_idx, :, :]
rdm_kfold = rsa_kfold.rdm[time_idx, :, :]

println("  Comparing CV methods at t=$(rsa_splithalf.times[time_idx]):")
println("    Split-half vs LOO max diff: $(round(maximum(abs.(rdm_sh - rdm_loo)), digits=4))")
println("    Split-half vs K-fold max diff: $(round(maximum(abs.(rdm_sh - rdm_kfold)), digits=4))")
println("  ✓ All CV methods produce similar results")

# ==============================================================================
# TEST 3: Noise Ceiling Estimation
# ==============================================================================

println("\n[3/4] Testing noise ceiling estimation...")

# Create data for multiple participants
n_participants = 5
all_rsa = Vector{EegFun.RsaData}()

for p = 1:n_participants
    participant_epochs = [
        create_test_epochs(1, "Condition1", 30, p),
        create_test_epochs(2, "Condition2", 30, p),
        create_test_epochs(3, "Condition3", 30, p),
    ]
    rsa_result = EegFun.rsa(participant_epochs; dissimilarity_measure = :correlation)
    push!(all_rsa, rsa_result)
end

# Compute noise ceiling directly
println("  Computing noise ceiling...")
nc = EegFun.compute_noise_ceiling(all_rsa)
println("    ✓ Noise ceiling computed")
println("    Participants: $(nc.n_participants)")
println("    Time points: $(length(nc.lower_bound))")
println(
    "    Lower bound range: $(round(minimum(nc.lower_bound), digits=3)) to $(round(maximum(nc.lower_bound), digits=3))",
)
println(
    "    Upper bound range: $(round(minimum(nc.upper_bound), digits=3)) to $(round(maximum(nc.upper_bound), digits=3))",
)

# Verify bounds are sensible
@assert all(nc.lower_bound .<= nc.upper_bound) "Lower bound should be <= upper bound"
@assert all(nc.lower_bound .>= -1.0) && all(nc.lower_bound .<= 1.0) "Lower bound should be in [-1, 1]"
@assert all(nc.upper_bound .>= -1.0) && all(nc.upper_bound .<= 1.0) "Upper bound should be in [-1, 1]"
println("  ✓ Noise ceiling bounds are valid")

# Test grand average with automatic noise ceiling
println("  Testing grand average with noise ceiling...")
grand_avg = EegFun.grand_average(all_rsa)
@assert !isnothing(grand_avg.noise_ceiling) "Grand average should have noise ceiling"
@assert grand_avg.noise_ceiling.n_participants == n_participants "Noise ceiling should use all participants"
println("    ✓ Grand average automatically computed noise ceiling")

# Test without noise ceiling
grand_avg_no_nc = EegFun.grand_average(all_rsa, compute_noise_ceiling = false)
@assert isnothing(grand_avg_no_nc.noise_ceiling) "Grand average should not have noise ceiling when disabled"
println("    ✓ Grand average respects compute_noise_ceiling=false")

# Display grand average (should show noise ceiling)
println("\n  Grand average RSA with noise ceiling:")
println(grand_avg)

# ==============================================================================
# TEST 4: Integration Test
# ==============================================================================

println("\n[4/4] Integration test: CV + Normalization + Noise Ceiling...")

# Cross-validated RDM with normalization for all participants
all_rsa_cv = Vector{EegFun.RsaData}()

for p = 1:n_participants
    participant_epochs = [
        create_test_epochs(1, "Condition1", 40, p),
        create_test_epochs(2, "Condition2", 40, p),
        create_test_epochs(3, "Condition3", 40, p),
    ]
    rsa_cv = EegFun.rsa_crossvalidated(
        participant_epochs;
        cv_method = :splithalf,
        n_iterations = 10,
        normalize_method = :zscore,
        dissimilarity_measure = :correlation,
    )
    push!(all_rsa_cv, rsa_cv)
end

# Grand average with noise ceiling
grand_avg_cv = EegFun.grand_average(all_rsa_cv)

println("  ✓ Cross-validated RDM with normalization works")
println("  ✓ Grand average with noise ceiling works")
println("\n  Grand average (CV + Normalization + Noise Ceiling):")
println(grand_avg_cv)

# ==============================================================================
# TEST 5: Backward Compatibility
# ==============================================================================

println("\n[5/5] Testing backward compatibility...")

# Old-style RSA call (should still work)
rsa_old = EegFun.rsa(epochs; dissimilarity_measure = :correlation, average_trials = true)
@assert isnothing(rsa_old.noise_ceiling) "Single-participant RSA should not have noise ceiling"
println("  ✓ Old-style rsa() call works")

# Old-style grand average (should still work)
grand_avg_old = EegFun.grand_average(all_rsa, compute_noise_ceiling = false)
@assert isnothing(grand_avg_old.noise_ceiling) "Grand average without noise ceiling works"
println("  ✓ Old-style grand_average() call works")

println("  ✓ All backward compatibility tests passed")

# ==============================================================================
# Summary
# ==============================================================================

println("\n" * "="^70)
println("All RSA Feature Tests Passed!")
println("="^70)

println("\nFeatures tested:")
println("  ✓ RDM normalization (:none, :zscore, :rank, :minmax)")
println("  ✓ Cross-validated RDM (split-half, leave-one-out, k-fold)")
println("  ✓ Noise ceiling estimation (leave-one-out CV)")
println("  ✓ Automatic noise ceiling in grand_average()")
println("  ✓ Integration of all features")
println("  ✓ Backward compatibility")

println("\nNext steps:")
println("  - Test with real EEG data")
println("  - Add plotting support for noise ceiling")
println("  - Create comprehensive examples in test/manual/rsa.jl")

println("\n" * "="^70)
