#!/usr/bin/env julia
"""
Example: RSA with Multiple Model Types

This example demonstrates how to handle different types of "other" data
when performing RSA analysis.
"""

using eegfun
using DataFrames
using Random

Random.seed!(42)

println("="^70)
println("RSA WITH MULTIPLE MODEL TYPES - EXAMPLE")
println("="^70)

# ============================================================================
# SETUP: Create example neural data
# ============================================================================

println("\n[1/4] Creating example neural data...")

# Create example epoch data (simplified - in practice, use your real data)
channels = [:Fz, :Cz, :Pz]
layout = eegfun.Layout(
    DataFrame(label = channels, inc = [90.0, 0.0, -90.0], azi = [0.0, 0.0, 0.0]),
    nothing,
    nothing,
)

# Create 3 conditions: Face, Object, Scene
n_epochs = 20
sample_rate = 250
times = collect(-0.2:(1/sample_rate):0.8)

epochs = eegfun.EpochData[]
for (cond_idx, cond_name) in enumerate(["Face", "Object", "Scene"])
    epoch_dfs = DataFrame[]
    for epoch_idx in 1:n_epochs
        df = DataFrame(time = times)
        for ch in channels
            # Create condition-specific signal
            signal = sin.(2π * (0.5 + cond_idx * 0.1) * times) .+ 0.1 * randn(length(times))
            df[!, ch] = signal
        end
        push!(epoch_dfs, df)
    end
    push!(epochs, eegfun.EpochData(
        "test_data",
        cond_idx,
        cond_name,
        epoch_dfs,
        layout,
        sample_rate,
        eegfun.AnalysisInfo(),
    ))
end

println("  ✓ Created epoch data for 3 conditions: $(join([e.condition_name for e in epochs], ", "))")

# ============================================================================
# STEP 1: Compute Neural RDM
# ============================================================================

println("\n[2/4] Computing neural RDM from EEG data...")

neural_rsa = eegfun.rsa(
    epochs,
    channels;
    time_range = (-0.2, 0.8),
    dissimilarity_measure = :correlation,
)

println("  ✓ Neural RDM computed")
println("    Conditions: $(join(neural_rsa.condition_names, ", "))")
println("    Time points: $(length(neural_rsa.times))")

# ============================================================================
# STEP 2: Create Model RDMs from Different Data Types
# ============================================================================

println("\n[3/4] Creating model RDMs from different data types...")

# Example 1: Feature vectors (e.g., word embeddings, image features)
println("\n  Model 1: Semantic embeddings (feature vectors)")
word_embeddings = [
    [0.1, 0.2, 0.3, 0.4],  # Face
    [0.2, 0.3, 0.4, 0.5],  # Object
    [0.5, 0.6, 0.7, 0.8],  # Scene
]
semantic_rdm = eegfun.create_rdm_from_vectors(word_embeddings, dissimilarity_measure=:euclidean)
println("    ✓ Created RDM from feature vectors")

# Example 2: Reaction times (behavioral data)
println("\n  Model 2: Reaction times")
rts = [0.35, 0.42, 0.38]  # RTs for Face, Object, Scene
rt_rdm = eegfun.create_rdm_from_reaction_times(rts)
println("    ✓ Created RDM from reaction times")

# Example 3: Similarity ratings (behavioral data)
println("\n  Model 3: Similarity ratings")
similarity_matrix = [
    1.0 2.5 5.0;  # Face-Face, Face-Object, Face-Scene
    2.5 1.0 4.5;  # Object-Face, Object-Object, Object-Scene
    5.0 4.5 1.0;  # Scene-Face, Scene-Object, Scene-Scene
]
similarity_rdm = eegfun.create_rdm_from_similarity_ratings(similarity_matrix)
println("    ✓ Created RDM from similarity ratings")

# Example 4: Categorical labels (theoretical model)
println("\n  Model 4: Category labels")
categories = [1, 1, 2]  # Face and Object in category 1, Scene in category 2
category_rdm = eegfun.create_rdm_from_categorical(categories)
println("    ✓ Created RDM from category labels")

# Example 5: Data matrix (e.g., image features as matrix)
println("\n  Model 5: Image features (matrix)")
image_features = [
    0.1 0.2 0.3 0.4;  # Face
    0.2 0.3 0.4 0.5;  # Object
    0.5 0.6 0.7 0.8;  # Scene
]
image_rdm = eegfun.create_rdm_from_matrix(image_features, dissimilarity_measure=:correlation)
println("    ✓ Created RDM from data matrix")

# ============================================================================
# STEP 3: Compare Neural RDM to All Models
# ============================================================================

println("\n[4/4] Comparing neural RDM to all models...")

# Option A: Compare individually
println("\n  Option A: Compare models individually")
rsa_with_semantic = eegfun.compare_models(
    deepcopy(neural_rsa),
    [semantic_rdm],
    model_names=["Semantic Embeddings"],
    n_permutations=100,
)
println("    ✓ Compared to semantic embeddings")

# Option B: Compare all models at once (recommended)
println("\n  Option B: Compare all models at once")
all_model_rdms = [semantic_rdm, rt_rdm, similarity_rdm, category_rdm, image_rdm]
all_model_names = [
    "Semantic Embeddings",
    "Reaction Times",
    "Similarity Ratings",
    "Category Labels",
    "Image Features",
]

rsa_with_all_models = eegfun.compare_models(
    neural_rsa,
    all_model_rdms,
    model_names=all_model_names,
    correlation_type=:spearman,
    n_permutations=100,
)

println("    ✓ Compared to all $(length(all_model_names)) models")

# Option C: Use convenience function for mixed data types
println("\n  Option C: Use convenience function for mixed data types")
models_dict = Dict(
    "Semantic Embeddings" => word_embeddings,
    "Reaction Times" => rts,
    "Similarity Ratings" => similarity_matrix,
    "Category Labels" => categories,
    "Image Features" => image_features,
)

auto_rdms, auto_names = eegfun.create_model_rdms(models_dict)
rsa_with_auto = eegfun.compare_models(
    deepcopy(neural_rsa),
    auto_rdms,
    model_names=auto_names,
    n_permutations=100,
)
println("    ✓ Automatically converted and compared $(length(auto_names)) models")

# ============================================================================
# STEP 4: Visualize Results
# ============================================================================

println("\n[5/5] Visualizing results...")

# Plot model correlations
fig = eegfun.plot_model_correlations(
    rsa_with_all_models,
    title="Neural RDM vs Multiple Models",
    colors=[:blue, :red, :green, :orange, :purple],
)
println("  ✓ Model correlation plot created")

# Print summary
println("\n" * "="^70)
println("SUMMARY")
println("="^70)
println("\nNeural RDM computed from EEG data:")
println("  - Conditions: $(join(neural_rsa.condition_names, ", "))")
println("  - Time points: $(length(neural_rsa.times))")

println("\nModel RDMs created from:")
println("  1. Feature vectors (semantic embeddings)")
println("  2. Behavioral data (reaction times)")
println("  3. Behavioral data (similarity ratings)")
println("  4. Categorical labels")
println("  5. Data matrix (image features)")

println("\nAll models compared to neural RDM using Spearman correlation.")
println("Check the plot to see which models best match neural patterns over time.")

println("\n" * "="^70)

