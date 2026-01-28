"""
Manual test for RSA (Representational Similarity Analysis) functionality

This script tests the RSA workflow:
1. Create test epoch data for multiple conditions
2. Compute RDMs for each participant separately
3. Compare to model RDMs
4. Average across participants (grand average)
5. Plot results

Tests both basic RSA and model comparison.
"""

using EegFun
using DataFrames
using Random
using Glob

# Set seed for reproducibility
Random.seed!(42)

println("\n" * "="^70)
println("RSA (REPRESENTATIONAL SIMILARITY ANALYSIS) MANUAL TEST")
println("="^70)

# ============================================================================
# SETUP: Create Test Data
# ============================================================================

println("\n[1/5] Creating test epoch data...")

# Create test epoch data for multiple participants and conditions
# Each participant will have epochs for 3 conditions
n_participants = 3
n_epochs_per_condition = 30
channels = [:Fz, :Cz, :Pz, :POz]

# Helper function to create epoch data for a participant and condition
function create_participant_condition_epochs(
    participant_id::Int,
    condition_id::Int,
    condition_name::String,
    n_epochs::Int,
    base_freq::Float64,
)
    # Create layout
    layout = EegFun.Layout(DataFrame(label = channels, inc = [90.0, 0.0, -90.0, -90.0], azi = [0.0, 0.0, 0.0, 0.0]), nothing, nothing)

    # Create epochs with condition-specific patterns
    epochs = Vector{DataFrame}()
    sample_rate = 250
    time_window = (-0.2, 0.8)
    times = collect(time_window[1]:(1/sample_rate):time_window[2])

    for epoch_idx = 1:n_epochs
        epoch_data = DataFrame(time = times)

        # Add condition-specific signal patterns
        for (ch_idx, ch_name) in enumerate(channels)
            # Base signal with condition-specific frequency
            signal = sin.(2π * base_freq * times)

            # Add condition-specific phase shift
            phase_shift = condition_id * π / 4
            signal = sin.(2π * base_freq * times .+ phase_shift)

            # Add channel-specific amplitude
            amplitude = 1.0 + 0.3 * ch_idx

            # Add noise
            noise = 0.2 * randn(length(times))

            # Add participant-specific offset
            participant_offset = 0.1 * participant_id

            epoch_data[!, ch_name] = amplitude .* signal .+ noise .+ participant_offset
        end

        push!(epochs, epoch_data)
    end

    # Create EpochData object
    epoch_data_obj = EegFun.EpochData(
        "test_p$(participant_id)_$(condition_name).jld2",
        condition_id,
        condition_name,
        epochs,
        layout,
        sample_rate,
        EegFun.AnalysisInfo(),
    )

    return epoch_data_obj
end

# Create epochs for all participants and conditions
all_participant_epochs = Vector{Vector{EegFun.EpochData}}()

for p = 1:n_participants
    participant_epochs = [
        create_participant_condition_epochs(p, 1, "Face", n_epochs_per_condition, 2.0),
        create_participant_condition_epochs(p, 2, "Object", n_epochs_per_condition, 1.8),
        create_participant_condition_epochs(p, 3, "Scene", n_epochs_per_condition, 1.6),
    ]
    push!(all_participant_epochs, participant_epochs)
end

println("  ✓ Created epoch data for $n_participants participants")
println("  ✓ Each participant has $n_epochs_per_condition epochs per condition")
println("  ✓ Conditions: Face, Object, Scene")

# ============================================================================
# TEST 1: Basic RSA Analysis
# ============================================================================

println("\n[2/5] Test 1: Basic RSA analysis with correlation-based dissimilarity...")

# Perform RSA for each participant
all_rsa_results = Vector{EegFun.RsaData}()

for (p_idx, participant_epochs) in enumerate(all_participant_epochs)
    rsa_result = EegFun.rsa(
        participant_epochs;
        channel_selection = EegFun.channels(channels),
        sample_selection = EegFun.samples((-0.2, 0.8)),
        dissimilarity_measure = :correlation,
        average_trials = true,
    )
    push!(all_rsa_results, rsa_result)
    println("  ✓ Participant $p_idx RSA complete")
end

# Compute grand average
println("\n[3/5] Computing grand average RSA...")
grand_avg_rsa = EegFun.grand_average(all_rsa_results)
println("  ✓ Grand average created")
println("    Conditions: $(join(grand_avg_rsa.condition_names, ", "))")
println("    Time range: $(round(grand_avg_rsa.times[1], digits=2)) to $(round(grand_avg_rsa.times[end], digits=2)) s")

# Plot grand average RDM at a specific time point
println("\n[4/5] Plotting RSA results...")
fig1 = EegFun.plot_rdm(grand_avg_rsa, time_point = 0.3, title = "Grand Average RDM at 0.3s")
println("  ✓ RDM plot created")

# Plot average RDM
fig2 = EegFun.plot_rdm(grand_avg_rsa, title = "Grand Average RDM (averaged across time)")
println("  ✓ Average RDM plot created")

# ============================================================================
# TEST 2: Model Comparison
# ============================================================================

println("\n[5/5] Test 2: Model comparison...")

if isnothing(grand_avg_rsa)
    println("  ⚠ Skipping model comparison - grand average RSA not available")
else
    try
        # Create example model RDMs
        # Model 1: Face and Object are similar, Scene is different
        model1_rdm = [
            0.0 0.3 0.8
            0.3 0.0 0.7
            0.8 0.7 0.0
        ]

        # Model 2: All conditions equally dissimilar
        model2_rdm = [
            0.0 0.5 0.5
            0.5 0.0 0.5
            0.5 0.5 0.0
        ]

        # Compare grand average RSA to models
        grand_avg_with_models = EegFun.compare_models(
            grand_avg_rsa,
            [model1_rdm, model2_rdm];
            model_names = ["Semantic Model", "Equal Distance Model"],
            correlation_type = :spearman,
            n_permutations = 100,  # Reduced for testing
        )

        println("  ✓ Model comparison complete")
        println("    Models compared: $(join(grand_avg_with_models.model_names, ", "))")

        # Plot model correlations
        fig3 = EegFun.plot_model_correlations(grand_avg_with_models, title = "Model Correlations Over Time", colors = [:red, :blue])
        println("  ✓ Model correlation plot created")

    catch e
        println("  ✗ Error: $e")
        rethrow(e)
    end
end

# ============================================================================
# TEST 3: Different Dissimilarity Measures
# ============================================================================

println("\n[6/7] Test 3: Different dissimilarity measures...")

# Test different dissimilarity measures on first participant
participant_epochs = all_participant_epochs[1]

dissimilarity_measures = [:correlation, :spearman, :euclidean]
rsa_results_different_measures = Dict{Symbol,EegFun.RsaData}()

for measure in dissimilarity_measures
    rsa_result = EegFun.rsa(
        participant_epochs;
        channel_selection = EegFun.channels(channels),
        sample_selection = EegFun.samples((-0.2, 0.8)),
        dissimilarity_measure = measure,
        average_trials = true,
    )
    rsa_results_different_measures[measure] = rsa_result
    println("  ✓ RSA with $measure dissimilarity complete")

    # Plot RDM at a specific time point for comparison
    fig = EegFun.plot_rdm(rsa_result, time_point = 0.3, title = "RDM at 0.3s using $measure dissimilarity")
end

println("  ✓ Compared different dissimilarity measures")
println("    Note: Different measures may give different RDM structures")

# ============================================================================
# TEST 4: Creating Model RDMs from Different Data Types
# ============================================================================

println("\n[7/7] Test 4: Creating model RDMs from different data types...")

# Get neural RSA for model comparison
neural_rsa = all_rsa_results[1]

# 4a. Create RDM from feature vectors (e.g., word embeddings)
println("\n  4a. Creating RDM from feature vectors...")
word_embeddings = [
    [0.1, 0.2, 0.3, 0.4],  # Face embedding
    [0.2, 0.3, 0.4, 0.5],  # Object embedding
    [0.5, 0.6, 0.7, 0.8],  # Scene embedding
]
model_rdm_vectors = EegFun.create_rdm_from_vectors(word_embeddings, dissimilarity_measure = :euclidean)
println("    ✓ Created RDM from word embeddings")

# 4b. Create RDM from reaction times
println("\n  4b. Creating RDM from reaction times...")
rts = [0.3, 0.5, 0.4]  # RTs for Face, Object, Scene
model_rdm_rts = EegFun.create_rdm_from_reaction_times(rts)
println("    ✓ Created RDM from reaction times")
println("      RTs: Face=$(rts[1])s, Object=$(rts[2])s, Scene=$(rts[3])s")

# 4c. Create RDM from similarity ratings
println("\n  4c. Creating RDM from similarity ratings...")
similarity_matrix = [
    1.0  2.0  5.0  # Face-Object=2, Face-Scene=5 (1=very similar, 7=very different)
    2.0  1.0  4.0  # Object-Face=2, Object-Scene=4
    5.0  4.0  1.0   # Scene-Face=5, Scene-Object=4
]
model_rdm_similarity = EegFun.create_rdm_from_similarity_ratings(similarity_matrix, convert_to_dissimilarity = true)
println("    ✓ Created RDM from similarity ratings")

# 4d. Create RDM from categorical labels
println("\n  4d. Creating RDM from categorical labels...")
categories = [1, 1, 2]  # Face=1, Object=1 (same category), Scene=2 (different)
model_rdm_categorical = EegFun.create_rdm_from_categorical(categories)
println("    ✓ Created RDM from categorical labels")
println("      Categories: Face=1, Object=1, Scene=2")

# 4e. Create RDM from data matrix
println("\n  4e. Creating RDM from data matrix...")
data_matrix = [
    0.1  0.2  0.3  0.4  # Face features
    0.2  0.3  0.4  0.5  # Object features
    0.5  0.6  0.7  0.8   # Scene features
]
model_rdm_matrix = EegFun.create_rdm_from_matrix(data_matrix, dissimilarity_measure = :correlation)
println("    ✓ Created RDM from data matrix")

# Compare neural RDM to all these models
println("\n  Comparing neural RDM to all model types...")
all_model_rdms = [model_rdm_vectors, model_rdm_rts, model_rdm_similarity, model_rdm_categorical, model_rdm_matrix]
model_names = ["Word Embeddings", "Reaction Times", "Similarity Ratings", "Categorical", "Feature Matrix"]

rsa_with_all_models = EegFun.compare_models(
    neural_rsa,
    all_model_rdms;
    model_names = model_names,
    correlation_type = :spearman,
    n_permutations = 100,  # Reduced for testing
)

println("    ✓ Compared neural RDM to $(length(model_names)) different model types")

# Plot model correlations
fig_models = EegFun.plot_model_correlations(
    rsa_with_all_models,
    title = "Neural RDM vs Multiple Model Types",
    colors = [:blue, :red, :green, :orange, :purple],
)
println("    ✓ Model correlation plot created")

# ============================================================================
# TEST 5: Temporal Model RDMs
# ============================================================================

println("\n[8/8] Test 5: Temporal model RDMs...")

# Create temporal model data (e.g., eye tracking or EDA)
# Format: [conditions × features × time]
n_conditions = 3
n_features = 2  # e.g., x and y position for eye tracking
n_timepoints = length(neural_rsa.times)

# Simulate temporal model data (e.g., eye tracking)
temporal_model_data = zeros(Float64, n_conditions, n_features, n_timepoints)
temporal_times = neural_rsa.times

for cond_idx = 1:n_conditions
    for feat_idx = 1:n_features
        # Create condition-specific temporal pattern
        base_freq = 1.0 + 0.2 * cond_idx
        signal = sin.(2π * base_freq * temporal_times)
        temporal_model_data[cond_idx, feat_idx, :] = signal .+ 0.1 * randn(n_timepoints)
    end
end

println("  Created temporal model data (simulated eye tracking)")
println("    Shape: [conditions × features × time] = [$n_conditions × $n_features × $n_timepoints]")

# Create temporal RDMs
temporal_rdms = EegFun.create_temporal_rdm(temporal_model_data, temporal_times; dissimilarity_measure = :correlation)

println("  ✓ Created temporal RDMs")
println("    Shape: [time × condition × condition] = [$(size(temporal_rdms, 1)) × $(size(temporal_rdms, 2)) × $(size(temporal_rdms, 3))]")

# Compare neural RDM to temporal model
rsa_with_temporal = EegFun.compare_models(
    neural_rsa,
    [temporal_rdms];
    model_names = ["Eye Tracking Model"],
    correlation_type = :spearman,
    n_permutations = 100,
)

println("  ✓ Compared neural RDM to temporal model")

# Plot temporal model correlations
fig_temporal = EegFun.plot_model_correlations(rsa_with_temporal, title = "Neural RDM vs Temporal Model (Eye Tracking)", colors = [:blue])
println("  ✓ Temporal model correlation plot created")

# ============================================================================
# TEST 6: RDM Visualization Options
# ============================================================================

println("\n[9/9] Test 6: RDM visualization options...")

# Plot RDM at different time points
time_points_to_plot = [-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
println("  Plotting RDMs at multiple time points: $time_points_to_plot")

for t in time_points_to_plot
    if t >= grand_avg_rsa.times[1] && t <= grand_avg_rsa.times[end]
        fig = EegFun.plot_rdm(grand_avg_rsa, time_point = t, title = "RDM at $(round(t, digits=2))s", colormap = :viridis)
    end
end

# Plot average RDM with different colormaps
colormaps = [:viridis, :plasma, :inferno, :magma]
println("  Plotting average RDM with different colormaps...")

for cmap in colormaps
    fig = EegFun.plot_rdm(grand_avg_rsa, title = "Average RDM (colormap: $cmap)", colormap = cmap)
end

println("  ✓ RDM visualization examples complete")

# ============================================================================
# TEST 7: Real EEG Data
# ============================================================================

println("\n" * "="^70)
println("TEST 7: Real EEG Data")
println("="^70)

# Check if data directory exists (optional - skip if not available)
data_dir = joinpath(@__DIR__, "..", "..", "..", "AttentionExp", "recoded")
layout_file = EegFun.read_layout("./data/layouts/biosemi/biosemi72.csv")
EegFun.polar_to_cartesian_xy!(layout_file)

# Epoch configuration - need at least 2 conditions for RSA (3+ is better)
epoch_cfg = [
    EegFun.EpochCondition(name = "Condition1", trigger_sequences = [[1]]),
    EegFun.EpochCondition(name = "Condition2", trigger_sequences = [[2]]),
    EegFun.EpochCondition(name = "Condition3", trigger_sequences = [[3]]),
    EegFun.EpochCondition(name = "Condition4", trigger_sequences = [[4]]),
]

# RSA parameters
dissimilarity_measure = :correlation  # or :spearman, :euclidean
average_trials = true
time_range = (-0.2, 0.8)  # Adjust based on your epoch window

using Glob

# Find all .bdf files in directory
bdf_files = glob("*.bdf", data_dir)

if isempty(bdf_files)
    println("\n⚠ No BDF files found in $data_dir")
    println("  Skipping real data test. To use real data:")
    println("  1. Place .bdf files in the data directory")
    println("  2. Update data_dir path if needed")
    println("  3. Update epoch_cfg with your trigger sequences")
else
    println("\n[1/4] Found $(length(bdf_files)) BDF file(s) to process")

    # Store RSA results for each participant
    all_rsa_results = EegFun.RsaData[]

    # Process each file
    for (file_idx, data_file) in enumerate(bdf_files)
        participant_id = basename(data_file)
        println("\n[$(file_idx)/$(length(bdf_files))] Processing: $participant_id")

        try
            # Preprocessing (same for all files)
            println("  Loading and preprocessing...")
            dat = EegFun.read_raw_data(data_file)
            dat = EegFun.create_eeg_dataframe(dat, layout_file)
            EegFun.rereference!(dat, :avg)
            EegFun.highpass_filter!(dat, 1)

            # Extract epochs
            println("  Extracting epochs...")
            epochs = EegFun.extract_epochs(dat, epoch_cfg, -0.5, 2.0)

            # Check we have enough conditions for RSA
            if length(epochs) < 2
                println("  ⚠ Skipping: Only $(length(epochs)) condition(s) found (need at least 2)")
                continue
            end

            # Perform RSA
            println("  Computing RSA...")
            rsa_result = EegFun.rsa(epochs; dissimilarity_measure = :correlation, average_trials = true)

            push!(all_rsa_results, rsa_result)
            println("  ✓ RSA complete")
            println("    Conditions: $(join(rsa_result.condition_names, ", "))")
            println("    Time range: $(round(rsa_result.times[1], digits=2)) to $(round(rsa_result.times[end], digits=2)) s")
            println("    Channels: $(length(rsa_result.channels))")

        catch e
            println("  ✗ Error processing $participant_id: $e")
            if isa(e, InterruptException)
                rethrow(e)
            end
            # Continue with next file
            continue
        end
    end

    # Create grand average
    if isempty(all_rsa_results)
        println("\n⚠ No participants successfully processed!")
        println("  Check that:")
        println("    - BDF files are in the correct directory")
        println("    - Epoch configuration matches your trigger sequences")
        println("    - Channels exist in your data")
    else
        println("\n[2/4] Creating grand average RSA...")
        grand_avg_rsa_real = EegFun.grand_average(all_rsa_results)

        println("  ✓ Grand average created")
        println("    Participants: $(length(all_rsa_results))")
        println("    Conditions: $(join(grand_avg_rsa_real.condition_names, ", "))")
        println("    Time range: $(round(grand_avg_rsa_real.times[1], digits=2)) to $(round(grand_avg_rsa_real.times[end], digits=2)) s")

        # Plot grand average RDM at key time points
        println("\n[3/4] Plotting grand average RDMs...")

        # Plot average RDM across all time
        fig1 = EegFun.plot_rdm(grand_avg_rsa_real, title = "Grand Average RDM (Real Data, Averaged Across Time)", colormap = :viridis)
        println("  ✓ Average RDM plot created")

        # Plot RDM at stimulus onset (0.0s)
        if 0.0 >= grand_avg_rsa_real.times[1] && 0.0 <= grand_avg_rsa_real.times[end]
            fig2 = EegFun.plot_rdm(
                grand_avg_rsa_real,
                time_point = 0.0,
                title = "Grand Average RDM at Stimulus Onset (0.0s)",
                colormap = :viridis,
            )
            println("  ✓ RDM at 0.0s plot created")
        end

        # Plot RDM at a later time point (e.g., 0.3s)
        if 0.3 >= grand_avg_rsa_real.times[1] && 0.3 <= grand_avg_rsa_real.times[end]
            fig3 = EegFun.plot_rdm(grand_avg_rsa_real, time_point = 0.3, title = "Grand Average RDM at 0.3s", colormap = :viridis)
            println("  ✓ RDM at 0.3s plot created")
        end

        # Model comparison example (if you have model data)
        println("\n[4/4] Model comparison example...")

        # Example: Create a simple model RDM based on condition similarity
        # In real analysis, you would use actual behavioral/computational models
        if length(grand_avg_rsa_real.condition_names) >= 2
            # Example model: Assume conditions 1 and 2 are more similar
            n_conds = length(grand_avg_rsa_real.condition_names)
            example_model_rdm = zeros(Float64, n_conds, n_conds)

            # Create a simple model where first two conditions are similar
            for i = 1:n_conds
                for j = 1:n_conds
                    if i == j
                        example_model_rdm[i, j] = 0.0
                    elseif (i == 1 && j == 2) || (i == 2 && j == 1)
                        example_model_rdm[i, j] = 0.3  # Similar
                    else
                        example_model_rdm[i, j] = 0.8  # Different
                    end
                end
            end

            # Compare to model
            rsa_with_model = EegFun.compare_models(
                grand_avg_rsa_real,
                [example_model_rdm];
                model_names = ["Example Model (Conditions 1&2 Similar)"],
                correlation_type = :spearman,
                n_permutations = 100,  # Reduced for testing
            )

            println("  ✓ Model comparison complete")

            # Plot model correlations
            fig4 = EegFun.plot_model_correlations(rsa_with_model, title = "Neural RDM vs Example Model (Real Data)", colors = [:blue])
            println("  ✓ Model correlation plot created")

            # Find time point with maximum correlation
            max_corr_idx = argmax(rsa_with_model.model_correlations[:, 1])
            max_corr_time = rsa_with_model.times[max_corr_idx]
            max_corr_val = rsa_with_model.model_correlations[max_corr_idx, 1]
            println("    Max correlation: $(round(max_corr_val, digits=3)) at $(round(max_corr_time, digits=3)) s")
        end

        println("\n" * "="^70)
        println("Real Data RSA Analysis Complete!")
        println("="^70)
        println("\nSummary:")
        println("  ✓ Processed $(length(all_rsa_results)) participant(s)")
        println("  ✓ Grand average RSA computed")
        println("  ✓ RDMs visualized at multiple time points")
        println("  ✓ Model comparison performed")
        println("\nNext steps:")
        println("  - Create model RDMs from your behavioral/computational data")
        println("  - Compare neural RDMs to your models using compare_models()")
        println("  - Explore different time windows and channel selections")
    end
end

# ============================================================================
# FINAL SUMMARY
# ============================================================================

println("\n" * "="^70)
println("ALL RSA TESTS COMPLETE")
println("="^70)

println("\nTests performed:")
println("  ✓ Basic RSA with correlation-based dissimilarity")
println("  ✓ Grand averaging across participants")
println("  ✓ Model comparison with permutation testing")
println("  ✓ Different dissimilarity measures (correlation, spearman, euclidean)")
println("  ✓ Creating model RDMs from:")
println("    - Feature vectors (word embeddings)")
println("    - Reaction times")
println("    - Similarity ratings")
println("    - Categorical labels")
println("    - Data matrices")
println("  ✓ Temporal model RDMs (eye tracking simulation)")
println("  ✓ RDM visualization at multiple time points")
println("  ✓ Different colormaps for RDM visualization")
# Check if real data was processed
real_data_processed = false
try
    if isdir(data_dir)
        try
            real_data_processed = !isempty(glob("*.bdf", data_dir))
        catch
            real_data_processed = !isempty(filter(f -> endswith(f, ".bdf"), readdir(data_dir)))
        end
    end
catch
end
if real_data_processed
    println("  ✓ Real EEG data analysis")
end
println("  ✓ RDM and model correlation visualization")

println("\nAll plots should be displayed above.")
println("Check that:")
println("  - RDMs show condition-specific patterns")
println("  - Different dissimilarity measures produce different RDM structures")
println("  - Model correlations vary over time")
println("  - Significant correlations are marked (if p-values computed)")
println("  - Temporal models show time-varying correlations")
if real_data_processed
    println("  - Real data RDMs show meaningful condition relationships")
end

println("\n" * "="^70)


