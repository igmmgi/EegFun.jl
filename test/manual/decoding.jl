#!/usr/bin/env julia
"""
Manual test for MVPA/Decoding functionality

This script tests the decoding workflow following erplab's approach:
1. Create test epoch data for multiple conditions
2. Decode for each participant separately
3. Average across participants (grand average)
4. Plot results

Tests both binary and multi-class classification with different models.
"""

using eegfun
using DataFrames
using Random

# Set seed for reproducibility
Random.seed!(42)

println("\n" * "="^70)
println("MVPA/DECODING MANUAL TEST")
println("="^70)

# ============================================================================
# INSTALL REQUIRED PACKAGES
# ============================================================================

println("\n[0/6] Installing required packages...")

try
    using Pkg
    packages_to_install = ["MLJ", "MLJLinearModels", "LIBSVM"]
    
    for pkg in packages_to_install
        try
            # Try to load it first
            eval(Meta.parse("using $pkg"))
            println("  ✓ $pkg already installed")
        catch
            # If not available, install it
            println("  Installing $pkg...")
            Pkg.add(pkg)
            println("  ✓ $pkg installed")
        end
    end
catch e
    println("  ⚠ Could not install packages automatically: $e")
    println("    Please install manually: using Pkg; Pkg.add([\"MLJ\", \"MLJLinearModels\", \"LIBSVM\"])")
end

# ============================================================================
# SETUP: Create Test Data
# ============================================================================

println("\n[1/6] Creating test epoch data...")

# Create test epoch data for multiple participants and conditions
# Each participant will have epochs for 2 conditions (binary classification)
n_participants = 3
n_epochs_per_condition = 30  # Enough for cross-validation
n_timepoints = 200  # -0.2 to 0.8s at 1000Hz
sample_rate = 1000
channels = [:Fz, :Cz, :Pz, :POz]  # 4 channels for MVPA

# Create layout (must match channel names exactly)
layout = eegfun.Layout(
    DataFrame(
        label = channels,  # [:Fz, :Cz, :Pz, :POz]
        inc = [90.0, 0.0, -90.0, -90.0],  # Rough positions
        azi = [0.0, 0.0, 0.0, 0.0],
    ),
    nothing,
    nothing,
)

analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)

# Create time vector
time = collect(range(-0.2, 0.8, length = n_timepoints))

# Function to create epochs for one participant and one condition
function create_participant_condition_epochs(
    participant_id::Int,
    condition_id::Int,
    condition_name::String,
    n_epochs::Int,
    signal_strength::Float64,
)
    epochs = DataFrame[]
    for epoch = 1:n_epochs
        df = DataFrame(time = time, epoch = fill(epoch, n_timepoints))

        # Create condition-specific signal patterns
        # Condition 1: Positive deflection around 200ms
        # Condition 2: Negative deflection around 200ms
        for (ch_idx, ch) in enumerate(channels)
            # Base signal with condition-specific pattern
            if condition_id == 1
                # Condition 1: Positive ERP component
                signal = signal_strength * exp.(-((time .- 0.2).^2) / (2 * 0.05^2)) .+
                        0.3 * signal_strength * exp.(-((time .- 0.4).^2) / (2 * 0.08^2))
            else
                # Condition 2: Negative ERP component
                signal = -signal_strength * exp.(-((time .- 0.2).^2) / (2 * 0.05^2)) .+
                        0.2 * signal_strength * exp.(-((time .- 0.5).^2) / (2 * 0.1^2))
            end

            # Add channel-specific variation and noise
            channel_factor = 1.0 + (ch_idx - 1) * 0.1
            noise = 0.3 * randn(n_timepoints)  # SNR ~ 3:1
            df[!, ch] = signal .* channel_factor .+ noise
        end

        push!(epochs, df)
    end

    return eegfun.EpochData(
        "participant_$(participant_id)",
        condition_id,
        condition_name,
        epochs,
        layout,
        sample_rate,
        analysis_info,
    )
end

# Create epochs for all participants
all_participant_epochs = []
for p in 1:n_participants
    # Participant p, Condition 1
    cond1 = create_participant_condition_epochs(p, 1, "Face", n_epochs_per_condition, 2.0)
    # Participant p, Condition 2
    cond2 = create_participant_condition_epochs(p, 2, "Object", n_epochs_per_condition, 2.0)
    push!(all_participant_epochs, [cond1, cond2])
end

println("  ✓ Created epoch data for $n_participants participants")
println("  ✓ Each participant has $(n_epochs_per_condition) epochs per condition")
println("  ✓ Conditions: Face vs Object")

# ============================================================================
# TEST 1: Binary Classification with Default LDA
# ============================================================================

println("\n[2/6] Test 1: Binary classification with default LDA model...")

try
    # MLJ is already loaded at eegfun module level
    # Load LogisticClassifier (linear classifier, similar to LDA) for default
    LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
    lda_model = LogisticClassifier()

    # Decode for Participant 1
    epochs_p1 = all_participant_epochs[1]
    decoded_p1 = eegfun.decode(
        epochs_p1,
        channels;
        model = lda_model,  # Use explicit model
        time_range = (-0.2, 0.8),
        n_iterations = 10,  # Reduced for testing
        n_folds = 3,
        equalize_trials = true,
    )

    println("  ✓ Participant 1 decoding complete")
    println("    Max accuracy: $(round(maximum(decoded_p1.average_score), digits=3))")
    println("    Chance level: $(round(decoded_p1.chance_level, digits=3))")

    # Decode for Participant 2
    epochs_p2 = all_participant_epochs[2]
    decoded_p2 = eegfun.decode(
        epochs_p2,
        channels;
        model = lda_model,
        time_range = (-0.2, 0.8),
        n_iterations = 10,
        n_folds = 3,
        equalize_trials = true,
    )

    println("  ✓ Participant 2 decoding complete")

    # Decode for Participant 3
    epochs_p3 = all_participant_epochs[3]
    decoded_p3 = eegfun.decode(
        epochs_p3,
        channels;
        model = lda_model,
        time_range = (-0.2, 0.8),
        n_iterations = 10,
        n_folds = 3,
        equalize_trials = true,
    )

    println("  ✓ Participant 3 decoding complete")

    # Grand average
    println("\n[3/6] Creating grand average across participants...")
    all_decoded = [decoded_p1, decoded_p2, decoded_p3]
    grand_avg = eegfun.grand_average(all_decoded)

    println("  ✓ Grand average created")
    println("    Max accuracy: $(round(maximum(grand_avg.average_score), digits=3))")
    println("    Time range: $(round(grand_avg.times[1], digits=2)) to $(round(grand_avg.times[end], digits=2)) s")

    # Plot grand average
    println("\n[4/6] Plotting grand average decoding results...")
    fig1 = eegfun.plot_decoding(grand_avg, title = "Grand Average: Face vs Object (LDA)")
    println("  ✓ Plot created")

catch e
    println("  ✗ Error: $e")
    println("    Make sure MLJ.jl is installed: using Pkg; Pkg.add(\"MLJ\")")
    rethrow(e)
end

# ============================================================================
# TEST 2: Binary Classification with SVM (ERPLAB-like)
# ============================================================================

println("\n[5/6] Test 2: Binary classification with SVM (ERPLAB-like)...")

try
    using MLJ
    using LIBSVM

    # Load SVM classifier
    SVMClassifier = @load SVMClassifier pkg=LIBSVM

    # Create SVM model with linear kernel (erplab default)
    svm_model = SVMClassifier(kernel = LIBSVM.Kernel.Linear)

    # Decode for Participant 1 with SVM
    epochs_p1 = all_participant_epochs[1]
    decoded_p1_svm = eegfun.decode(
        epochs_p1,
        channels;
        model = svm_model,
        time_range = (-0.2, 0.8),
        n_iterations = 10,  # Reduced for testing
        n_folds = 3,
        equalize_trials = true,
    )

    println("  ✓ SVM decoding complete")
    println("    Max accuracy: $(round(maximum(decoded_p1_svm.average_score), digits=3))")
    println("    Method: $(decoded_p1_svm.method)")

    # Plot SVM results
    fig2 = eegfun.plot_decoding(decoded_p1_svm, title = "Participant 1: Face vs Object (SVM)")

catch e
    if occursin("LIBSVM", string(e))
        println("  ⚠ LIBSVM.jl not available, skipping SVM test")
        println("    Install with: using Pkg; Pkg.add(\"LIBSVM\")")
    else
        println("  ✗ Error: $e")
        rethrow(e)
    end
end

# ============================================================================
# TEST 3: Multi-class Classification (3 conditions)
# ============================================================================

println("\n[6/6] Test 3: Multi-class classification (3 conditions)...")

try
    using MLJ

    # Create 3 conditions for one participant
    cond1 = create_participant_condition_epochs(1, 1, "Face", n_epochs_per_condition, 2.0)
    cond2 = create_participant_condition_epochs(1, 2, "Object", n_epochs_per_condition, 2.0)
    cond3 = create_participant_condition_epochs(1, 3, "Scene", n_epochs_per_condition, 1.8)

    epochs_3way = [cond1, cond2, cond3]

    # Load model for 3-way classification
    LogisticClassifier = @load LogisticClassifier pkg=MLJLinearModels
    lda_model_3way = LogisticClassifier()
    
    # Decode with 3 conditions
    decoded_3way = eegfun.decode(
        epochs_3way,
        channels;
        model = lda_model_3way,
        time_range = (-0.2, 0.8),
        n_iterations = 10,
        n_folds = 3,
        equalize_trials = true,
    )

    println("  ✓ 3-way classification complete")
    println("    Max accuracy: $(round(maximum(decoded_3way.average_score), digits=3))")
    println("    Chance level: $(round(decoded_3way.chance_level, digits=3)) (1/3)")
    println("    Conditions: $(join(decoded_3way.condition_names, ", "))")

    # Plot confusion matrix at peak time
    max_idx = argmax(decoded_3way.average_score)
    peak_time = decoded_3way.times[max_idx]
    println("\n    Plotting confusion matrix at peak time ($(round(peak_time, digits=3)) s)...")
    fig3 = eegfun.plot_confusion_matrix(decoded_3way, time_point = peak_time)

catch e
    println("  ✗ Error: $e")
    rethrow(e)
end

# ============================================================================
# SUMMARY
# ============================================================================

println("\n" * "="^70)
println("DECODING TESTS COMPLETE")
println("="^70)
println("\nTests performed:")
println("  ✓ Binary classification with LDA (default)")
println("  ✓ Grand averaging across participants")
println("  ✓ Binary classification with SVM (ERPLAB-like)")
println("  ✓ Multi-class classification (3 conditions)")
println("\nAll plots should be displayed above.")
println("Check that:")
println("  - Decoding accuracy is above chance level")
println("  - Accuracy varies over time (time-point-by-time-point)")
println("  - Grand average shows group-level pattern")
println("  - Confusion matrices show class-specific performance")

