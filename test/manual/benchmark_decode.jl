using eegfun
using DataFrames
using Random
using BenchmarkTools

# Set seed for reproducibility
Random.seed!(42)

# Create test epoch data
n_epochs_per_condition = 30
n_timepoints = 200
sample_rate = 1000
channels = [:Fz, :Cz, :Pz, :POz]

# Create layout
layout = eegfun.Layout(
    DataFrame(label = channels, inc = [90.0, 0.0, -90.0, -90.0], azi = [0.0, 0.0, 0.0, 0.0]),
    nothing,
    nothing,
)

analysis_info = eegfun.AnalysisInfo(:none, 0.0, 0.0)
time = collect(range(-0.2, 0.8, length = n_timepoints))

# Function to create epochs
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

        for (ch_idx, ch) in enumerate(channels)
            if condition_id == 1
                signal =
                    signal_strength * exp.(-((time .- 0.2) .^ 2) / (2 * 0.05^2)) .+
                    0.3 * signal_strength * exp.(-((time .- 0.4) .^ 2) / (2 * 0.08^2))
            else
                signal =
                    -signal_strength * exp.(-((time .- 0.2) .^ 2) / (2 * 0.05^2)) .+
                    0.2 * signal_strength * exp.(-((time .- 0.5) .^ 2) / (2 * 0.1^2))
            end

            channel_factor = 1.0 + (ch_idx - 1) * 0.1
            noise = 0.3 * randn(n_timepoints)
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

# Create epochs for participant 1
cond1 = create_participant_condition_epochs(1, 1, "Face", n_epochs_per_condition, 0.1)
cond2 = create_participant_condition_epochs(1, 2, "Object", n_epochs_per_condition, 0.1)
epochs_p1 = [cond1, cond2]

println("Running benchmark...")
println("Configuration: 10 iterations, 3 folds, 200 timepoints, 4 channels, 60 trials")
println()

# Benchmark
result = @btime eegfun.decode_libsvm($epochs_p1; n_iterations = 10, n_folds = 3, equalize_trials = true)

println()
println("Decoding complete!")
println("Max accuracy: $(round(maximum(result.average_score), digits=3))")
