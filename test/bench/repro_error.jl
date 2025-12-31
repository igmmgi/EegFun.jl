using eegfun
using DataFrames
using FFTW

# Mock EpochData
n_samples = 1250
n_channels = 72
n_trials = 160
times = range(-0.5, stop=2, length=n_samples)
sr = 500.0

data = [DataFrame(Dict(Symbol("CH$i") => rand(n_samples) for i in 1:n_channels)) for _ in 1:n_trials]

struct MockEpochs
    file::String
    condition::Int
    condition_name::String
    data::Vector{DataFrame}
    layout::Any
    sample_rate::Float64
    times::StepRangeLen
    analysis_info::Dict
end

# We can't easily mock the struct because eegfun uses it in dispatch.
# Let's try to use the load_data if it works, or just create a real one.

# Use existing data if available to be realistic
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
if !isempty(epoch_files)
    println("Loading real data for repro...")
    epochs = eegfun.load_data(joinpath(data_dir, epoch_files[1]))[1]
    foi = 1:1:40
    toi = -0.5:0.01:2
    try
        println("Calling tf_analysis (FORCE SEQUENTIAL)...")
        @time tf_data = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)
        println("Success!")
    catch e
        println("Error caught: $e")
        Base.showerror(stdout, e, catch_backtrace())
    end
else
    println("No real data found, skipping repro.")
end
