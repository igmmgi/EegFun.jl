using eegfun
using BenchmarkTools

# Load epoched data
data_dir = "/home/ian/Documents/Julia/output_data"
epoch_files = filter(f -> endswith(f, "_epochs_cleaned.jld2"), readdir(data_dir))
if isempty(epoch_files)
    println("No epoch files found in $data_dir")
    exit(1)
end
epoch_file = joinpath(data_dir, epoch_files[1])
println("Loading $epoch_file...")
epochs = eegfun.load_data(epoch_file)[1]

# Time and frequency parameters
toi = -0.5:0.01:2
foi = 1:1:40

println("Running benchmark for :wavelet method...")
# First run for compilation
eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)

# Actual benchmark
@btime tf_data1 = eegfun.tf_analysis(epochs, foi, toi; method=:wavelet, width=7)

println("\nNumber of channels: $(length(eegfun.get_selected_channels(epochs, eegfun.channels(); include_meta=false, include_extra=false)))")
println("Number of trials: $(length(epochs.data))")
println("Threads used: $(Threads.nthreads())")
