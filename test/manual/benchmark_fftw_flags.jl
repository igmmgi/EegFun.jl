using BenchmarkTools
using FFTW
using LinearAlgebra

# Simulation parameters matching tf_hanning usage
timewinidx = 128
n_trials = 160
n_executions = 251

# Input data
tmpdat = randn(Float64, timewinidx, n_trials)

println("Benchmarking FFTW.ESTIMATE vs FFTW.MEASURE")
println("Plan dimensions: ($timewinidx, $n_trials), Executions: $n_executions")

# 1. ESTIMATE (Default)
println("\n--- ESTIMATE ---")
print("Planning time: ")
@btime p_est = plan_rfft($tmpdat, 1; flags=FFTW.ESTIMATE)

p_est = plan_rfft(tmpdat, 1; flags=FFTW.ESTIMATE)
out = p_est * tmpdat # Warmup
out_buf = zeros(ComplexF64, div(timewinidx, 2) + 1, n_trials)

print("Execution time (per call): ")
t_est = @btime mul!($out_buf, $p_est, $tmpdat)
total_est = n_executions * t_est

# 2. MEASURE
println("\n--- MEASURE ---")
print("Planning time: ")
@btime p_meas = plan_rfft($tmpdat, 1; flags=FFTW.MEASURE)

p_meas = plan_rfft(tmpdat, 1; flags=FFTW.MEASURE)

print("Execution time (per call): ")
t_meas = @btime mul!($out_buf, $p_meas, $tmpdat)
total_meas = n_executions * t_meas

println("\n--- RESULTS ---")
println("Total Execution Time for $n_executions calls:")
println("ESTIMATE: $(total_est) ns")
println("MEASURE:  $(total_meas) ns")
println("Ratio (est/meas): $(total_est/total_meas)")

if total_meas < total_est
    println("MEASURE is faster!")
else
    println("ESTIMATE is faster (or equal)!")
end
