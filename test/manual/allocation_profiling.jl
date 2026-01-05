# Quick guide to profiling allocations in Julia

# Method 1: @timev - shows detailed allocation breakdown
# Usage:
# @timev tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01))

# Method 2: @allocated - measure allocations of a specific expression
# Usage:
# allocs = @allocated tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01))
# println("Allocations: ", allocs)

# Method 3: Profile.Allocs - BEST for finding WHERE allocations occur
# First, install ProfileView if needed: using Pkg; Pkg.add("ProfileView")
using Profile
using ProfileView

# Clear previous profiling data
Profile.clear()
Profile.Allocs.clear()

# Profile allocations (sample_rate=1 means sample every allocation)
Profile.Allocs.@profile sample_rate=1 begin
    tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01))
end

# View the allocation profile (interactive flame graph)
ProfileView.view_allocs()

# Or print a text summary:
# Profile.Allocs.print()

# Method 4: @code_warntype - check for type instability (causes allocations)
# Usage:
# @code_warntype eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01))
# Look for red text indicating type instability

# Method 5: @btime with detailed output
# Usage:
# @btime tf_data = eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01)) evals=1

# Method 6: JET.jl for static analysis (catches allocations at compile time)
# using JET
# @report_opt eegfun.tf_stft(epochs_synthetic, lin_freqs = (1, 40, 1), window_length = 0.5, time_steps = (-0.5, 1.5, 0.01))

