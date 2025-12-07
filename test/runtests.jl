using Test

# Import the package
using eegfun

# Include shared test utilities first
include("test_utils.jl")

# Clean up old coverage files at the start
if Base.JLOptions().code_coverage != 0
    println("\nCleaning up old coverage files...")
    using Coverage
    Coverage.clean_folder(joinpath(@__DIR__, "..", "src"))
    Coverage.clean_folder(@__DIR__)
end

println("Running eegfun.jl Test Suite")
println("=" ^ 40)

@testset "eegfun" begin
 
    include("automatic/apply.jl")
    include("automatic/apply.jl")
    include("automatic/artifact_detection.jl")
    include("automatic/baseline.jl")
    include("automatic/batch_utils.jl") # TODO: tests seem iffy here!
    include("automatic/channel_average.jl")
    include("automatic/channel_difference.jl")
    include("automatic/channel_metrics.jl")
    include("automatic/channel_repair.jl")
    include("automatic/channel_summary.jl")
    include("automatic/condition_combine.jl")
    include("automatic/condition_difference.jl")
    include("automatic/config.jl")
    include("automatic/data.jl")
    include("automatic/epochs.jl")
    include("automatic/erp_measurements.jl")
    include("automatic/files.jl")
    include("automatic/filter.jl")
    include("automatic/gfp.jl")
    include("automatic/grand_average.jl")
    include("automatic/ica.jl")
    include("automatic/jackknife_average.jl")
    include("automatic/layout.jl")
    include("automatic/layout_system.jl")
    include("automatic/logging.jl")
    include("automatic/lrp.jl")
    include("automatic/mirror.jl")
    include("automatic/misc.jl")
    include("automatic/plot_channel_summary.jl")
    include("automatic/plot_layout.jl")
    include("automatic/plot_misc.jl")
    include("automatic/plot_triggers.jl")
    include("automatic/print.jl")
    include("automatic/realign.jl")
    include("automatic/rereference.jl")
    include("automatic/resample.jl")
    include("automatic/triggers.jl")
    include("automatic/types.jl")
    include("automatic/utils.jl")

end

println("\nAll tests completed!")
