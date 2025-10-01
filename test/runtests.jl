using Test

# Import the package
using eegfun

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

    include("tests/artifact_detection.jl")
    include("tests/baseline.jl")
    include("tests/channel_average.jl")
    include("tests/channel_difference.jl")
    include("tests/channel_metrics.jl")
    include("tests/channel_summary.jl")
    include("tests/config.jl")
    include("tests/data.jl")
    include("tests/epochs.jl")
    include("tests/filter.jl")
    include("tests/ica.jl")
    include("tests/layout.jl")
    include("tests/layout_system.jl")
    include("tests/files.jl")
    include("tests/misc.jl")
    include("tests/logging.jl")
    include("tests/types.jl")
    include("tests/data.jl")
    include("tests/plot_channel_summary.jl")
    include("tests/plot_layout.jl")
    include("tests/plot_triggers.jl")
    include("tests/print.jl")
    include("tests/rereference.jl")
    include("tests/triggers.jl")
    include("tests/types.jl")
    include("tests/utils.jl")

end

println("\nAll tests completed!")
