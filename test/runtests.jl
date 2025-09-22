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
    include("test_artifact_detection.jl")
    include("test_baseline.jl")
    include("test_channel_average.jl")
    include("test_channel_summary.jl")
    include("test_channel_difference.jl")
    include("test_config.jl")
    include("test_channel_metrics.jl")
    include("test_data.jl")
    include("test_layout.jl")  
    include("test_layout_system.jl")
    include("test_plot_channel_summary.jl")
    include("test_plot_layout.jl")
    include("test_print.jl")
    include("test_filter.jl")
    include("test_epochs.jl")
    include("test_ica.jl")
    include("test_rereference.jl")
    include("test_types.jl")   
    include("test_triggers.jl")   
    include("test_plot_triggers.jl")
end

println("\nAll tests completed!")
