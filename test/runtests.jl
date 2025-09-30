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

    include("artifact_detection.jl")
    include("baseline.jl")
    include("channel_average.jl")
    include("channel_summary.jl")
    include("channel_difference.jl")
    include("config.jl")
    include("channel_metrics.jl")
    include("data.jl")
    include("layout.jl")  
    include("layout_system.jl")
    include("plot_channel_summary.jl")
    include("plot_layout.jl")
    include("print.jl")
    include("filter.jl")
    include("epochs.jl")
    include("ica.jl")
    include("rereference.jl")
    include("types.jl")   
    include("triggers.jl")   
    include("plot_triggers.jl")
    
    # Batch processing tests
    include("batch_filter.jl")
    include("batch_average_epochs.jl")
    include("batch_channel_summary.jl")
    include("batch_combine_channels.jl")
    include("batch_combine_conditions.jl")
    include("batch_difference_conditions.jl")
    include("batch_erp_measurements.jl")
    include("batch_grandaverage.jl")
    include("batch_rereference.jl")
    include("batch_utils.jl")
end

println("\nAll tests completed!")
