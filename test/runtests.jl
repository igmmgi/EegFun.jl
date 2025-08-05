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
    include("test_config.jl")
    include("test_layout.jl")  # Commented out due to failures
    include("test_types.jl")   # Commented out due to failures
    include("test_data.jl")    # Commented out due to failures
    include("test_print.jl")
    # include("test_io.jl")
    # include("test_processing.jl")
    # include("test_visualization.jl")
end

println("\nAll tests completed!")
