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
    # include("test_config.jl")
    # include("test_layout.jl")  
    # include("test_types.jl")   
    # include("test_data.jl")    
    # include("test_print.jl")
    include("test_channel_average.jl")
end

println("\nAll tests completed!")
