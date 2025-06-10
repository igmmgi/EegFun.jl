using Test

# Import the package
using eegfun

println("Running eegfun.jl Test Suite")
println("=" ^ 40)

@testset "eegfun.jl Tests" begin
    
    # Configuration system tests
    @testset "Configuration System" begin
        include("test_config.jl")
    end
    
    # Note: Add other test modules here as they are created
    # @testset "Data Processing" begin
    #     include("test_processing.jl")
    # end
    
    # @testset "Filtering" begin
    #     include("test_filtering.jl") 
    # end
    
    # @testset "ICA" begin
    #     include("test_ica.jl")
    # end
    
    # @testset "Plotting" begin
    #     include("test_plotting.jl")
    # end
    
end

println("\nAll tests completed!") 