#!/usr/bin/env julia

# Simple script to run configuration tests
# Usage: julia test/run_config_tests.jl

using Pkg

# Activate the project environment
Pkg.activate(".")

# Add Test package if not already available
try
    using Test
catch
    Pkg.add("Test")
    using Test
end

# Run the configuration tests
println("Running Configuration System Tests...")
println("=" ^ 50)

include("test_config.jl")

println("\nConfiguration tests completed!") 