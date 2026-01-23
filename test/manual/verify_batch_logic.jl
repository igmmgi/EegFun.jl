using Test

# Mocking necessary environment to test internal functions
module MockEnv
# Mocking what's needed from eegfun
abstract type EegFunData end
abstract type EegData <: EegFunData end
struct InfoIca <: EegFunData end
struct ErpData <: EegData
    file::String
    condition::Int
    condition_name::String
end
end

@testset "batch.jl EegFunData Verification" begin

    # Test _load_data and _single_or_vector logic with EegFunData
    @testset "load_all_data core logic with EegFunData" begin

        abstract type MockEegFunData end
        abstract type MockEegData <: MockEegFunData end
        struct MockErpData <: MockEegData end
        struct MockInfoIca <: MockEegFunData end

        function mock_load_data(path)
            if path == "erp.jld2"
                return MockErpData()
            elseif path == "ica.jld2"
                return MockInfoIca()
            elseif path == "mixed.jld2"
                return [MockErpData(), MockInfoIca()]
            else
                return nothing
            end
        end

        # New simplified logic for _load_data(Dict) would look like this:
        function test_load_all_data_core(::Type{T}, files) where {T}
            all_data = T[]
            for file in files
                file_data = mock_load_data(file)
                isnothing(file_data) && continue
                if file_data isa Vector{<:T}
                    append!(all_data, file_data)
                elseif file_data isa T
                    push!(all_data, file_data)
                end
            end
            return all_data
        end

        files = ["erp.jld2", "ica.jld2"]
        data = test_load_all_data_core(MockEegFunData, files)
        @test length(data) == 2
        @test data[1] isa MockErpData
        @test data[2] isa MockInfoIca

        # Test mixed case
        mixed_data = test_load_all_data_core(MockEegFunData, ["mixed.jld2"])
        @test length(mixed_data) == 2
        @test mixed_data[1] isa MockErpData
        @test mixed_data[2] isa MockInfoIca
    end
end
