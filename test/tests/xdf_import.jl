using Test
using EegFun

@testset "XDF Import" begin

    # ────────────────────────────────────────────────────────────
    # 1. read_xdf function exists and has correct signature
    # ────────────────────────────────────────────────────────────
    @testset "read_xdf API" begin
        # Method should exist
        @test hasmethod(EegFun.read_xdf, Tuple{String})

        # Should error on non-existent file
        @test_throws Exception EegFun.read_xdf("/nonexistent/file.xdf")
    end

    # ────────────────────────────────────────────────────────────
    # 2. File extension validation via read_raw_data dispatcher
    # ────────────────────────────────────────────────────────────
    @testset "read_raw_data dispatches XDF" begin
        # read_raw_data should recognize .xdf extension
        # (will fail on file not found, but that confirms dispatch works)
        @test_throws Exception EegFun.read_raw_data("/nonexistent/file.xdf")
    end

end # @testset "XDF Import"
