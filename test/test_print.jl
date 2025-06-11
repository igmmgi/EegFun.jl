using Test
using Dates
using OrderedCollections
using eegfun

@testset "Print Utilities" begin
    # Test vector printing
    @testset "Vector Printing" begin
        # Test short vector
        v = [1, 2, 3]
        @test eegfun.print_vector(v) == "1, 2, 3"
        
        # Test long vector
        v = collect(1:20)
        result = eegfun.print_vector(v)
        @test startswith(result, "1, 2, 3, 4, 5, ...")
        @test endswith(result, "16, 17, 18, 19, 20")
        
        # Test UnitRange
        v = 1:20
        result = eegfun.print_vector(v)
        @test startswith(result, "1, 2, 3, 4, 5, ...")
        @test endswith(result, "16, 17, 18, 19, 20")
        
        # Test custom max_length and n_ends
        v = collect(1:20)
        result = eegfun.print_vector(v, max_length=15, n_ends=3)
        @test startswith(result, "1, 2, 3, ...")
        @test endswith(result, "18, 19, 20")
    end

    # Test git commit
    commit = eegfun.get_git_commit()
    @test typeof(commit) === String
    @test !isempty(commit)

    # Test version
    version = eegfun.get_eegfun_version()
    @test typeof(version) === String
    @test !isempty(version)

    # Test config printing
    @testset "Config Printing" begin
        # Test basic config
        config = Dict(
            "test" => Dict(
                "value" => 1,
                "array" => [1, 2, 3],
                "nested" => Dict("key" => "value")
            )
        )
        
        # Test printing to string
        io = IOBuffer()
        eegfun.print_config(config, io)
        output = String(take!(io))
        @test contains(output, "test")
        @test contains(output, "value = 1")
        @test contains(output, "array = [1, 2, 3]")
        @test contains(output, "nested")
        @test contains(output, "key = \"value\"")
        @test contains(output, "metadata")
        @test contains(output, "generated_at")
        @test contains(output, "eegfun_version")
        @test contains(output, "git_commit")
        
        # Test printing to file
        test_file = "test_config.toml"
        eegfun.print_config(config, test_file)
        @test isfile(test_file)
        
        # Verify file contents
        file_content = read(test_file, String)
        @test contains(file_content, "test")
        @test contains(file_content, "value = 1")
        @test contains(file_content, "array = [1, 2, 3]")
        @test contains(file_content, "nested")
        @test contains(file_content, "key = \"value\"")
        @test contains(file_content, "metadata")
        
        # Clean up
        rm(test_file)
        
        # Test with empty config
        empty_config = Dict{String, Any}()
        io = IOBuffer()
        eegfun.print_config(empty_config, io)
        output = String(take!(io))
        @test contains(output, "metadata")
        @test !contains(output, "test")
        
        # Test with config containing metadata
        config_with_meta = Dict(
            "metadata" => Dict("old" => "value"),
            "test" => Dict("value" => 1)
        )
        io = IOBuffer()
        eegfun.print_config(config_with_meta, io)
        output = String(take!(io))
        @test contains(output, "metadata")
        @test !contains(output, "old = \"value\"")  # Old metadata should be replaced
        @test contains(output, "test")
        @test contains(output, "value = 1")
    end
end 