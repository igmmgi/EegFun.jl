using Test
using Dates
using OrderedCollections

@testset "Print Utilities" begin
    @testset "Vector Printing" begin
        # Test short vector
        v = [1, 2, 3]
        @test EegFun.print_vector(v) == "1, 2, 3"

        # Test long vector
        v = collect(1:20)
        result = EegFun.print_vector(v)
        @test startswith(result, "1, 2, 3, 4, 5, ...")
        @test endswith(result, "16, 17, 18, 19, 20")

        # Test UnitRange
        v = 1:20
        result = EegFun.print_vector(v)
        @test startswith(result, "1, 2, 3, 4, 5, ...")
        @test endswith(result, "16, 17, 18, 19, 20")

        # Test custom max_length and n_ends
        v = collect(1:20)
        result = EegFun.print_vector(v, max_length = 15, n_ends = 3)
        @test startswith(result, "1, 2, 3, ...")
        @test endswith(result, "18, 19, 20")

        # Test edge cases
        @test EegFun.print_vector(Int[]) == "[]"  # Empty vector
        @test EegFun.print_vector([1]) == "1"  # Single element
        @test EegFun.print_vector([1, 2]) == "1, 2"  # Two elements
        @test EegFun.print_vector([1, 2, 3, 4, 5]) == "1, 2, 3, 4, 5"  # Exactly max_length

        # Test with different data types
        @test EegFun.print_vector([:a, :b, :c]) == "a, b, c"  # Symbols
        @test EegFun.print_vector(["a", "b", "c"]) == "a, b, c"  # Strings (no quotes added)
        @test EegFun.print_vector([1.5, 2.5, 3.5]) == "1.5, 2.5, 3.5"  # Floats
        @test EegFun.print_vector([true, false, true]) == "true, false, true"  # Booleans

        # Test with very long vectors
        v = collect(1:100)
        result = EegFun.print_vector(v, max_length = 5, n_ends = 2)
        @test startswith(result, "1, 2, ...")
        @test endswith(result, "99, 100")
        @test count(==("..."), split(result, ", ")) == 1  # Only one ellipsis

        # Test with max_length smaller than n_ends
        v = collect(1:10)
        result = EegFun.print_vector(v, max_length = 5, n_ends = 10)
        # When n_ends > max_length, it shows all elements since vector is shorter than n_ends
        @test result == "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ..., 1, 2, 3, 4, 5, 6, 7, 8, 9, 10"

        # Test with zero n_ends (should show only first part)
        v = collect(1:10)
        result = EegFun.print_vector(v, max_length = 5, n_ends = 0)
        @test result == "..."  # When n_ends is 0, only ellipsis is shown

        # Test with negative n_ends (should throw error)
        v = collect(1:10)
        @test_throws ArgumentError EegFun.print_vector(v, max_length = 5, n_ends = -2)
    end

    @testset "Version" begin
        version = EegFun.get_package_version(package_name = "EegFun")
        @test typeof(version) === String
        @test !isempty(version)

        # Test that it returns a valid version format or "unknown"
        if version != "unknown"
            @test occursin(r"^\d+\.\d+\.\d+", version)  # Should match semver pattern
        end
    end

    @testset "Config Printing" begin
        # Test basic config
        config = Dict("test" => Dict("value" => 1, "array" => [1, 2, 3], "nested" => Dict("key" => "value")))

        # Test printing to string
        io = IOBuffer()
        EegFun.print_config(config, io)
        output = String(take!(io))
        @test contains(output, "test")
        @test contains(output, "value = 1")
        @test contains(output, "array = [1, 2, 3]")
        @test contains(output, "nested")
        @test contains(output, "key = \"value\"")
        @test contains(output, "metadata")
        @test contains(output, "date")
        @test contains(output, "julia_version")
        @test contains(output, "EegFun_version")

        # Test printing to file
        test_file = "test_config.toml"
        EegFun.print_config(config, test_file)
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
        empty_config = Dict{String,Any}()
        io = IOBuffer()
        EegFun.print_config(empty_config, io)
        output = String(take!(io))
        @test contains(output, "metadata")
        @test !contains(output, "test")

        # Test with config containing metadata
        config_with_meta = Dict("metadata" => Dict("old" => "value"), "test" => Dict("value" => 1))
        io = IOBuffer()
        EegFun.print_config(config_with_meta, io)
        output = String(take!(io))
        @test contains(output, "metadata")
        @test !contains(output, "old = \"value\"")  # Old metadata should be replaced
        @test contains(output, "test")
        @test contains(output, "value = 1")

        # Test with complex nested structures
        complex_config = Dict(
            "section1" => Dict(
                "string_val" => "hello world",
                "int_val" => 42,
                "float_val" => 3.14,
                "bool_val" => true,
                "array_val" => [1, 2, 3, 4, 5],
                "nested" => Dict("deep_key" => "deep_value", "numbers" => [10, 20, 30]),
            ),
            "section2" => Dict("empty_array" => Int[], "empty_dict" => Dict{String,Any}(), "mixed_array" => [1, "string", 3.14, true]),
        )

        io = IOBuffer()
        EegFun.print_config(complex_config, io)
        output = String(take!(io))

        # Test that all values are properly formatted
        @test contains(output, "string_val = \"hello world\"")
        @test contains(output, "int_val = 42")
        @test contains(output, "float_val = 3.14")
        @test contains(output, "bool_val = true")
        @test contains(output, "array_val = [1, 2, 3, 4, 5]")
        @test contains(output, "deep_key = \"deep_value\"")
        @test contains(output, "numbers = [10, 20, 30]")
        @test contains(output, "empty_array = []")
        @test contains(output, "[section2.empty_dict]")  # Empty dicts become empty sections in TOML
        @test contains(output, "mixed_array = [1, \"string\", 3.14, true]")

        # Test metadata structure
        @test contains(output, "[metadata]")
        @test contains(output, "date = ")
        @test contains(output, "julia_version = ")
        @test contains(output, "EegFun_version = ")

        # Test with Symbol keys (should be converted to strings)
        symbol_config = Dict(:symbol_key => Dict(:nested_symbol => "value"))
        io = IOBuffer()
        EegFun.print_config(symbol_config, io)
        output = String(take!(io))
        @test contains(output, "symbol_key")
        @test contains(output, "nested_symbol = \"value\"")

        # Test error handling for file operations
        # Test with invalid filename (directory that doesn't exist)
        invalid_file = "/nonexistent/directory/test.toml"
        @test_throws SystemError EegFun.print_config(config, invalid_file)

        # Test with read-only directory (if possible)
        # This would require creating a read-only directory, which is complex
        # and platform-dependent, so we'll skip this test for now
    end

    @testset "Edge Cases and Error Handling" begin
        # Test print_vector with various edge cases
        @test EegFun.print_vector([NaN]) == "NaN"
        @test EegFun.print_vector([Inf]) == "Inf"
        @test EegFun.print_vector([-Inf]) == "-Inf"
        @test EegFun.print_vector([missing]) == "missing"
        @test EegFun.print_vector([nothing]) == "nothing"

        # Test with very large numbers
        @test EegFun.print_vector([1e10, 1e-10]) == "1.0e10, 1.0e-10"

        # Test with special characters in strings
        @test EegFun.print_vector(["a\nb", "c\td"]) == "a\nb, c\td"  # No escaping in print_vector

        # Test config printing with TOML-compatible data types
        edge_config = Dict(
            "edge_cases" => Dict(
                "nan_val" => NaN,
                "inf_val" => Inf,
                "neg_inf_val" => -Inf,
                "large_int" => typemax(Int64),
                "small_int" => typemin(Int64),
                "large_float" => 1e308,
                "small_float" => 1e-308,
            ),
        )

        io = IOBuffer()
        EegFun.print_config(edge_config, io)
        output = String(take!(io))

        # Test that special values are handled properly
        @test contains(output, "nan_val = nan")  # TOML uses lowercase
        @test contains(output, "inf_val = +inf")  # TOML uses +inf
        @test contains(output, "neg_inf_val = -inf")  # TOML uses -inf
        @test contains(output, "large_int = ")
        @test contains(output, "small_int = ")
        @test contains(output, "large_float = ")
        @test contains(output, "small_float = ")

    end

    @testset "Performance and Memory" begin
        # Test print_vector with very large vectors
        large_vector = collect(1:10000)
        result = EegFun.print_vector(large_vector, max_length = 10, n_ends = 5)
        @test length(result) < 100  # Should be truncated
        @test startswith(result, "1, 2, 3, 4, 5, ...")
        @test endswith(result, "9996, 9997, 9998, 9999, 10000")

        # Test config printing with large nested structures
        large_config = Dict("large_section" => Dict("array" => collect(1:100)))
        io = IOBuffer()
        EegFun.print_config(large_config, io)
        output = String(take!(io))
        @test contains(output, "large_section")
        @test contains(
            output,
            "array = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]",
        )
    end
end
