using Test
using EegFun

@testset "File Utilities" begin
    # Create temporary test directory
    test_dir = mktempdir()

    # Create test files for various scenarios
    test_files = [
        "1_epochs_cleaned.jld2",
        "2_epochs_cleaned.jld2",
        "3_epochs_cleaned.jld2",
        "1_erps_cleaned.jld2",
        "2_erps_cleaned.jld2",
        "Flank_C_3_epochs_cleaned.jld2",
        "Flank_C_4_epochs_cleaned.jld2",
        "Flank_C_5_epochs_cleaned.jld2",
        "test_file.csv",
        "another_file.txt",
    ]

    # Create the test files
    for file in test_files
        touch(joinpath(test_dir, file))
    end

    # Create subdirectories for recursive search tests
    subdir1 = joinpath(test_dir, "subdir1")
    subdir2 = joinpath(test_dir, "subdir2", "nested")
    mkpath(subdir1)
    mkpath(subdir2)

    # Create files in subdirectories
    touch(joinpath(subdir1, "biosemi64.csv"))
    touch(joinpath(subdir1, "config.toml"))
    touch(joinpath(subdir2, "biosemi32.csv"))
    touch(joinpath(subdir2, "layout.csv"))


    @testset "check_files_exist with Vector{String}" begin
        # Test with existing files
        existing_files = [joinpath(test_dir, "1_epochs_cleaned.jld2"), joinpath(test_dir, "2_epochs_cleaned.jld2")]
        result = EegFun.check_files_exist(existing_files)
        @test result == true

        # Test with non-existing files
        missing_files = [joinpath(test_dir, "nonexistent1.jld2"), joinpath(test_dir, "nonexistent2.jld2")]
        result = EegFun.check_files_exist(missing_files)
        @test result == false

        # Test with mixed existing and missing files
        mixed_files = [joinpath(test_dir, "1_epochs_cleaned.jld2"), joinpath(test_dir, "nonexistent.jld2")]
        result = EegFun.check_files_exist(mixed_files)
        @test result == false
    end


    @testset "get_files with String pattern" begin
        # Test with regex pattern
        files = EegFun.get_files(test_dir, ".*epochs_cleaned.*")
        @test length(files) == 6  # 3 numbered files + 3 Flank files
        @test all(contains.(files, "epochs_cleaned"))

        # Test with specific pattern
        files = EegFun.get_files(test_dir, "1_.*")
        @test length(files) == 2  # 1_epochs_cleaned.jld2 and 1_erps_cleaned.jld2
        @test all(contains.(files, "1_"))

        # Test with no matches
        files = EegFun.get_files(test_dir, "nonexistent.*")
        @test isempty(files)
    end

    @testset "get_files with Vector{String}" begin
        # Test with specific filenames
        specific_files = ["1_epochs_cleaned.jld2", "2_epochs_cleaned.jld2"]
        files = EegFun.get_files(test_dir, specific_files)
        @test length(files) == 2
        @test all(isfile.(files))

        # Test with non-existing files
        missing_files = ["nonexistent1.jld2", "nonexistent2.jld2"]
        files = EegFun.get_files(test_dir, missing_files)
        @test length(files) == 2  # Still returns paths, even if files don't exist
        @test all(contains.(files, "nonexistent"))
    end

    @testset "find_file" begin
        # Test exact match in root directory
        result = EegFun.find_file("test_file.csv", test_dir)
        @test result == joinpath(test_dir, "test_file.csv")

        # Test exact match in subdirectory (recursive)
        result = EegFun.find_file("biosemi64.csv", test_dir)
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test with extensions
        result = EegFun.find_file("biosemi64", test_dir, extensions = [".csv"])
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test with multiple extensions
        result = EegFun.find_file("config", test_dir, extensions = [".toml", ".csv"])
        @test result == joinpath(subdir1, "config.toml")

        # Test non-recursive search
        result = EegFun.find_file("biosemi64.csv", test_dir, recursive = false)
        @test result === nothing  # File is in subdirectory

        # Test non-recursive search in subdirectory
        result = EegFun.find_file("biosemi64.csv", subdir1, recursive = false)
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test file not found
        result = EegFun.find_file("nonexistent.csv", test_dir)
        @test result === nothing

        # Test directory doesn't exist
        result = EegFun.find_file("test.csv", "/nonexistent/directory")
        @test result === nothing

        # Test with empty extensions
        result = EegFun.find_file("test_file.csv", test_dir, extensions = String[])
        @test result == joinpath(test_dir, "test_file.csv")
    end


    @testset "Edge cases and error handling" begin
        # Test with empty vectors
        result = EegFun.check_files_exist(String[])
        @test result == true  # Empty list should return true


        # Test find_file with empty filename
        result = EegFun.find_file("", test_dir)
        @test result === nothing

        # Test get_files with empty directory
        empty_dir = mktempdir()
        files = EegFun.get_files(empty_dir, ".*")
        @test isempty(files)

        # Clean up empty directory
        rm(empty_dir)
    end

    # Clean up test directory
    rm(test_dir, recursive = true, force = true)
end
