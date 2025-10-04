using Test
using eegfun

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

    @testset "check_files_exist with Vector{Int}" begin
        # Create test files in current directory for these tests
        test_files_current = ["1_epochs_cleaned.jld2", "2_epochs_cleaned.jld2", "3_epochs_cleaned.jld2"]
        for file in test_files_current
            touch(file)
        end

        try
            # Test with existing files
            existing_files = [1, 2, 3]
            result = eegfun.check_files_exist(existing_files, "epochs_cleaned")
            @test result == true

            # Test with non-existing files
            missing_files = [10, 11, 12]
            result = eegfun.check_files_exist(missing_files, "epochs_cleaned")
            @test result == false

            # Test with mixed existing and missing files
            mixed_files = [1, 10, 2]
            result = eegfun.check_files_exist(mixed_files, "epochs_cleaned")
            @test result == false
        finally
            # Clean up test files
            for file in test_files_current
                rm(file, force = true)
            end
        end
    end

    @testset "check_files_exist with Int" begin
        # Create test file in current directory for this test
        test_file = "1_epochs_cleaned.jld2"
        touch(test_file)

        try
            # Test with existing file
            result = eegfun.check_files_exist(1, "epochs_cleaned")
            @test result == true

            # Test with non-existing file
            result = eegfun.check_files_exist(10, "epochs_cleaned")
            @test result == false
        finally
            # Clean up test file
            rm(test_file, force = true)
        end
    end

    @testset "check_files_exist with Vector{String}" begin
        # Test with existing files
        existing_files = [joinpath(test_dir, "1_epochs_cleaned.jld2"), joinpath(test_dir, "2_epochs_cleaned.jld2")]
        result = eegfun.check_files_exist(existing_files)
        @test result == true

        # Test with non-existing files
        missing_files = [joinpath(test_dir, "nonexistent1.jld2"), joinpath(test_dir, "nonexistent2.jld2")]
        result = eegfun.check_files_exist(missing_files)
        @test result == false

        # Test with mixed existing and missing files
        mixed_files = [joinpath(test_dir, "1_epochs_cleaned.jld2"), joinpath(test_dir, "nonexistent.jld2")]
        result = eegfun.check_files_exist(mixed_files)
        @test result == false
    end

    @testset "check_files_exist with subjects and conditions" begin
        # Create test files in current directory for these tests
        test_files_subj_cond =
            ["1_1_epochs_cleaned.jld2", "1_2_epochs_cleaned.jld2", "2_1_epochs_cleaned.jld2", "2_2_epochs_cleaned.jld2"]
        for file in test_files_subj_cond
            touch(file)
        end

        try
            # Test with existing files
            result = eegfun.check_files_exist([1, 2], [1, 2], "epochs_cleaned")
            @test result == true

            # Test with non-existing files
            result = eegfun.check_files_exist([10, 11], [1, 2], "epochs_cleaned")
            @test result == false

            # Test with single subject and condition
            result = eegfun.check_files_exist(1, 1, "epochs_cleaned")
            @test result == true

            result = eegfun.check_files_exist(10, 1, "epochs_cleaned")
            @test result == false
        finally
            # Clean up test files
            for file in test_files_subj_cond
                rm(file, force = true)
            end
        end
    end

    @testset "get_files with String pattern" begin
        # Test with regex pattern
        files = eegfun.get_files(test_dir, ".*epochs_cleaned.*")
        @test length(files) == 6  # 3 numbered files + 3 Flank files
        @test all(contains.(files, "epochs_cleaned"))

        # Test with specific pattern
        files = eegfun.get_files(test_dir, "1_.*")
        @test length(files) == 2  # 1_epochs_cleaned.jld2 and 1_erps_cleaned.jld2
        @test all(contains.(files, "1_"))

        # Test with no matches
        files = eegfun.get_files(test_dir, "nonexistent.*")
        @test isempty(files)
    end

    @testset "get_files with Vector{String}" begin
        # Test with specific filenames
        specific_files = ["1_epochs_cleaned.jld2", "2_epochs_cleaned.jld2"]
        files = eegfun.get_files(test_dir, specific_files)
        @test length(files) == 2
        @test all(isfile.(files))

        # Test with non-existing files
        missing_files = ["nonexistent1.jld2", "nonexistent2.jld2"]
        files = eegfun.get_files(test_dir, missing_files)
        @test length(files) == 2  # Still returns paths, even if files don't exist
        @test all(contains.(files, "nonexistent"))
    end

    @testset "find_file" begin
        # Test exact match in root directory
        result = eegfun.find_file("test_file.csv", test_dir)
        @test result == joinpath(test_dir, "test_file.csv")

        # Test exact match in subdirectory (recursive)
        result = eegfun.find_file("biosemi64.csv", test_dir)
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test with extensions
        result = eegfun.find_file("biosemi64", test_dir, extensions = [".csv"])
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test with multiple extensions
        result = eegfun.find_file("config", test_dir, extensions = [".toml", ".csv"])
        @test result == joinpath(subdir1, "config.toml")

        # Test non-recursive search
        result = eegfun.find_file("biosemi64.csv", test_dir, recursive = false)
        @test result === nothing  # File is in subdirectory

        # Test non-recursive search in subdirectory
        result = eegfun.find_file("biosemi64.csv", subdir1, recursive = false)
        @test result == joinpath(subdir1, "biosemi64.csv")

        # Test file not found
        result = eegfun.find_file("nonexistent.csv", test_dir)
        @test result === nothing

        # Test directory doesn't exist
        result = eegfun.find_file("test.csv", "/nonexistent/directory")
        @test result === nothing

        # Test with empty extensions
        result = eegfun.find_file("test_file.csv", test_dir, extensions = String[])
        @test result == joinpath(test_dir, "test_file.csv")
    end

    @testset "_filter_files" begin
        # Test files with participant numbers
        test_files_with_participants = [
            "Flank_C_3_epochs_cleaned.jld2",
            "Flank_C_4_epochs_cleaned.jld2",
            "Flank_C_5_epochs_cleaned.jld2",
            "1_epochs_cleaned.jld2",
            "2_epochs_cleaned.jld2",
            "3_epochs_cleaned.jld2",
            "no_participant_file.jld2",
        ]

        # Test include filter
        filtered = eegfun._filter_files(test_files_with_participants, include = [3, 4])
        @test length(filtered) == 4  # Flank_C_3, Flank_C_4, 3_epochs_cleaned, no_participant_file
        @test "Flank_C_3_epochs_cleaned.jld2" in filtered
        @test "Flank_C_4_epochs_cleaned.jld2" in filtered
        @test "3_epochs_cleaned.jld2" in filtered
        @test "no_participant_file.jld2" in filtered

        # Test exclude filter
        filtered = eegfun._filter_files(test_files_with_participants, exclude = [3, 4])
        @test length(filtered) == 4  # Flank_C_5, 1_, 2_, no_participant
        @test "Flank_C_5_epochs_cleaned.jld2" in filtered
        @test "1_epochs_cleaned.jld2" in filtered
        @test "2_epochs_cleaned.jld2" in filtered
        @test "no_participant_file.jld2" in filtered

        # Test both include and exclude
        filtered = eegfun._filter_files(test_files_with_participants, include = [3, 4, 5], exclude = [4])
        @test length(filtered) == 4  # Flank_C_3, 3_epochs_cleaned, Flank_C_5, no_participant_file
        @test "Flank_C_3_epochs_cleaned.jld2" in filtered
        @test "3_epochs_cleaned.jld2" in filtered
        @test "Flank_C_5_epochs_cleaned.jld2" in filtered
        @test "no_participant_file.jld2" in filtered

        # Test with single Int values
        filtered = eegfun._filter_files(test_files_with_participants, include = 3)
        @test length(filtered) == 3  # Flank_C_3, 3_epochs_cleaned, no_participant_file
        @test "Flank_C_3_epochs_cleaned.jld2" in filtered
        @test "3_epochs_cleaned.jld2" in filtered
        @test "no_participant_file.jld2" in filtered

        filtered = eegfun._filter_files(test_files_with_participants, exclude = 3)
        @test length(filtered) == 5  # All others except 3, but including no_participant

        # Test with nothing values (should include all)
        filtered = eegfun._filter_files(test_files_with_participants, include = nothing, exclude = nothing)
        @test length(filtered) == length(test_files_with_participants)

        # Test files with no participant number
        files_no_participant = ["no_participant.jld2", "another_file.jld2"]
        filtered = eegfun._filter_files(files_no_participant, include = [1, 2])
        @test length(filtered) == 2  # Files without participant numbers are included by default

        filtered = eegfun._filter_files(files_no_participant, exclude = [1, 2])
        @test length(filtered) == 2  # All files included since they have no participant numbers
    end

    @testset "Edge cases and error handling" begin
        # Test with empty vectors
        result = eegfun.check_files_exist(Int[], "epochs_cleaned")
        @test result == true  # Empty list should return true

        result = eegfun.check_files_exist(String[])
        @test result == true  # Empty list should return true

        # Test _filter_files with empty input
        filtered = eegfun._filter_files(String[])
        @test isempty(filtered)

        # Test find_file with empty filename
        result = eegfun.find_file("", test_dir)
        @test result === nothing

        # Test get_files with empty directory
        empty_dir = mktempdir()
        files = eegfun.get_files(empty_dir, ".*")
        @test isempty(files)

        # Clean up empty directory
        rm(empty_dir)
    end

    # Clean up test directory
    rm(test_dir, recursive = true, force = true)
end
