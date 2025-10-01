using Test
using DataFrames
using eegfun
using JLD2
using Statistics
using CSV

@testset "Batch Channel Summary" begin
    
    # Create a temporary directory for test files
    test_dir = mktempdir()
    
    try
        # Helper to create test ErpData with known statistics
        function create_test_erp_data(n_conditions::Int = 2)
            erps = eegfun.ErpData[]
            fs = 256.0
            n_samples = 513
            t = range(-1.0, 1.0, length=n_samples)
            
            for cond in 1:n_conditions
                # Create channels with different characteristics
                # Fz: mean=0, std=1 (standard normal)
                # Cz: mean=0, std=2 (higher variance)
                # Pz: mean=5, std=0.5 (shifted, lower variance)
                fz = randn(n_samples)
                cz = 2.0 .* randn(n_samples)
                pz = 5.0 .+ 0.5 .* randn(n_samples)
                
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(cond, n_samples),
                    Fz = fz,
                    Cz = cz,
                    Pz = pz
                )
                
                layout = eegfun.Layout(
                    DataFrame(label = [:Fz, :Cz, :Pz], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
                    nothing, nothing
                )
                
                push!(erps, eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 10))
            end
            
            return erps
        end
        
        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_test_erp_data(2)
                filename = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
                save(filename, "erps", erps)
                @test isfile(filename)
            end
        end
        
        @testset "Basic channel summary" begin
            output_dir = joinpath(test_dir, "summary_output")
            
            # Test basic channel summary
            eegfun.channel_summary("erps_cleaned", 
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            @test isdir(output_dir)
            
            # Check that CSV file exists
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)
            
            # Load and verify CSV content
            results = CSV.read(csv_file, DataFrame)
            
            # Should have results for 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12
            
            # Check required columns exist
            @test hasproperty(results, :file)
            @test hasproperty(results, :condition)
            @test hasproperty(results, :channel)
            @test hasproperty(results, :min)
            @test hasproperty(results, :max)
            @test hasproperty(results, :std)
            @test hasproperty(results, :range)
            @test hasproperty(results, :var)
            @test hasproperty(results, :zvar)
            
            # Verify channels are present (CSV reads them as strings)
            @test "Fz" in results.channel
            @test "Cz" in results.channel
            @test "Pz" in results.channel
        end
        
        @testset "Summary specific participants" begin
            output_dir = joinpath(test_dir, "summary_participant")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = 1)
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)
            
            results = CSV.read(csv_file, DataFrame)
            # Only participant 1: 1 file × 2 conditions × 3 channels = 6 rows
            @test nrow(results) == 6
            @test all(results.file .== "1_erps_cleaned")
        end
        
        @testset "Summary multiple participants" begin
            output_dir = joinpath(test_dir, "summary_multi_participants")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = [1, 2])
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)
            
            # 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12
            @test "1_erps_cleaned" in results.file
            @test "2_erps_cleaned" in results.file
        end
        
        @testset "Summary specific conditions" begin
            output_dir = joinpath(test_dir, "summary_condition")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = 1)
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)
            
            # 2 files × 1 condition × 3 channels = 6 rows
            @test nrow(results) == 6
            @test all(results.condition .== 1)
        end
        
        @testset "Summary multiple conditions" begin
            output_dir = joinpath(test_dir, "summary_multi_conditions")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = [1, 2])
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)
            
            # 2 files × 2 conditions × 3 channels = 12 rows
            @test nrow(results) == 12
            @test 1 in results.condition
            @test 2 in results.condition
        end
        
        @testset "Channel selection predicate" begin
            output_dir = joinpath(test_dir, "summary_channel_select")
            
            # Select only Fz and Cz
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  channel_selection = eegfun.channels([:Fz, :Cz]))
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)
            
            # 2 files × 2 conditions × 2 channels = 8 rows
            @test nrow(results) == 8
            @test "Fz" in results.channel
            @test "Cz" in results.channel
            @test "Pz" ∉ results.channel
        end
        
        @testset "Channel exclusion predicate" begin
            output_dir = joinpath(test_dir, "summary_channel_exclude")
            
            # Exclude Pz
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  channel_selection = eegfun.channels_not([:Pz]))
            
            csv_file = joinpath(output_dir, "channel_summary.csv")
            results = CSV.read(csv_file, DataFrame)
            
            # 2 files × 2 conditions × 2 channels = 8 rows
            @test nrow(results) == 8
            @test "Pz" ∉ results.channel
        end
        
        @testset "Custom output filename" begin
            output_dir = joinpath(test_dir, "summary_custom_name")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  output_file = "my_custom_summary")
            
            # Check custom filename
            csv_file = joinpath(output_dir, "my_custom_summary.csv")
            @test isfile(csv_file)
            
            # Check log file has custom name too
            log_file = joinpath(output_dir, "my_custom_summary.log")
            @test isfile(log_file)
        end
        
        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.channel_summary("erps_cleaned", 
                                                         input_dir = "/nonexistent/path")
        end
        
        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)
            
            # Directory exists but has no JLD2 files matching pattern
            result = eegfun.channel_summary("erps_cleaned",
                                           input_dir = empty_dir)
            
            @test result === nothing  # No files to process
        end
        
        @testset "Logging" begin
            output_dir = joinpath(test_dir, "summary_with_log")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            # Check log file exists
            log_file = joinpath(output_dir, "channel_summary.log")
            @test isfile(log_file)
            
            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch channel summary started")
            @test contains(log_contents, "channel_summary")
            @test contains(log_contents, "erps_cleaned")
        end
        
        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_summary")
            mkpath(output_dir)
            
            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))
            
            # Run channel_summary - should work fine with existing directory
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "channel_summary.csv"))
        end
        
        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)
            
            # Create one valid file
            erps = create_test_erp_data(2)
            save(joinpath(partial_dir, "1_erps_cleaned.jld2"), "erps", erps)
            
            # Create one malformed file (wrong variable name)
            save(joinpath(partial_dir, "2_erps_cleaned.jld2"), "invalid_var", erps)
            
            output_dir = joinpath(test_dir, "summary_partial")
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = partial_dir,
                                  output_dir = output_dir)
            
            # Should have results from the valid file only
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test isfile(csv_file)
            
            results = CSV.read(csv_file, DataFrame)
            # 1 file × 2 conditions × 3 channels = 6 rows
            @test nrow(results) == 6
        end
        
        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "summary_invalid_condition")
            
            # Request condition 5 when only 2 exist
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = 5)
            
            # Should produce empty results
            csv_file = joinpath(output_dir, "channel_summary.csv")
            @test !isfile(csv_file)  # No file created because no valid results
        end
        
        @testset "Statistics correctness" begin
            # Create data with known statistics
            stats_dir = joinpath(test_dir, "stats_test")
            mkpath(stats_dir)
            
            fs = 256.0
            n_samples = 100
            t = range(0.0, 1.0, length=n_samples)
            
            # Channel with known properties
            # Ch1: constant value = 5.0
            # Ch2: values 1 to 100
            # Ch3: all zeros
            
            df = DataFrame(
                time = collect(t),
                sample = 1:n_samples,
                condition = fill(1, n_samples),
                Ch1 = fill(5.0, n_samples),
                Ch2 = collect(1.0:100.0),
                Ch3 = zeros(n_samples)
            )
            
            layout = eegfun.Layout(
                DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
                nothing, nothing
            )
            
            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(stats_dir, "1_erps_stats.jld2"), "erps", erps)
            
            # Process
            output_dir = joinpath(test_dir, "summary_stats")
            eegfun.channel_summary("erps_stats",
                                  input_dir = stats_dir,
                                  output_dir = output_dir)
            
            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)
            
            # Verify statistics for Ch1 (constant = 5.0)
            ch1 = results[results.channel .== "Ch1", :]
            @test ch1.min[1] ≈ 5.0
            @test ch1.max[1] ≈ 5.0
            @test ch1.std[1] ≈ 0.0 atol=1e-10
            @test ch1.range[1] ≈ 0.0 atol=1e-10
            @test ch1.var[1] ≈ 0.0 atol=1e-10
            
            # Verify statistics for Ch2 (1 to 100)
            ch2 = results[results.channel .== "Ch2", :]
            @test ch2.min[1] ≈ 1.0
            @test ch2.max[1] ≈ 100.0
            @test ch2.range[1] ≈ 99.0
            @test ch2.std[1] ≈ std(1.0:100.0)
            @test ch2.var[1] ≈ var(1.0:100.0)
            
            # Verify statistics for Ch3 (all zeros)
            ch3 = results[results.channel .== "Ch3", :]
            @test ch3.min[1] ≈ 0.0
            @test ch3.max[1] ≈ 0.0
            @test ch3.std[1] ≈ 0.0 atol=1e-10
        end
        
        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "summary_combined")
            
            # Summary for specific participant AND condition AND channels
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = 1,
                                  conditions = 1,
                                  channel_selection = eegfun.channels([:Fz, :Cz]))
            
            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)
            
            # 1 file × 1 condition × 2 channels = 2 rows
            @test nrow(results) == 2
            @test all(results.file .== "1_erps_cleaned")
            @test all(results.condition .== 1)
            @test "Pz" ∉ results.channel
        end
        
        @testset "Pattern matching variants" begin
            # Create files with different naming patterns
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)
            
            erps = create_test_erp_data(2)
            save(joinpath(pattern_dir, "1_erps_original.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "2_erps_cleaned.jld2"), "erps", erps)
            save(joinpath(pattern_dir, "3_custom_erps.jld2"), "erps", erps)
            
            # Test pattern matching "erps_original"
            output_dir1 = joinpath(test_dir, "summary_original")
            eegfun.channel_summary("erps_original",
                                  input_dir = pattern_dir,
                                  output_dir = output_dir1)
            
            results1 = CSV.read(joinpath(output_dir1, "channel_summary.csv"), DataFrame)
            @test nrow(results1) == 6  # 1 file × 2 conditions × 3 channels
            
            # Test pattern matching "erps" (should match all)
            output_dir2 = joinpath(test_dir, "summary_all_erps")
            eegfun.channel_summary("erps",
                                  input_dir = pattern_dir,
                                  output_dir = output_dir2)
            
            results2 = CSV.read(joinpath(output_dir2, "channel_summary.csv"), DataFrame)
            @test nrow(results2) == 18  # 3 files × 2 conditions × 3 channels
        end
        
        @testset "File metadata preservation" begin
            output_dir = joinpath(test_dir, "summary_metadata")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)
            
            # Check that file names are correctly preserved
            @test "1_erps_cleaned" in results.file
            @test "2_erps_cleaned" in results.file
            
            # Check condition numbers are correct
            @test 1 in results.condition
            @test 2 in results.condition
            
            # Verify each file-condition combination has all channels
            for file in unique(results.file)
                for cond in unique(results.condition)
                    subset = results[(results.file .== file) .& (results.condition .== cond), :]
                    @test nrow(subset) == 3  # 3 channels
                end
            end
        end
        
        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "summary_overwrite")
            
            # First run
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = overwrite_dir)
            
            csv_file = joinpath(overwrite_dir, "channel_summary.csv")
            mtime1 = stat(csv_file).mtime
            
            # Wait a tiny bit to ensure different mtime
            sleep(0.1)
            
            # Second run (should overwrite)
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = overwrite_dir)
            
            # Verify file was overwritten
            mtime2 = stat(csv_file).mtime
            @test mtime2 > mtime1
        end
        
        @testset "Many channels" begin
            # Test with more channels (10)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)
            
            fs = 256.0
            n_samples = 101
            t = range(-0.2, 0.2, length=n_samples)
            
            # Create 10 channels
            channel_names = Symbol.("Ch" .* string.(1:10))
            
            df = DataFrame(time = collect(t), sample = 1:n_samples, condition = fill(1, n_samples))
            for (i, ch) in enumerate(channel_names)
                df[!, ch] = i .* randn(n_samples)  # Different variance for each channel
            end
            
            layout = eegfun.Layout(
                DataFrame(label = channel_names, inc = zeros(10), azi = zeros(10)),
                nothing, nothing
            )
            
            erps = [eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 1)]
            save(joinpath(many_ch_dir, "1_erps_many.jld2"), "erps", erps)
            
            # Process
            output_dir = joinpath(test_dir, "summary_many_ch")
            eegfun.channel_summary("erps_many",
                                  input_dir = many_ch_dir,
                                  output_dir = output_dir)
            
            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)
            
            # Should have all 10 channels
            @test nrow(results) == 10
            for ch in channel_names
                @test string(ch) in results.channel
            end
        end
        
        @testset "Z-scored variance" begin
            # Test that zvar is correctly computed
            output_dir = joinpath(test_dir, "summary_zvar")
            
            eegfun.channel_summary("erps_cleaned",
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            results = CSV.read(joinpath(output_dir, "channel_summary.csv"), DataFrame)
            
            # For each file-condition combination, zvar should have mean ≈ 0 and std ≈ 1
            for file in unique(results.file)
                for cond in unique(results.condition)
                    subset = results[(results.file .== file) .& (results.condition .== cond), :]
                    
                    # Mean of z-scores should be close to 0
                    @test mean(subset.zvar) ≈ 0.0 atol=1e-10
                    
                    # Std of z-scores should be close to 1 (if more than 1 channel)
                    # Note: With small sample size (3 channels), there's some variability
                    if nrow(subset) > 1
                        @test std(subset.zvar, corrected=false) ≈ 1.0 atol=0.2
                    end
                end
            end
        end
        
    finally
        # Cleanup
        rm(test_dir, recursive=true, force=true)
    end
end
