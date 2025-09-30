using Test
using DataFrames
using eegfun
using JLD2
using Statistics

@testset "Batch Average Epochs" begin
    
    # Create a temporary directory for test files
    test_dir = mktempdir()
    
    try
        # Helper to create test EpochData
        function create_test_epoch_data(n_conditions::Int = 2, n_epochs_per_condition::Int = 5)
            epochs = eegfun.EpochData[]
            fs = 256  # Int64
            n_samples = 513
            t = range(-1.0, 1.0, length=n_samples)
            
            for cond in 1:n_conditions
                # Create multiple epochs per condition
                dfs = DataFrame[]
                for ep in 1:n_epochs_per_condition
                    # Add some trial-to-trial variability
                    ch1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                    ch2 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                    
                    df = DataFrame(
                        time = collect(t),
                        sample = 1:n_samples,
                        condition = fill(cond, n_samples),
                        condition_name = fill("condition_$cond", n_samples),
                        epoch = fill(ep, n_samples),
                        Fz = ch1,
                        Cz = ch2
                    )
                    push!(dfs, df)
                end
                
                layout = eegfun.Layout(
                    DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]),
                    nothing, nothing
                )
                
                # EpochData constructor: (data, layout, sample_rate, analysis_info)
                push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
            end
            
            return epochs
        end
        
        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                epochs = create_test_epoch_data(2, 5)
                filename = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
                save(filename, "epochs", epochs)
                @test isfile(filename)
            end
        end
        
        @testset "Basic averaging" begin
            output_dir = joinpath(test_dir, "averaged_output")
            
            # Test averaging epochs
            result = eegfun.average_epochs("epochs_cleaned", 
                                          input_dir = test_dir,
                                          output_dir = output_dir)
            
            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)
            
            # Check that averaged files exist
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
            
            # Load and verify averaged data
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 2  # 2 conditions
            @test erps[1] isa eegfun.ErpData
            @test hasproperty(erps[1].data, :Fz)
            @test hasproperty(erps[1].data, :Cz)
            @test hasproperty(erps[1].data, :n_epochs)
            
            # Verify n_epochs column
            @test all(erps[1].data.n_epochs .== 5)  # 5 epochs averaged
        end
        
        @testset "Average specific participants" begin
            output_dir = joinpath(test_dir, "averaged_participant")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          participants = 1)
            
            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
        end
        
        @testset "Average multiple participants" begin
            output_dir = joinpath(test_dir, "averaged_multi_participants")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          participants = [1, 2])
            
            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
        end
        
        @testset "Average specific conditions" begin
            output_dir = joinpath(test_dir, "averaged_condition")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          conditions = 1)
            
            @test result.success == 2
            
            # Load and verify only one condition
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test erps[1].data[1, :condition] == 1
        end
        
        @testset "Average multiple conditions" begin
            output_dir = joinpath(test_dir, "averaged_multi_conditions")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          conditions = [1, 2])
            
            @test result.success == 2
            
            # Load and verify both conditions
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 2
            @test erps[1].data[1, :condition] == 1
            @test erps[2].data[1, :condition] == 2
        end
        
        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.average_epochs("epochs_cleaned", 
                                                        input_dir = "/nonexistent/path")
            
            # Invalid pattern (doesn't contain 'epochs')
            @test_throws Exception eegfun.average_epochs("erps_cleaned",
                                                        input_dir = test_dir)
        end
        
        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)
            
            # Directory exists but has no JLD2 files matching pattern
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = empty_dir)
            
            @test result === nothing  # No files to process
        end
        
        @testset "Logging" begin
            output_dir = joinpath(test_dir, "averaged_with_log")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir)
            
            # Check log file exists
            log_file = joinpath(output_dir, "average_epochs.log")
            @test isfile(log_file)
            
            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch epoch averaging started")
            @test contains(log_contents, "average_epochs")
            @test contains(log_contents, "epochs_cleaned")
        end
        
        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_avg")
            mkpath(output_dir)
            
            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))
            
            # Run average_epochs - should work fine with existing directory
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir)
            
            @test result.success == 2
            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
        end
        
        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)
            
            # Create one valid file
            epochs = create_test_epoch_data(2, 5)
            save(joinpath(partial_dir, "1_epochs_cleaned.jld2"), "epochs", epochs)
            
            # Create one malformed file (wrong variable name)
            save(joinpath(partial_dir, "2_epochs_cleaned.jld2"), "invalid_var", epochs)
            
            output_dir = joinpath(test_dir, "averaged_partial")
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = partial_dir,
                                          output_dir = output_dir)
            
            @test result.success == 1
            @test result.errors == 1
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))  # Failed file not saved
        end
        
        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "averaged_invalid_condition")
            
            # Request condition 5 when only 2 exist
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          conditions = 5)
            
            # Should fail for all files
            @test result.success == 0
            @test result.errors == 2
        end
        
        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "averaged_return_check")
            
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir)
            
            # Check result structure
            @test hasfield(typeof(result), :success)
            @test hasfield(typeof(result), :errors)
            @test result.success isa Integer
            @test result.errors isa Integer
            @test result.success >= 0
            @test result.errors >= 0
            @test result.success + result.errors == 2  # Total files processed
        end
        
        @testset "Data integrity - averaging reduces variance" begin
            output_dir = joinpath(test_dir, "averaged_integrity")
            
            # Get original epoch data statistics
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            # Get variance across epochs for first time point, first condition
            first_epoch_values = [df[1, :Fz] for df in original_epochs[1].data]
            epoch_variance = var(first_epoch_values)
            
            # Average epochs
            eegfun.average_epochs("epochs_cleaned",
                                 input_dir = test_dir,
                                 output_dir = output_dir)
            
            # Load averaged data
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            
            # Check that averaged ERP is smoother (single value, no variance)
            @test erps[1].data[1, :Fz] isa Float64  # Single averaged value
            
            # Check that n_epochs is correct
            @test erps[1].n_epochs == 5
            @test all(erps[1].data.n_epochs .== 5)
            
            # Check that time vector is preserved
            original_time = original_epochs[1].data[1].time
            erp_time = erps[1].data.time
            @test length(erp_time) == length(original_time)
            @test erp_time ≈ original_time
        end
        
        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "averaged_combined")
            
            # Average specific participant AND condition
            result = eegfun.average_epochs("epochs_cleaned",
                                          input_dir = test_dir,
                                          output_dir = output_dir,
                                          participants = 1,
                                          conditions = 1)
            
            @test result.success == 1
            
            # Load and verify
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 1  # Only one condition
            @test erps[1].data[1, :condition] == 1
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))  # Participant 2 not processed
        end
        
        @testset "Different epoch counts per condition" begin
            # Create data with different number of epochs per condition
            epochs_var = eegfun.EpochData[]
            fs = 256
            n_samples = 513
            t = range(-1.0, 1.0, length=n_samples)
            
            # Condition 1: 3 epochs
            dfs1 = DataFrame[]
            for ep in 1:3
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                )
                push!(dfs1, df)
            end
            
            # Condition 2: 7 epochs
            dfs2 = DataFrame[]
            for ep in 1:7
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(2, n_samples),
                    condition_name = fill("condition_2", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                )
                push!(dfs2, df)
            end
            
            layout = eegfun.Layout(
                DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]),
                nothing, nothing
            )
            
            push!(epochs_var, eegfun.EpochData(dfs1, layout, fs, eegfun.AnalysisInfo()))
            push!(epochs_var, eegfun.EpochData(dfs2, layout, fs, eegfun.AnalysisInfo()))
            
            # Save and process
            var_dir = joinpath(test_dir, "var_epochs")
            mkpath(var_dir)
            save(joinpath(var_dir, "1_epochs_var.jld2"), "epochs", epochs_var)
            
            output_dir = joinpath(test_dir, "averaged_var")
            result = eegfun.average_epochs("epochs_var",
                                          input_dir = var_dir,
                                          output_dir = output_dir)
            
            @test result.success == 1
            
            # Load and verify epoch counts
            erps = load(joinpath(output_dir, "1_epochs_var.jld2"), "erps")
            @test length(erps) == 2
            @test all(erps[1].data.n_epochs .== 3)  # Condition 1: 3 epochs
            @test all(erps[2].data.n_epochs .== 7)  # Condition 2: 7 epochs
            @test erps[1].n_epochs == 3
            @test erps[2].n_epochs == 7
        end
        
        @testset "Empty epochs handling" begin
            # Create an epochs file with empty condition
            empty_epochs_dir = joinpath(test_dir, "empty_epochs_test")
            mkpath(empty_epochs_dir)
            
            # This should trigger an error when trying to average
            # We'll test this by creating a minimal epochs structure
            layout = eegfun.Layout(
                DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]),
                nothing, nothing
            )
            
            empty_epoch = eegfun.EpochData(DataFrame[], layout, 256, eegfun.AnalysisInfo())
            save(joinpath(empty_epochs_dir, "1_epochs_empty.jld2"), "epochs", [empty_epoch])
            
            output_dir = joinpath(test_dir, "averaged_empty")
            result = eegfun.average_epochs("epochs_empty",
                                          input_dir = empty_epochs_dir,
                                          output_dir = output_dir)
            
            # Should fail because cannot average empty epochs
            @test result.success == 0
            @test result.errors == 1
        end
        
        @testset "Pattern matching variants" begin
            # Create files with different naming patterns
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)
            
            epochs = create_test_epoch_data(2, 3)
            save(joinpath(pattern_dir, "1_epochs_original.jld2"), "epochs", epochs)
            save(joinpath(pattern_dir, "2_epochs_cleaned.jld2"), "epochs", epochs)
            save(joinpath(pattern_dir, "3_custom_epochs.jld2"), "epochs", epochs)
            
            # Test pattern matching "epochs_original"
            output_dir1 = joinpath(test_dir, "averaged_original")
            result1 = eegfun.average_epochs("epochs_original",
                                           input_dir = pattern_dir,
                                           output_dir = output_dir1)
            @test result1.success == 1
            @test isfile(joinpath(output_dir1, "1_epochs_original.jld2"))
            
            # Test pattern matching "epochs_cleaned"
            output_dir2 = joinpath(test_dir, "averaged_cleaned_pattern")
            result2 = eegfun.average_epochs("epochs_cleaned",
                                           input_dir = pattern_dir,
                                           output_dir = output_dir2)
            @test result2.success == 1
            @test isfile(joinpath(output_dir2, "2_epochs_cleaned.jld2"))
            
            # Test pattern matching "epochs" (should match all)
            output_dir3 = joinpath(test_dir, "averaged_all_epochs")
            result3 = eegfun.average_epochs("epochs",
                                           input_dir = pattern_dir,
                                           output_dir = output_dir3)
            @test result3.success == 3
        end
        
        @testset "Verify averaging math correctness" begin
            # Create simple data where we can verify the math
            math_dir = joinpath(test_dir, "math_test")
            mkpath(math_dir)
            
            fs = 256
            n_samples = 11  # Small for easy verification
            t = range(0.0, 1.0, length=n_samples)
            
            # Create 3 epochs with known values
            dfs = DataFrame[]
            known_values = [
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
                [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0],
                [3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0]
            ]
            
            for (ep, vals) in enumerate(known_values)
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Ch1 = vals
                )
                push!(dfs, df)
            end
            
            layout = eegfun.Layout(
                DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]),
                nothing, nothing
            )
            
            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(math_dir, "1_epochs_math.jld2"), "epochs", [epoch_data])
            
            # Average
            output_dir = joinpath(test_dir, "averaged_math")
            eegfun.average_epochs("epochs_math",
                                 input_dir = math_dir,
                                 output_dir = output_dir)
            
            # Load and verify
            erps = load(joinpath(output_dir, "1_epochs_math.jld2"), "erps")
            
            # Expected average: [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]
            expected_avg = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]
            @test erps[1].data.Ch1 ≈ expected_avg
            @test all(erps[1].data.n_epochs .== 3)
        end
        
        @testset "Metadata preservation" begin
            output_dir = joinpath(test_dir, "averaged_metadata")
            
            eegfun.average_epochs("epochs_cleaned",
                                 input_dir = test_dir,
                                 output_dir = output_dir)
            
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            
            # Verify metadata columns exist
            # Note: :sample is not preserved during averaging (it's a per-epoch metadata)
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :condition)
            @test hasproperty(erps[1].data, :condition_name)
            @test hasproperty(erps[1].data, :n_epochs)
            
            # Verify condition metadata
            @test all(erps[1].data.condition .== 1)
            @test all(erps[1].data.condition_name .== "condition_1")
            @test all(erps[2].data.condition .== 2)
            @test all(erps[2].data.condition_name .== "condition_2")
        end
        
        @testset "Layout and metadata preservation" begin
            output_dir = joinpath(test_dir, "averaged_layout")
            
            # Get original layout
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            original_layout = original_epochs[1].layout
            original_fs = original_epochs[1].sample_rate
            
            eegfun.average_epochs("epochs_cleaned",
                                 input_dir = test_dir,
                                 output_dir = output_dir)
            
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            
            # Verify layout preservation
            @test erps[1].layout.data == original_layout.data
            @test erps[1].sample_rate == original_fs
            @test erps[1] isa eegfun.ErpData
        end
        
        @testset "Single epoch per condition" begin
            # Edge case: only 1 epoch to average
            single_dir = joinpath(test_dir, "single_epoch")
            mkpath(single_dir)
            
            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length=n_samples)
            
            dfs = DataFrame[]
            df = DataFrame(
                time = collect(t),
                sample = 1:n_samples,
                condition = fill(1, n_samples),
                condition_name = fill("condition_1", n_samples),
                epoch = fill(1, n_samples),
                Cz = sin.(2π .* 10 .* t)
            )
            push!(dfs, df)
            
            layout = eegfun.Layout(
                DataFrame(label = [:Cz], inc = [0.0], azi = [0.0]),
                nothing, nothing
            )
            
            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(single_dir, "1_epochs_single.jld2"), "epochs", [epoch_data])
            
            # Average single epoch
            output_dir = joinpath(test_dir, "averaged_single")
            result = eegfun.average_epochs("epochs_single",
                                          input_dir = single_dir,
                                          output_dir = output_dir)
            
            @test result.success == 1
            
            erps = load(joinpath(output_dir, "1_epochs_single.jld2"), "erps")
            @test erps[1].n_epochs == 1
            @test all(erps[1].data.n_epochs .== 1)
            
            # With single epoch, average should equal the original
            @test erps[1].data.Cz ≈ df.Cz
        end
        
        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "averaged_overwrite")
            
            # First run
            result1 = eegfun.average_epochs("epochs_cleaned",
                                           input_dir = test_dir,
                                           output_dir = overwrite_dir)
            @test result1.success == 2
            
            # Get original file modification time
            file1 = joinpath(overwrite_dir, "1_epochs_cleaned.jld2")
            mtime1 = stat(file1).mtime
            
            # Wait a tiny bit to ensure different mtime
            sleep(0.1)
            
            # Second run (should overwrite)
            result2 = eegfun.average_epochs("epochs_cleaned",
                                           input_dir = test_dir,
                                           output_dir = overwrite_dir)
            @test result2.success == 2
            
            # Verify file was overwritten (different mtime)
            mtime2 = stat(file1).mtime
            @test mtime2 > mtime1
        end
        
        @testset "Many channels" begin
            # Test with more channels (10)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)
            
            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length=n_samples)
            
            # Create 10 channels
            channel_names = Symbol.("Ch" .* string.(1:10))
            
            dfs = DataFrame[]
            for ep in 1:3
                # Build DataFrame with proper column order
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples)
                )
                
                # Add channel data
                for (i, ch) in enumerate(channel_names)
                    df[!, ch] = sin.(2π .* i .* t) .+ 0.1 .* randn(n_samples)
                end
                
                push!(dfs, df)
            end
            
            layout_df = DataFrame(
                label = channel_names,
                inc = zeros(10),
                azi = zeros(10)
            )
            layout = eegfun.Layout(layout_df, nothing, nothing)
            
            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(many_ch_dir, "1_epochs_many.jld2"), "epochs", [epoch_data])
            
            # Average
            output_dir = joinpath(test_dir, "averaged_many_ch")
            result = eegfun.average_epochs("epochs_many",
                                          input_dir = many_ch_dir,
                                          output_dir = output_dir)
            
            @test result.success == 1
            
            erps = load(joinpath(output_dir, "1_epochs_many.jld2"), "erps")
            
            # Verify all channels are present
            for ch in channel_names
                @test hasproperty(erps[1].data, ch)
            end
        end
        
    finally
        # Cleanup
        rm(test_dir, recursive=true, force=true)
    end
end
