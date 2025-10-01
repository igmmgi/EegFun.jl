using Test
using DataFrames
using eegfun
using JLD2
using Statistics

@testset "Batch Filter" begin
    
    # Create a temporary directory for test files
    test_dir = mktempdir()
    
    try
        # Helper to create test ERP data
        function create_test_erp_data(n_conditions::Int = 2)
            erps = eegfun.ErpData[]
            fs = 256.0
            n_samples = 513
            t = range(-1.0, 1.0, length=n_samples)
            
            for cond in 1:n_conditions
                # Create signal with low and high frequency components
                ch1 = sin.(2π .* 5 .* t) .+ 0.2 .* sin.(2π .* 50 .* t) .+ 0.1 .* randn(n_samples)
                ch2 = cos.(2π .* 5 .* t) .+ 0.2 .* cos.(2π .* 50 .* t) .+ 0.1 .* randn(n_samples)
                
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(cond, n_samples),
                    Fz = ch1,
                    Cz = ch2
                )
                
                layout = eegfun.Layout(
                    DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]),
                    nothing, nothing
                )
                
                push!(erps, eegfun.ErpData(df, layout, fs, eegfun.AnalysisInfo(), 50))
            end
            
            return erps
        end
        
        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                erps = create_test_erp_data(2)
                # Use filename format consistent with codebase (numeric participant ID)
                filename = joinpath(test_dir, "$(participant)_erps.jld2")
                save(filename, "erps", erps)
                @test isfile(filename)
            end
        end
        
        @testset "Basic filtering" begin
            output_dir = joinpath(test_dir, "filtered_output")
            
            # Test low-pass filtering
            result = eegfun.filter("erps", 30.0, 
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  filter_type = "lp")
            
            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)
            
            # Check that filtered files exist
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test isfile(joinpath(output_dir, "2_erps.jld2"))
            
            # Load and verify filtered data
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            @test length(filtered_data) == 2  # 2 conditions
            @test filtered_data[1] isa eegfun.ErpData
            @test hasproperty(filtered_data[1].data, :Fz)
            @test hasproperty(filtered_data[1].data, :Cz)
            
            # Verify high frequencies are attenuated (not a rigorous test, just sanity check)
            @test std(filtered_data[1].data.Fz) < 2.0  # Should be reduced
        end
        
        @testset "Filter specific participants" begin
            output_dir = joinpath(test_dir, "filtered_participant")
            
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = 1)
            
            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))
        end
        
        @testset "Filter specific conditions" begin
            output_dir = joinpath(test_dir, "filtered_condition")
            
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = 1)
            
            @test result.success == 2
            
            # Load and verify only one condition
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            @test length(filtered_data) == 1
            @test filtered_data[1].data[1, :condition] == 1
        end
        
        @testset "High-pass filter" begin
            output_dir = joinpath(test_dir, "filtered_hp")
            
            result = eegfun.filter("erps", 1.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  filter_type = "hp")
            
            @test result.success == 2
            @test result.errors == 0
            
            # Load and verify
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            @test filtered_data[1] isa eegfun.ErpData
        end
        
        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.filter("erps", 30.0, 
                                                input_dir = "/nonexistent/path")
            
            # Invalid filter type
            @test_throws Exception eegfun.filter("erps", 30.0,
                                                input_dir = test_dir,
                                                filter_type = "invalid")
            
            # Invalid cutoff frequency
            @test_throws Exception eegfun.filter("erps", -10.0,
                                                input_dir = test_dir)
        end
        
        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "filtered_nomatch")
            
            # Pattern that won't match any files
            result = eegfun.filter("nonexistent_pattern", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            @test result === nothing  # Function returns nothing when no files found
        end
        
        @testset "Logging" begin
            output_dir = joinpath(test_dir, "filtered_with_log")
            
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            # Check log file exists
            log_file = joinpath(output_dir, "filter.log")
            @test isfile(log_file)
            
            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch filtering started")
            @test contains(log_contents, "filter_type=\"lp\"")
            @test contains(log_contents, "cutoff: 30.0 Hz")
        end
        
        @testset "Multiple participants filter" begin
            output_dir = joinpath(test_dir, "filtered_multi_participants")
            
            # Filter both participants
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = [1, 2])
            
            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test isfile(joinpath(output_dir, "2_erps.jld2"))
        end
        
        @testset "Multiple conditions filter" begin
            output_dir = joinpath(test_dir, "filtered_multi_conditions")
            
            # Filter both conditions
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = [1, 2])
            
            @test result.success == 2
            
            # Load and verify both conditions are present
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            @test length(filtered_data) == 2
            @test filtered_data[1].data[1, :condition] == 1
            @test filtered_data[2].data[1, :condition] == 2
        end
        
        @testset "EpochData support" begin
            # Create EpochData test files
            epochs_dir = joinpath(test_dir, "epochs_test")
            mkpath(epochs_dir)
            
            function create_test_epoch_data()
                epochs = eegfun.EpochData[]
                fs = 256  # Int64, not Float64
                n_samples = 513
                t = range(-1.0, 1.0, length=n_samples)
                
                for cond in 1:2
                    # Create multiple epochs per condition
                    dfs = DataFrame[]
                    for ep in 1:3
                        ch1 = sin.(2π .* 5 .* t) .+ 0.2 .* sin.(2π .* 50 .* t) .+ 0.1 .* randn(n_samples)
                        ch2 = cos.(2π .* 5 .* t) .+ 0.2 .* cos.(2π .* 50 .* t) .+ 0.1 .* randn(n_samples)
                        
                        df = DataFrame(
                            time = collect(t),
                            sample = 1:n_samples,
                            condition = fill(cond, n_samples),
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
            
            # Save epoch data
            epochs = create_test_epoch_data()
            save(joinpath(epochs_dir, "1_epochs.jld2"), "epochs", epochs)
            
            # Filter epoch data
            output_dir = joinpath(test_dir, "filtered_epochs")
            result = eegfun.filter("epochs", 30.0,
                                  input_dir = epochs_dir,
                                  output_dir = output_dir)
            
            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs.jld2"))
            
            # Load and verify
            filtered_epochs = load(joinpath(output_dir, "1_epochs.jld2"), "epochs")
            @test filtered_epochs[1] isa eegfun.EpochData
            @test length(filtered_epochs[1].data) == 3  # 3 epochs
        end
        
        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output")
            mkpath(output_dir)
            
            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))
            
            # Run filter - should work fine with existing directory
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir)
            
            @test result.success == 2
            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
        end
        
        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)
            
            # Create one valid file
            erps = create_test_erp_data(2)
            save(joinpath(partial_dir, "1_erps.jld2"), "erps", erps)
            
            # Create one malformed file (wrong variable name)
            save(joinpath(partial_dir, "2_erps.jld2"), "invalid_var", erps)
            
            output_dir = joinpath(test_dir, "filtered_partial")
            result = eegfun.filter("erps", 30.0,
                                  input_dir = partial_dir,
                                  output_dir = output_dir)
            
            @test result.success == 1
            @test result.errors == 1
            @test isfile(joinpath(output_dir, "1_erps.jld2"))
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))  # Failed file not saved
        end
        
        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "filtered_invalid_condition")
            
            # Request condition 5 when only 2 exist
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  conditions = 5)
            
            # Should fail for all files
            @test result.success == 0
            @test result.errors == 2
        end
        
        @testset "Empty pattern match" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)
            
            # Directory exists but has no JLD2 files
            result = eegfun.filter("erps", 30.0,
                                  input_dir = empty_dir)
            
            @test result === nothing  # No files to process
        end
        
        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "filtered_return_check")
            
            result = eegfun.filter("erps", 30.0,
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
        
        @testset "Data integrity - frequency attenuation" begin
            output_dir = joinpath(test_dir, "filtered_integrity")
            
            # Get original data statistics
            original_data = load(joinpath(test_dir, "1_erps.jld2"), "erps")
            original_signal = original_data[1].data.Fz
            
            # Apply low-pass filter
            eegfun.filter("erps", 30.0,
                         input_dir = test_dir,
                         output_dir = output_dir)
            
            # Load filtered data
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            filtered_signal = filtered_data[1].data.Fz
            
            # Check that high-frequency noise is reduced (lower std deviation)
            @test std(filtered_signal) < std(original_signal)
            
            # Check that data length is preserved
            @test length(filtered_signal) == length(original_signal)
            
            # Check that time vector is preserved
            @test filtered_data[1].data.time == original_data[1].data.time
        end
        
        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "filtered_combined")
            
            # Filter specific participant AND condition
            result = eegfun.filter("erps", 30.0,
                                  input_dir = test_dir,
                                  output_dir = output_dir,
                                  participants = 1,
                                  conditions = 1)
            
            @test result.success == 1
            
            # Load and verify
            filtered_data = load(joinpath(output_dir, "1_erps.jld2"), "erps")
            @test length(filtered_data) == 1  # Only one condition
            @test filtered_data[1].data[1, :condition] == 1
            @test !isfile(joinpath(output_dir, "2_erps.jld2"))  # Participant 2 not processed
        end
        
    finally
        # Cleanup
        rm(test_dir, recursive=true, force=true)
    end
end
