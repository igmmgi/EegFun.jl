"""
Test suite for src/analysis/batch/difference_conditions.jl
"""

using Test
using JLD2
using DataFrames
using CSV

# Helper function to create test ERP data
function create_test_erp_data(participant::Int, condition::Int, n_timepoints::Int=100, n_channels::Int=3)
    # Create time vector
    time = collect(range(-0.2, 0.8, length=n_timepoints))
    
    # Create channel data with some condition-specific patterns
    channel_data = Dict{Symbol, Any}()
    channel_data[:time] = time
    
    # Add metadata columns
    channel_data[:condition] = fill(condition, n_timepoints)
    channel_data[:condition_name] = fill("condition_$condition", n_timepoints)
    channel_data[:participant] = fill(participant, n_timepoints)
    
    # Add EEG channels with condition-specific patterns
    for (i, ch) in enumerate([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
        # Create some signal with condition-specific amplitude
        signal = sin.(2π * 0.1 * time) .* (condition * 0.5) .+ randn(n_timepoints) * 0.1
        channel_data[ch] = signal
    end
    
    # Create DataFrame with columns in correct order
    df = DataFrame()
    df.time = channel_data[:time]
    df.condition = channel_data[:condition]
    df.condition_name = channel_data[:condition_name]
    df.participant = channel_data[:participant]
    for ch in [:Fz, :Cz, :Pz][1:min(n_channels, 3)]
        df[!, ch] = channel_data[ch]
    end
    
    # Create ErpData
    layout = eegfun.ElectrodeLayout([:Fz, :Cz, :Pz][1:min(n_channels, 3)])
    return eegfun.ErpData(df, layout, 250.0, "test_analysis", 10)
end

@testset "Batch Difference Conditions" begin
    # Create temporary test directory
    test_dir = mktempdir()
    
    @testset "Basic difference wave creation" begin
        # Create test ERP files
        for participant in 1:3
            erps = [
                create_test_erp_data(participant, 1),
                create_test_erp_data(participant, 2),
                create_test_erp_data(participant, 3),
                create_test_erp_data(participant, 4)
            ]
            
            file_path = joinpath(test_dir, "$(participant)_erps_cleaned.jld2")
            save(file_path, "erps", erps)
        end
        
        output_dir = joinpath(test_dir, "differences")
        
        # Test basic difference creation
        result = eegfun.difference_conditions("erps_cleaned", [(1, 2), (3, 4)],
                                            input_dir = test_dir,
                                            output_dir = output_dir)
        
        # Verify output files were created
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 3  # One for each participant
        
        # Load and verify one output file
        output_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        @test isfile(output_file)
        
        differences = load(output_file, "differences")
        @test length(differences) == 2  # Two difference waves
        
        # Verify difference wave structure
        diff1 = differences[1]
        @test diff1.data.condition[1] == 1  # First difference labeled as condition 1
        @test diff1.data.condition_name[1] == "difference_1_2"
        
        diff2 = differences[2]
        @test diff2.data.condition[1] == 2  # Second difference labeled as condition 2
        @test diff2.data.condition_name[1] == "difference_3_4"
    end
    
    @testset "Single participant processing" begin
        output_dir = joinpath(test_dir, "differences_single")
        
        result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                            input_dir = test_dir,
                                            participants = 2,
                                            output_dir = output_dir)
        
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 1  # Only participant 2
        @test "2_erps_cleaned.jld2" in output_files
    end
    
    @testset "Vector condition pairs" begin
        output_dir = joinpath(test_dir, "differences_vector")
        
        result = eegfun.difference_conditions("erps_cleaned", [[1, 2], [3, 4]],
                                            input_dir = test_dir,
                                            output_dir = output_dir)
        
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 3
    end
    
    @testset "Missing conditions handling" begin
        # Create file with only some conditions
        erps = [
            create_test_erp_data(99, 1),
            create_test_erp_data(99, 2)
            # Missing conditions 3 and 4
        ]
        
        file_path = joinpath(test_dir, "99_erps_cleaned.jld2")
        save(file_path, "erps", erps)
        
        output_dir = joinpath(test_dir, "differences_missing")
        
        result = eegfun.difference_conditions("erps_cleaned", [(1, 2), (3, 4)],
                                            input_dir = test_dir,
                                            participants = 99,
                                            output_dir = output_dir)
        
        # Should still create file but only with available pairs
        @test isdir(output_dir)
        output_files = readdir(output_dir)
        @test length(output_files) == 1
        
        differences = load(joinpath(output_dir, "99_erps_cleaned.jld2"), "differences")
        @test length(differences) == 1  # Only one difference wave (1-2)
    end
    
    @testset "Error handling" begin
        @testset "Invalid input directory" begin
            @test_throws Exception eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                                               input_dir = "/nonexistent/dir")
        end
        
        @testset "Non-ERP pattern" begin
            @test_throws Exception eegfun.difference_conditions("epochs_cleaned", [(1, 2)],
                                                               input_dir = test_dir)
        end
        
        @testset "Empty condition pairs" begin
            @test_throws Exception eegfun.difference_conditions("erps_cleaned", [],
                                                               input_dir = test_dir)
        end
        
        @testset "Invalid condition pairs" begin
            @test_throws Exception eegfun.difference_conditions("erps_cleaned", [(1, "invalid")],
                                                               input_dir = test_dir)
        end
    end
    
    @testset "Data integrity" begin
        output_dir = joinpath(test_dir, "differences_integrity")
        
        result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                            input_dir = test_dir,
                                            participants = 1,
                                            output_dir = output_dir)
        
        # Load original and difference data
        original_file = joinpath(test_dir, "1_erps_cleaned.jld2")
        original_erps = load(original_file, "erps")
        
        diff_file = joinpath(output_dir, "1_erps_cleaned.jld2")
        differences = load(diff_file, "differences")
        
        # Verify difference calculation
        erp1 = original_erps[1]  # Condition 1
        erp2 = original_erps[2]  # Condition 2
        diff = differences[1]    # Difference 1-2
        
        # Check that difference is actually erp1 - erp2
        for ch in [:Fz, :Cz, :Pz]
            if hasproperty(erp1.data, ch) && hasproperty(erp2.data, ch) && hasproperty(diff.data, ch)
                expected_diff = erp1.data[!, ch] .- erp2.data[!, ch]
                @test all(abs.(diff.data[!, ch] .- expected_diff) .< 1e-10)
            end
        end
        
        # Verify metadata preservation
        @test diff.sample_rate == erp1.sample_rate
        @test nrow(diff.data) == nrow(erp1.data)
        @test diff.n_epochs == min(erp1.n_epochs, erp2.n_epochs)
    end
    
    @testset "Edge cases" begin
        @testset "Identical condition pairs" begin
            output_dir = joinpath(test_dir, "differences_identical")
            
            # Should work but create zero differences
            result = eegfun.difference_conditions("erps_cleaned", [(1, 1)],
                                                input_dir = test_dir,
                                                participants = 1,
                                                output_dir = output_dir)
            
            @test isdir(output_dir)
            differences = load(joinpath(output_dir, "1_erps_cleaned.jld2"), "differences")
            @test length(differences) == 1
            
            # Verify difference is zero
            diff = differences[1]
            for ch in [:Fz, :Cz, :Pz]
                if hasproperty(diff.data, ch)
                    @test all(abs.(diff.data[!, ch]) .< 1e-10)
                end
            end
        end
        
        @testset "No matching files" begin
            output_dir = joinpath(test_dir, "differences_none")
            
            result = eegfun.difference_conditions("nonexistent_pattern", [(1, 2)],
                                                input_dir = test_dir,
                                                output_dir = output_dir)
            
            @test result === nothing
        end
        
        @testset "Empty ERP data" begin
            # Create file with empty ERP list
            empty_file = joinpath(test_dir, "empty_erps_cleaned.jld2")
            save(empty_file, "erps", eegfun.ErpData[])
            
            output_dir = joinpath(test_dir, "differences_empty")
            
            result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                                input_dir = test_dir,
                                                participants = 999,  # Non-existent participant
                                                output_dir = output_dir)
            
            @test result === nothing
        end
    end
    
    @testset "Output directory handling" begin
        @testset "Custom output directory" begin
            custom_dir = joinpath(test_dir, "custom_output")
            
            result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                                input_dir = test_dir,
                                                output_dir = custom_dir)
            
            @test isdir(custom_dir)
            @test length(readdir(custom_dir)) == 3
        end
        
        @testset "Auto-generated output directory" begin
            result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                                input_dir = test_dir)
            
            # Should create directory with pattern-based name
            expected_dir = joinpath(test_dir, "differences_erps_cleaned_1-2")
            @test isdir(expected_dir)
        end
    end
    
    @testset "Logging and return values" begin
        output_dir = joinpath(test_dir, "differences_logging")
        
        result = eegfun.difference_conditions("erps_cleaned", [(1, 2)],
                                            input_dir = test_dir,
                                            output_dir = output_dir)
        
        # Check that log file was created
        log_file = joinpath(output_dir, "difference_conditions.log")
        @test isfile(log_file)
        
        # Verify log content contains expected information
        log_content = read(log_file, String)
        @test occursin("Batch condition differencing started", log_content)
        @test occursin("Found 3 JLD2 files", log_content)
        @test occursin("Batch operation complete", log_content)
    end
    
    # Cleanup
    rm(test_dir, recursive=true)
end
