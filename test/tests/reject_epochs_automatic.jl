"""
Test suite for src/analysis/reject_epochs_automatic.jl
"""

using Test
using JLD2
using DataFrames
using Statistics
using Random

# Helper function to create test epoched data with varying artifact levels
function create_test_epochs_with_artifacts(
    n_epochs::Int = 20,
    n_timepoints::Int = 100,
    n_channels::Int = 3;
    n_bad_epochs::Int = 3,
    artifact_scale::Float64 = 10.0
)
    time = collect(range(-0.2, 0.8, length = n_timepoints))
    sample_rate = Int(round(n_timepoints / (time[end] - time[1])))
    
    epochs = DataFrame[]
    bad_epoch_indices = sort(randperm(n_epochs)[1:n_bad_epochs])
    
    for i = 1:n_epochs
        epoch_df = DataFrame()
        epoch_df.time = copy(time)
        epoch_df.trial = fill(i, n_timepoints)
        epoch_df.condition = fill(1, n_timepoints)
        
        # Determine if this is a bad epoch
        is_bad = i in bad_epoch_indices
        
        # Add EEG channels
        for ch_idx = 1:n_channels
            ch_name = Symbol("Ch$ch_idx")
            
            if is_bad
                # Add large artifacts to bad epochs
                signal = randn(n_timepoints) .* artifact_scale .+ sin.(2π * 2.0 * time) .* 2.0
            else
                # Normal signal for good epochs
                signal = randn(n_timepoints) .* 0.5 .+ sin.(2π * 2.0 * time) .* 2.0
            end
            
            epoch_df[!, ch_name] = signal
        end
        
        push!(epochs, epoch_df)
    end
    
    # Create layout
    channel_labels = [Symbol("Ch$i") for i = 1:n_channels]
    layout = eegfun.Layout(
        DataFrame(
            label = channel_labels,
            inc = zeros(n_channels),
            azi = collect(range(0, 2π, length = n_channels + 1))[1:n_channels],
        ),
        nothing,
        nothing,
    )
    
    return eegfun.EpochData(epochs, layout, sample_rate, eegfun.AnalysisInfo()), bad_epoch_indices
end

@testset "Automatic Epoch Rejection" begin
    @testset "Basic rejection (non-mutating)" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(20, 100, 3, n_bad_epochs = 3)
        original_n_epochs = length(epoch_data.data)
        
        # Apply rejection
        clean_data, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Check that original data is unchanged
        @test length(epoch_data.data) == original_n_epochs
        
        # Check that some epochs were rejected
        @test length(clean_data.data) < original_n_epochs
        @test rejection_info.n_original == original_n_epochs
        @test rejection_info.n_remaining == length(clean_data.data)
        @test rejection_info.n_remaining == original_n_epochs - length(rejection_info.rejected_epochs)
    end
    
    @testset "In-place rejection" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(20, 100, 3, n_bad_epochs = 3)
        original_n_epochs = length(epoch_data.data)
        
        # Apply rejection in-place
        rejection_info = eegfun.reject_epochs_automatic!(epoch_data, 2.0)
        
        # Check that data was modified
        @test length(epoch_data.data) < original_n_epochs
        @test length(epoch_data.data) == rejection_info.n_remaining
    end
    
    @testset "Rejection with different z-criteria" begin
        epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 3, n_bad_epochs = 4)
        
        # More aggressive (z=1.5) should reject more
        clean_aggressive, info_aggressive = eegfun.reject_epochs_automatic(epoch_data, 1.5)
        
        # More conservative (z=3.0) should reject fewer
        clean_conservative, info_conservative = eegfun.reject_epochs_automatic(epoch_data, 3.0)
        
        # Check relationship
        @test info_aggressive.n_remaining <= info_conservative.n_remaining
        @test length(info_aggressive.rejected_epochs) >= length(info_conservative.rejected_epochs)
    end
    
    @testset "EpochRejectionInfo structure" begin
        epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 3, n_bad_epochs = 3)
        
        _, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Check structure
        @test rejection_info isa eegfun.EpochRejectionInfo
        @test rejection_info.z_criterion == 2.0
        @test rejection_info.n_original > 0
        @test rejection_info.n_remaining >= 0
        @test rejection_info.n_remaining <= rejection_info.n_original
        
        # Check that rejected epochs is the union of all criteria
        all_rejected = unique(vcat(
            rejection_info.rejected_by_variance,
            rejection_info.rejected_by_max,
            rejection_info.rejected_by_min,
            rejection_info.rejected_by_abs,
            rejection_info.rejected_by_range,
            rejection_info.rejected_by_kurtosis
        ))
        @test sort(all_rejected) == sort(rejection_info.rejected_epochs)
    end
    
    @testset "Channel selection" begin
        epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 5, n_bad_epochs = 3)
        
        # Rejection using all channels
        _, info_all = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Rejection using subset of channels
        _, info_subset = eegfun.reject_epochs_automatic(epoch_data, 2.0,
                                                        channel_selection = eegfun.channels([:Ch1, :Ch2]))
        
        # Results may differ depending on which channels have artifacts
        @test info_all.n_original == info_subset.n_original
    end
    
    @testset "Metrics calculation" begin
        # Create epochs with known characteristics
        epoch_data, _ = create_test_epochs_with_artifacts(10, 100, 3, n_bad_epochs = 2)
        
        _, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # At least one metric should have identified epochs for rejection
        n_criteria_triggered = sum([
            !isempty(rejection_info.rejected_by_variance),
            !isempty(rejection_info.rejected_by_max),
            !isempty(rejection_info.rejected_by_min),
            !isempty(rejection_info.rejected_by_abs),
            !isempty(rejection_info.rejected_by_range),
            !isempty(rejection_info.rejected_by_kurtosis)
        ])
        
        if !isempty(rejection_info.rejected_epochs)
            @test n_criteria_triggered > 0
        end
    end
    
    @testset "Data structure preservation" begin
        epoch_data, _ = create_test_epochs_with_artifacts(15, 100, 3, n_bad_epochs = 2)
        
        clean_data, _ = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Check that data structure is preserved
        @test clean_data isa eegfun.EpochData
        @test clean_data.sample_rate == epoch_data.sample_rate
        @test clean_data.layout.data == epoch_data.layout.data
        
        # Check that remaining epochs have correct structure
        for epoch in clean_data.data
            @test hasproperty(epoch, :time)
            @test hasproperty(epoch, :trial)
            @test hasproperty(epoch, :condition)
            @test hasproperty(epoch, :Ch1)
            @test hasproperty(epoch, :Ch2)
            @test hasproperty(epoch, :Ch3)
        end
    end
    
    @testset "Error handling: empty data" begin
        # Create empty epoch data
        empty_epochs = eegfun.EpochData(
            DataFrame[],
            eegfun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing),
            1000,
            eegfun.AnalysisInfo()
        )
        
        @test_throws Exception eegfun.reject_epochs_automatic!(empty_epochs, 2.0)
    end
    
    @testset "Error handling: invalid z-criterion" begin
        epoch_data, _ = create_test_epochs_with_artifacts(10, 100, 3)
        
        # Negative criterion
        @test_throws Exception eegfun.reject_epochs_automatic(epoch_data, -1.0)
        
        # Zero criterion
        @test_throws Exception eegfun.reject_epochs_automatic(epoch_data, 0.0)
    end
    
    @testset "Error handling: no channels selected" begin
        epoch_data, _ = create_test_epochs_with_artifacts(10, 100, 3)
        
        # Empty channel selection
        @test_throws Exception eegfun.reject_epochs_automatic(epoch_data, 2.0,
                                                               channel_selection = eegfun.channels([]))
    end
    
    @testset "Edge case: no epochs rejected" begin
        # Create very clean data where no epochs should be rejected with conservative criterion
        epoch_data, _ = create_test_epochs_with_artifacts(10, 100, 3, n_bad_epochs = 0)
        
        clean_data, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 10.0)  # Very high criterion
        
        # Should keep all epochs (or most of them)
        @test length(clean_data.data) == length(epoch_data.data) || 
              length(clean_data.data) >= length(epoch_data.data) - 1
        @test isempty(rejection_info.rejected_epochs) || length(rejection_info.rejected_epochs) <= 1
    end
    
    @testset "Edge case: all epochs rejected" begin
        # This is unlikely but test that it handles it gracefully
        epoch_data, _ = create_test_epochs_with_artifacts(5, 100, 3, n_bad_epochs = 5, artifact_scale = 100.0)
        
        clean_data, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 0.5)  # Very low criterion
        
        # Should handle having few or no epochs remaining
        @test clean_data.data isa Vector{DataFrame}
        @test rejection_info.n_remaining >= 0
    end
    
    @testset "Realistic scenario: mixed artifacts" begin
        # Create more realistic test data
        n_epochs = 30
        time = collect(range(-0.2, 0.8, length = 150))
        epochs = DataFrame[]
        
        # Add various types of epochs
        for i = 1:n_epochs
            epoch_df = DataFrame()
            epoch_df.time = copy(time)
            epoch_df.trial = fill(i, length(time))
            epoch_df.condition = fill(1, length(time))
            
            if i <= 20
                # Good epochs (most of them)
                epoch_df.Ch1 = randn(length(time)) .* 0.5
                epoch_df.Ch2 = randn(length(time)) .* 0.5
                epoch_df.Ch3 = randn(length(time)) .* 0.5
            elseif i <= 23
                # High variance epochs
                epoch_df.Ch1 = randn(length(time)) .* 5.0
                epoch_df.Ch2 = randn(length(time)) .* 5.0
                epoch_df.Ch3 = randn(length(time)) .* 5.0
            elseif i <= 26
                # High amplitude epochs
                epoch_df.Ch1 = randn(length(time)) .* 0.5 .+ 50.0
                epoch_df.Ch2 = randn(length(time)) .* 0.5 .+ 50.0
                epoch_df.Ch3 = randn(length(time)) .* 0.5
            else
                # High kurtosis epochs (spiky)
                signal = randn(length(time)) .* 0.5
                signal[rand(1:length(time), 5)] .= 20.0  # Add spikes
                epoch_df.Ch1 = signal
                epoch_df.Ch2 = signal .+ randn(length(time)) .* 0.1
                epoch_df.Ch3 = signal .+ randn(length(time)) .* 0.1
            end
            
            push!(epochs, epoch_df)
        end
        
        layout = eegfun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
            nothing,
            nothing
        )
        
        epoch_data = eegfun.EpochData(epochs, layout, 1000, eegfun.AnalysisInfo())
        
        # Apply rejection
        clean_data, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Should reject some but not all epochs
        @test 0 < length(rejection_info.rejected_epochs) < n_epochs
        @test rejection_info.n_remaining > n_epochs / 2  # Keep at least half
        
        # Different criteria should catch different types of artifacts
        @test !isempty(rejection_info.rejected_by_variance) ||
              !isempty(rejection_info.rejected_by_max) ||
              !isempty(rejection_info.rejected_by_kurtosis)
    end
    
    @testset "Show method for EpochRejectionInfo" begin
        epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 3)
        _, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Test that show method works without error
        io = IOBuffer()
        show(io, rejection_info)
        output = String(take!(io))
        
        @test occursin("EpochRejectionInfo", output)
        @test occursin("Z-criterion", output)
        @test occursin("Original epochs", output)
        @test occursin("Remaining epochs", output)
    end
    
    @testset "Rejection report saving" begin
        epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 3)
        _, rejection_info = eegfun.reject_epochs_automatic(epoch_data, 2.0)
        
        # Create temporary file
        test_file = tempname() * ".txt"
        
        # Save report
        eegfun.save_rejection_report(rejection_info, test_file)
        
        # Check file was created and has content
        @test isfile(test_file)
        content = read(test_file, String)
        @test !isempty(content)
        @test occursin("AUTOMATIC EPOCH REJECTION REPORT", content)
        @test occursin("Z-criterion", content)
        
        # Cleanup
        rm(test_file)
    end
end

@testset "Batch Automatic Epoch Rejection" begin
    # Create temporary test directory
    test_dir = mktempdir()
    
    @testset "Basic batch processing" begin
        # Create test epoch files for multiple participants
        for participant = 1:3
            epoch_data, _ = create_test_epochs_with_artifacts(20, 100, 3, n_bad_epochs = 3)
            file_path = joinpath(test_dir, "$(participant)_epochs.jld2")
            save(file_path, "epochs", epoch_data)
        end
        
        output_dir = joinpath(test_dir, "rejected_test")
        
        # Test basic rejection
        result = eegfun.reject_epochs_automatic("epochs", 2.0,
                                                input_dir = test_dir,
                                                output_dir = output_dir)
        
        # Verify output directory was created
        @test isdir(output_dir)
        
        # Verify output files
        output_files = readdir(output_dir)
        @test "1_epochs.jld2" in output_files
        @test "2_epochs.jld2" in output_files
        @test "3_epochs.jld2" in output_files
        
        # Check for rejection reports
        @test "1_epochs_rejection_report.txt" in output_files
        @test "2_epochs_rejection_report.txt" in output_files
        @test "3_epochs_rejection_report.txt" in output_files
        
        # Load and verify one file
        clean_epochs = load(joinpath(output_dir, "1_epochs.jld2"), "epochs")
        @test clean_epochs isa eegfun.EpochData
        
        # Should have rejected some epochs
        original_epochs = load(joinpath(test_dir, "1_epochs.jld2"), "epochs")
        @test length(clean_epochs.data) <= length(original_epochs.data)
    end
    
    @testset "Batch with participant filtering" begin
        output_dir = joinpath(test_dir, "rejected_filtered")
        
        # Process only participants 1 and 2
        result = eegfun.reject_epochs_automatic("epochs", 2.0,
                                                input_dir = test_dir,
                                                participants = [1, 2],
                                                output_dir = output_dir)
        
        @test isdir(output_dir)
        
        # Should only have 2 output files (plus reports and log)
        jld2_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir))
        @test length(jld2_files) == 2
        @test "1_epochs.jld2" in jld2_files
        @test "2_epochs.jld2" in jld2_files
        @test !("3_epochs.jld2" in jld2_files)
    end
    
    @testset "Batch error handling: no matching files" begin
        output_dir = joinpath(test_dir, "rejected_nomatch")
        
        result = eegfun.reject_epochs_automatic("nonexistent", 2.0,
                                                input_dir = test_dir,
                                                output_dir = output_dir)
        
        # Should return nothing when no files found
        @test result === nothing
    end
    
    @testset "Batch logging" begin
        output_dir = joinpath(test_dir, "rejected_logging")
        
        result = eegfun.reject_epochs_automatic("epochs", 2.0,
                                                input_dir = test_dir,
                                                output_dir = output_dir)
        
        # Check that log file was created
        log_file = joinpath(output_dir, "reject_epochs_automatic.log")
        @test isfile(log_file)
        
        # Verify log content
        log_content = read(log_file, String)
        @test occursin("Batch automatic epoch rejection started", log_content)
        @test occursin("Z-criterion", log_content)
        @test occursin("Found", log_content)
    end
    
    # Cleanup
    rm(test_dir, recursive = true)
end

