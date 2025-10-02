"""
Test suite for src/analysis/realign.jl
"""

using Test
using JLD2
using DataFrames

# Helper function to create test epoched data with realignment times
function create_test_epoch_data_with_rt(
    n_epochs::Int = 10,
    n_timepoints::Int = 200,
    n_channels::Int = 3;
    epoch_start::Float64 = -0.5,
    epoch_end::Float64 = 1.5,
    rt_range::Tuple{Float64,Float64} = (0.3, 0.7)
)
    # Create time vector (stimulus-locked)
    time = collect(range(epoch_start, epoch_end, length = n_timepoints))
    sample_rate = Int(round(n_timepoints / (epoch_end - epoch_start)))
    
    # Generate random RTs for each epoch
    rts = sort(rand(rt_range[1]:0.01:rt_range[2], n_epochs))
    
    # Create epochs
    epochs = DataFrame[]
    for i = 1:n_epochs
        epoch_df = DataFrame()
        epoch_df.time = copy(time)
        epoch_df.trial = fill(i, n_timepoints)
        epoch_df.condition = fill(1, n_timepoints)
        epoch_df.rt = fill(rts[i], n_timepoints)  # RT constant within epoch
        
        # Add EEG channels with some pattern
        for ch_idx = 1:n_channels
            ch_name = Symbol("Ch$ch_idx")
            # Simple signal that peaks around RT
            signal = exp.(-((time .- rts[i]) .^ 2) ./ 0.05) .* 10.0 .+ randn(n_timepoints) .* 0.5
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
    
    return eegfun.EpochData(epochs, layout, sample_rate, eegfun.AnalysisInfo())
end

@testset "Realignment" begin
    @testset "Basic realignment (non-mutating)" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Store original time range
        original_time_min = minimum(epoch_data.data[1].time)
        original_time_max = maximum(epoch_data.data[1].time)
        
        # Realign to RT
        realigned = eegfun.realign(epoch_data, :rt)
        
        # Check that original data is unchanged
        @test minimum(epoch_data.data[1].time) ≈ original_time_min
        @test maximum(epoch_data.data[1].time) ≈ original_time_max
        
        # Check that realigned data has time 0 near RT
        # After realignment, the RT column should now be approximately 0
        # (actually, it should be exactly the negative of the original RT, but after cropping it should be near 0)
        for epoch in realigned.data
            # The rt column value should be 0 after realignment
            # (since we subtracted it from the time column)
            @test all(epoch.rt .≈ 0.0)
        end
        
        # Check that all epochs have the same time range after realignment
        for i = 2:length(realigned.data)
            @test realigned.data[i].time ≈ realigned.data[1].time
        end
        
        # Check structure
        @test realigned isa eegfun.EpochData
        @test length(realigned.data) == 10
    end
    
    @testset "In-place realignment" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Get original RT values for verification
        original_rts = [epoch.rt[1] for epoch in epoch_data.data]
        
        # Realign in-place
        eegfun.realign!(epoch_data, :rt)
        
        # Check that data was modified
        # RT column should now be 0 (or very close to 0)
        for epoch in epoch_data.data
            @test all(abs.(epoch.rt) .< 1e-10)
        end
        
        # Check that all epochs have the same time range
        for i = 2:length(epoch_data.data)
            @test epoch_data.data[i].time ≈ epoch_data.data[1].time
        end
    end
    
    @testset "Time window cropping" begin
        # Create epochs with varying RTs
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3, 
                                                     epoch_start = -0.5, 
                                                     epoch_end = 1.5,
                                                     rt_range = (0.3, 0.7))
        
        # Get time ranges before realignment
        original_lengths = [nrow(epoch) for epoch in epoch_data.data]
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # Check that epochs are cropped
        realigned_lengths = [nrow(epoch) for epoch in realigned.data]
        
        # All realigned epochs should have the same length
        @test all(realigned_lengths .== realigned_lengths[1])
        
        # Realigned length should be less than or equal to original
        # (equality happens if all RTs are identical)
        @test realigned_lengths[1] <= original_lengths[1]
        
        # Check that time windows are identical across epochs
        for i = 1:length(realigned.data)
            @test minimum(realigned.data[i].time) ≈ minimum(realigned.data[1].time)
            @test maximum(realigned.data[i].time) ≈ maximum(realigned.data[1].time)
        end
    end
    
    @testset "Common time window calculation" begin
        # Create epochs with RTs that span a wide range
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3,
                                                     epoch_start = -0.5,
                                                     epoch_end = 1.5,
                                                     rt_range = (0.2, 1.0))
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # The common window should accommodate all trials
        # Early RTs (e.g., 0.2s) limit how much pre-response data we can have
        # Late RTs (e.g., 1.0s) limit how much post-response data we can have
        
        # With epoch from -0.5 to 1.5 and RT at 0.2:
        #   After realignment: -0.7 to 1.3
        # With epoch from -0.5 to 1.5 and RT at 1.0:
        #   After realignment: -1.5 to 0.5
        # Common window: max(-0.7, -1.5) to min(1.3, 0.5) = -0.7 to 0.5
        
        # Just verify that we got some reasonable window
        time_range = maximum(realigned.data[1].time) - minimum(realigned.data[1].time)
        @test time_range > 0
        @test time_range < 2.0  # Less than original epoch length
    end
    
    @testset "Channel data preservation" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Get channel values at a specific time point before realignment
        # Find time closest to 0.5s in first epoch
        time_idx = argmin(abs.(epoch_data.data[1].time .- 0.5))
        original_value = epoch_data.data[1].Ch1[time_idx]
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # Channel data should be preserved (just shifted in time)
        # We can't easily verify the exact values without knowing the RT,
        # but we can verify that the data range is similar
        @test minimum(realigned.data[1].Ch1) ≈ minimum(epoch_data.data[1].Ch1) atol = 1.0
        @test maximum(realigned.data[1].Ch1) ≈ maximum(epoch_data.data[1].Ch1) atol = 1.0
    end
    
    @testset "Metadata preservation" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # Check that metadata columns are preserved
        for i = 1:length(realigned.data)
            @test hasproperty(realigned.data[i], :trial)
            @test hasproperty(realigned.data[i], :condition)
            @test hasproperty(realigned.data[i], :rt)
            
            # Trial numbers should be preserved
            @test all(realigned.data[i].trial .== epoch_data.data[i].trial[1])
        end
        
        # Layout and other metadata should be preserved
        @test realigned.layout.data == epoch_data.layout.data
        @test realigned.sample_rate == epoch_data.sample_rate
    end
    
    @testset "Error handling: missing column" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Try to realign to non-existent column
        @test_throws Exception eegfun.realign(epoch_data, :nonexistent_column)
    end
    
    @testset "Error handling: varying values within epoch" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Modify RT column to have varying values within first epoch
        epoch_data.data[1].rt .= collect(1:nrow(epoch_data.data[1]))
        
        # Should error because RT should be constant within epoch
        @test_throws Exception eegfun.realign(epoch_data, :rt)
    end
    
    @testset "Error handling: non-finite values" begin
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Set RT to NaN in first epoch
        epoch_data.data[1].rt .= NaN
        
        # Should error
        @test_throws Exception eegfun.realign(epoch_data, :rt)
    end
    
    @testset "Error handling: insufficient epoch length" begin
        # Create epochs that are too short for the RT values
        epoch_data = create_test_epoch_data_with_rt(10, 50, 3,
                                                     epoch_start = 0.0,
                                                     epoch_end = 0.3,
                                                     rt_range = (0.25, 0.28))
        
        # This might work or fail depending on exact RTs
        # If it fails, it should fail gracefully with a clear message
        try
            realigned = eegfun.realign(epoch_data, :rt)
            # If it succeeds, check that result is valid
            @test all(length(realigned.data[i].time) > 0 for i = 1:length(realigned.data))
        catch e
            # Should be a clear error message about insufficient epoch length
            @test occursin("common time window", sprint(showerror, e)) ||
                  occursin("No samples found", sprint(showerror, e))
        end
    end
    
    @testset "Multiple channels preserved" begin
        # Create data with more channels
        epoch_data = create_test_epoch_data_with_rt(5, 100, 10)
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # All channels should be present
        for i = 1:10
            ch_name = Symbol("Ch$i")
            @test hasproperty(realigned.data[1], ch_name)
        end
        
        # Number of channels should match
        n_channels_original = length(eegfun.channel_labels(epoch_data))
        n_channels_realigned = length(eegfun.channel_labels(realigned))
        @test n_channels_original == n_channels_realigned
    end
    
    @testset "Edge case: identical RTs" begin
        # Create epochs where all RTs are identical
        epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
        
        # Set all RTs to the same value
        for epoch in epoch_data.data
            epoch.rt .= 0.5
        end
        
        # Realign
        realigned = eegfun.realign(epoch_data, :rt)
        
        # With identical RTs, all epochs should have the same length as before
        # (no cropping needed)
        @test nrow(realigned.data[1]) == nrow(epoch_data.data[1])
        
        # All epochs should have identical time vectors
        for i = 2:length(realigned.data)
            @test realigned.data[i].time == realigned.data[1].time
        end
    end
    
    @testset "Realistic use case: response-locked LRP" begin
        # Simulate a realistic scenario for response-locked LRP
        # Stimulus-locked epochs from -0.5 to 2.0s
        # RTs varying from 0.4 to 1.2s
        epoch_data = create_test_epoch_data_with_rt(20, 300, 4,
                                                     epoch_start = -0.5,
                                                     epoch_end = 2.0,
                                                     rt_range = (0.4, 1.2))
        
        # Store original RTs for verification
        original_rts = [epoch.rt[1] for epoch in epoch_data.data]
        
        # Realign to response
        realigned = eegfun.realign(epoch_data, :rt)
        
        # Verify that:
        # 1. All epochs have the same time vector
        for i = 2:length(realigned.data)
            @test realigned.data[i].time ≈ realigned.data[1].time
        end
        
        # 2. RT column is now 0 (or very close)
        for epoch in realigned.data
            @test all(abs.(epoch.rt) .< 1e-10)
        end
        
        # 3. We have some pre-response and post-response data
        @test minimum(realigned.data[1].time) < 0
        @test maximum(realigned.data[1].time) > 0
        
        # 4. The resulting epoch is shorter than the original
        # (because we need to accommodate varying RTs)
        @test nrow(realigned.data[1]) < nrow(epoch_data.data[1])
    end
end

@testset "Batch realignment" begin
    # Create temporary test directory
    test_dir = mktempdir()
    
    @testset "Basic batch processing" begin
        # Create test epoch files for multiple participants
        for participant = 1:3
            epoch_data = create_test_epoch_data_with_rt(10, 200, 3)
            file_path = joinpath(test_dir, "$(participant)_epochs.jld2")
            save(file_path, "epochs", epoch_data)
        end
        
        output_dir = joinpath(test_dir, "realigned_test")
        
        # Test basic realignment
        result = eegfun.realign("epochs", :rt, input_dir = test_dir, output_dir = output_dir)
        
        # Verify output directory was created
        @test isdir(output_dir)
        
        # Verify output files
        output_files = readdir(output_dir)
        @test "1_epochs.jld2" in output_files
        @test "2_epochs.jld2" in output_files
        @test "3_epochs.jld2" in output_files
        
        # Load and verify one file
        realigned = load(joinpath(output_dir, "1_epochs.jld2"), "epochs")
        @test realigned isa eegfun.EpochData
        
        # Check that RT is now 0
        for epoch in realigned.data
            @test all(abs.(epoch.rt) .< 1e-10)
        end
    end
    
    @testset "Batch with participant filtering" begin
        output_dir = joinpath(test_dir, "realigned_filtered")
        
        # Process only participants 1 and 2
        result = eegfun.realign("epochs", :rt, 
                               input_dir = test_dir, 
                               participants = [1, 2],
                               output_dir = output_dir)
        
        @test isdir(output_dir)
        
        # Should only have 2 output files
        output_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir))
        @test length(output_files) == 2
        @test "1_epochs.jld2" in output_files
        @test "2_epochs.jld2" in output_files
        @test !("3_epochs.jld2" in output_files)
    end
    
    @testset "Batch error handling: no matching files" begin
        output_dir = joinpath(test_dir, "realigned_nomatch")
        
        result = eegfun.realign("nonexistent", :rt,
                               input_dir = test_dir,
                               output_dir = output_dir)
        
        # Should return nothing when no files found
        @test result === nothing
    end
    
    @testset "Batch logging" begin
        output_dir = joinpath(test_dir, "realigned_logging")
        
        result = eegfun.realign("epochs", :rt,
                               input_dir = test_dir,
                               output_dir = output_dir)
        
        # Check that log file was created
        log_file = joinpath(output_dir, "realign.log")
        @test isfile(log_file)
        
        # Verify log content
        log_content = read(log_file, String)
        @test occursin("Batch realignment started", log_content)
        @test occursin("Found", log_content)
        @test occursin("Realigning to column: :rt", log_content)
    end
    
    # Cleanup
    rm(test_dir, recursive = true)
end

