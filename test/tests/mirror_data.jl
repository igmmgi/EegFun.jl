"""
Tests for data mirroring functions.

Tests mirror_data() and unmirror_data() for both EpochData and ErpData.
"""

using DataFrames

@testset "Data Mirroring" begin
    
    @testset "Mirror EpochData - :pre" begin
        # Create simple test epoch
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1)
        )
        
        epoch2 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples) .* 2,
            Pz = collect(n_samples:-1.0:1) .* 2
        )
        
        epochs = eegfun.EpochData([epoch1, epoch2], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        original_length = nrow(epochs.data[1])
        
        # Mirror :pre
        eegfun.mirror_data!(epochs, :pre)
        
        # Check length increased
        @test nrow(epochs.data[1]) > original_length
        
        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])  # Uniform spacing
        
        # Unmirror
        eegfun.unmirror_data!(epochs, :pre)
        
        # Check restored to original
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ collect(time)
        @test epochs.data[1].Cz ≈ collect(1.0:n_samples)
    end
    
    
    @testset "Mirror EpochData - :post" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1)
        )
        
        epochs = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        original_length = nrow(epochs.data[1])
        original_data = copy(epochs.data[1])
        
        # Mirror :post
        eegfun.mirror_data!(epochs, :post)
        
        # Check length increased
        @test nrow(epochs.data[1]) > original_length
        
        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])
        
        # Unmirror
        eegfun.unmirror_data!(epochs, :post)
        
        # Check restored
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ original_data.time
        @test epochs.data[1].Cz ≈ original_data.Cz
    end
    
    
    @testset "Mirror EpochData - :both" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1),
            trial = fill(1, n_samples),
            condition = fill(1, n_samples)
        )
        
        epochs = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        original_length = nrow(epochs.data[1])
        original_time = copy(epochs.data[1].time)
        original_cz = copy(epochs.data[1].Cz)
        
        # Mirror :both
        eegfun.mirror_data!(epochs, :both)
        
        # Check length (should be roughly 3× original)
        mirrored_length = nrow(epochs.data[1])
        @test mirrored_length > 2 * original_length
        
        # Check time vector continuity
        time_diffs = diff(epochs.data[1].time)
        @test all(time_diffs .≈ time_diffs[1])
        
        # Unmirror
        eegfun.unmirror_data!(epochs, :both)
        
        # Check fully restored
        @test nrow(epochs.data[1]) == original_length
        @test epochs.data[1].time ≈ original_time
        @test epochs.data[1].Cz ≈ original_cz
    end
    
    
    @testset "Mirror EpochData - Non-mutating" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples)
        )
        
        epochs_original = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        original_length = nrow(epochs_original.data[1])
        
        # Non-mutating mirror
        epochs_mirrored = eegfun.mirror_data(epochs_original, :both)
        
        # Check original unchanged
        @test nrow(epochs_original.data[1]) == original_length
        
        # Check mirrored is different
        @test nrow(epochs_mirrored.data[1]) > original_length
        
        # Non-mutating unmirror
        epochs_unmirrored = eegfun.unmirror_data(epochs_mirrored, :both)
        
        # Check unmirrored matches original
        @test nrow(epochs_unmirrored.data[1]) == original_length
        @test epochs_unmirrored.data[1].time ≈ epochs_original.data[1].time
        @test epochs_unmirrored.data[1].Cz ≈ epochs_original.data[1].Cz
    end
    
    
    @testset "Mirror ErpData - :both" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        erp_df = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1)
        )
        
        erp = eegfun.ErpData(erp_df, eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo(), 10)
        
        original_length = nrow(erp.data)
        original_time = copy(erp.data.time)
        original_cz = copy(erp.data.Cz)
        
        # Mirror
        eegfun.mirror_data!(erp, :both)
        
        # Check length increased
        @test nrow(erp.data) > 2 * original_length
        
        # Check time continuity
        time_diffs = diff(erp.data.time)
        @test all(time_diffs .≈ time_diffs[1])
        
        # Unmirror
        eegfun.unmirror_data!(erp, :both)
        
        # Check restored
        @test nrow(erp.data) == original_length
        @test erp.data.time ≈ original_time
        @test erp.data.Cz ≈ original_cz
    end
    
    
    @testset "Mirror ErpData - Non-mutating" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        erp_df = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples)
        )
        
        erp_original = eegfun.ErpData(erp_df, eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo(), 10)
        original_length = nrow(erp_original.data)
        
        # Non-mutating mirror
        erp_mirrored = eegfun.mirror_data(erp_original, :both)
        
        # Check original unchanged
        @test nrow(erp_original.data) == original_length
        
        # Check mirrored is different
        @test nrow(erp_mirrored.data) > original_length
        
        # Non-mutating unmirror
        erp_unmirrored = eegfun.unmirror_data(erp_mirrored, :both)
        
        # Check matches original
        @test nrow(erp_unmirrored.data) == original_length
        @test erp_unmirrored.data.time ≈ erp_original.data.time
        @test erp_unmirrored.data.Cz ≈ erp_original.data.Cz
    end
    
    
    @testset "Invalid side parameter" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples)
        )
        
        epochs = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        # Test invalid side
        @test_throws Exception mirror_data!(epochs, :invalid)
        @test_throws Exception unmirror_data!(epochs, :middle)
    end
    
    
    @testset "Multiple epochs" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        # Create 3 different epochs
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples)
        )
        
        epoch2 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples) .* 2
        )
        
        epoch3 = DataFrame(
            time = collect(time),
            Cz = collect(n_samples:-1.0:1)
        )
        
        epochs = eegfun.EpochData([epoch1, epoch2, epoch3], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        original_cz1 = copy(epochs.data[1].Cz)
        original_cz2 = copy(epochs.data[2].Cz)
        original_cz3 = copy(epochs.data[3].Cz)
        
        # Mirror and unmirror
        eegfun.mirror_data!(epochs, :both)
        eegfun.unmirror_data!(epochs, :both)
        
        # Check all epochs restored correctly
        @test epochs.data[1].Cz ≈ original_cz1
        @test epochs.data[2].Cz ≈ original_cz2
        @test epochs.data[3].Cz ≈ original_cz3
    end
    
    
    @testset "Mirroring preserves metadata columns" begin
        time = -0.2:0.1:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = collect(1.0:n_samples),
            Pz = collect(n_samples:-1.0:1),
            trial = fill(1, n_samples),
            condition = fill(2, n_samples),
            response = fill("left", n_samples)
        )
        
        epochs = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        # Mirror
        eegfun.mirror_data!(epochs, :both)
        
        # Check metadata preserved
        @test all(epochs.data[1].trial .== 1)
        @test all(epochs.data[1].condition .== 2)
        @test all(epochs.data[1].response .== "left")
        
        # Unmirror
        eegfun.unmirror_data!(epochs, :both)
        
        # Check metadata still there
        @test all(epochs.data[1].trial .== 1)
        @test all(epochs.data[1].condition .== 2)
        @test all(epochs.data[1].response .== "left")
    end
    
    
    @testset "Integration with filter" begin
        # This test would require the filter function
        # Commented out for now, but shows intended usage
        
        # time = -0.5:0.001:0.5
        # n_samples = length(time)
        # 
        # # Create epoch with known frequency
        # epoch1 = DataFrame(
        #     time = collect(time),
        #     Cz = sin.(2π * 10 .* time)  # 10 Hz sine wave
        # )
        # 
        # epochs = EpochData([epoch1], Layout(), 1000.0, AnalysisInfo())
        # 
        # # Mirror, filter, unmirror
        # mirror_data!(epochs, :both)
        # filter!(epochs, 5.0, 15.0)  # Bandpass around 10 Hz
        # unmirror_data!(epochs, :both)
        # 
        # # Check epoch restored to original length
        # @test nrow(epochs.data[1]) == n_samples
    end
    
    
    @testset "Roundtrip consistency" begin
        # Test that mirror → unmirror → mirror → unmirror works
        time = -0.2:0.05:0.2
        n_samples = length(time)
        
        epoch1 = DataFrame(
            time = collect(time),
            Cz = randn(n_samples),
            Pz = randn(n_samples)
        )
        
        epochs = eegfun.EpochData([epoch1], eegfun.Layout(DataFrame(), nothing, nothing), 100.0, eegfun.AnalysisInfo())
        
        original_cz = copy(epochs.data[1].Cz)
        original_pz = copy(epochs.data[1].Pz)
        original_time = copy(epochs.data[1].time)
        
        # First roundtrip
        eegfun.mirror_data!(epochs, :both)
        eegfun.unmirror_data!(epochs, :both)
        
        @test epochs.data[1].Cz ≈ original_cz
        @test epochs.data[1].Pz ≈ original_pz
        @test epochs.data[1].time ≈ original_time
        
        # Second roundtrip
        eegfun.mirror_data!(epochs, :both)
        eegfun.unmirror_data!(epochs, :both)
        
        @test epochs.data[1].Cz ≈ original_cz
        @test epochs.data[1].Pz ≈ original_pz
        @test epochs.data[1].time ≈ original_time
    end
end

