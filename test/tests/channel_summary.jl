using Test
using DataFrames
using Statistics
using eegfun

@testset "eegfun" begin
    @testset "channel_summary" begin

    # Helper function to create test data
    function create_continuous_test_data()
        fs = 100.0
        t = 0.0:1/fs:1.0  # 1 second of data
        n_samples = length(t)
        
        # Create test signals with known properties
        A = 2.0 .* sin.(2π .* 5 .* t) .+ randn(n_samples) * 0.1  # 5 Hz sine + noise
        B = 3.0 .* cos.(2π .* 10 .* t) .+ randn(n_samples) * 0.2  # 10 Hz cosine + noise  
        C = ones(n_samples) * 1.5  # Constant signal
        triggers = vcat(fill(0, 30), fill(1, 20), fill(0, 51))  # Simple trigger pattern
        
        df = DataFrame(time = collect(t), triggers = triggers, A = A, B = B, C = C)
        layout = eegfun.Layout(DataFrame(label = [:A, :B, :C], inc = [0.0, 0.0, 0.0], azi = [0.0, 90.0, 180.0]), nothing, nothing)
        dat = eegfun.ContinuousData(copy(df, copycols=true), layout, fs, eegfun.AnalysisInfo())
        return dat
    end
    
    function create_epoch_test_data()
        # Create continuous data first
        continuous_dat = create_continuous_test_data()
        
        # Create epoch condition
        ec = eegfun.EpochCondition(name = "test", trigger_sequences = [[1]], reference_index = 1)
        
        # Extract epochs
        epoch_dat = eegfun.extract_epochs(continuous_dat, 1, ec, -0.1, 0.2)
        return epoch_dat
    end

    @testset "_channel_summary_impl core function" begin
        dat = create_continuous_test_data()
        
        # Test basic functionality
        sample_indices = 1:50  # First half of data
        channel_names = [:A, :B]
        result = eegfun._channel_summary_impl(dat.data, collect(sample_indices), channel_names)
        
        @test result isa DataFrame
        @test nrow(result) == 2  # Two channels
        @test :channel in propertynames(result)
        @test :min in propertynames(result)
        @test :max in propertynames(result)
        @test :std in propertynames(result)
        @test :range in propertynames(result)
        @test :var in propertynames(result)
        @test :zvar in propertynames(result)
        
        # Test column order
        @test propertynames(result)[1] == :channel
        
        # Test range calculation
        @test all(result.range .≈ result.max .- result.min)
        
        # Test variance calculation
        selected_data = dat.data[collect(sample_indices), channel_names]
        expected_var = var.(eachcol(selected_data))
        @test result.var ≈ expected_var
        
        # Test that channels are correctly identified
        @test Set(result.channel) == Set(channel_names)
    end
    
    @testset "_channel_summary_impl input validation" begin
        dat = create_continuous_test_data()
        
        # Test empty samples
        @test_throws Exception eegfun._channel_summary_impl(dat.data, Int[], [:A, :B])
        
        # Test empty channels
        @test_throws Exception eegfun._channel_summary_impl(dat.data, [1, 2, 3], Symbol[])
        
        # Test invalid channel names
        @test_throws Exception eegfun._channel_summary_impl(dat.data, [1, 2, 3], [:NonExistent])
        
        # Test invalid sample indices
        @test_throws Exception eegfun._channel_summary_impl(dat.data, [1000], [:A])  # Out of bounds
        @test_throws Exception eegfun._channel_summary_impl(dat.data, [0], [:A])     # Below bounds
    end
    
    @testset "_channel_summary_impl edge cases" begin
        dat = create_continuous_test_data()
        
        # Test single sample
        result = eegfun._channel_summary_impl(dat.data, [1], [:A])
        @test nrow(result) == 1
        @test result.range[1] == 0.0  # min == max for single sample
        @test isnan(result.std[1]) || result.std[1] == 0.0  # std is NaN for single sample in Julia
        
        # Test single channel
        result = eegfun._channel_summary_impl(dat.data, collect(1:10), [:A])
        @test nrow(result) == 1
        @test result.channel[1] == :A
        
        # Test with constant signal (zero variance case)
        result = eegfun._channel_summary_impl(dat.data, collect(1:50), [:C])  # C is constant
        @test result.var[1] ≈ 0.0
        @test result.zvar[1] == 0.0  # Should handle zero variance case
    end

    @testset "channel_summary SingleDataFrameEeg" begin
        dat = create_continuous_test_data()
        
        # Test basic functionality with defaults
        result = eegfun.channel_summary(dat)
        @test result isa DataFrame
        @test nrow(result) == 3  # Three layout channels (A, B, C)
        @test :channel in propertynames(result)
        
        # Test that all layout channels are included by default
        @test Set(result.channel) == Set([:A, :B, :C])
        
        # Test specific channel selection
        result_subset = eegfun.channel_summary(dat, channel_selection = eegfun.channels([:A, :B]))
        @test nrow(result_subset) == 2
        @test Set(result_subset.channel) == Set([:A, :B])
        
        # Test sample selection (first half)
        n_samples = nrow(dat.data)
        result_samples = eegfun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
        @test result_samples isa DataFrame
        @test nrow(result_samples) == 3
        
        # Test include_meta and include_extra flags
        result_meta = eegfun.channel_summary(dat, include_meta = true)
        @test :time in result_meta.channel || :triggers in result_meta.channel  # Should include meta columns
        
        # Test statistics are reasonable
        result = eegfun.channel_summary(dat)
        @test all(result.range .>= 0)  # Range should be non-negative
        @test all(result.var .>= 0)    # Variance should be non-negative
        @test all(result.std .>= 0)    # Standard deviation should be non-negative
    end
    
    @testset "channel_summary SingleDataFrameEeg input validation" begin
        # Test empty data
        empty_df = DataFrame()
        layout = eegfun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing)
        empty_dat = eegfun.ContinuousData(empty_df, layout, 100.0, eegfun.AnalysisInfo())
        @test_throws Exception eegfun.channel_summary(empty_dat)
    end

    @testset "channel_summary MultiDataFrameEeg" begin
        epoch_dat = create_epoch_test_data()
        
        # Test basic functionality
        result = eegfun.channel_summary(epoch_dat)
        @test result isa DataFrame
        @test :epoch in propertynames(result)
        @test :channel in propertynames(result)
        
        # Test that we have results for each epoch and channel
        n_epochs = eegfun.n_epochs(epoch_dat)
        n_channels = 3  # A, B, C
        @test nrow(result) == n_epochs * n_channels
        
        # Test epoch numbers are preserved
        original_epoch_numbers = [df.epoch[1] for df in epoch_dat.data]
        result_epoch_numbers = unique(result.epoch)
        @test Set(result_epoch_numbers) == Set(original_epoch_numbers)
        
        # Test channel selection works across epochs
        result_subset = eegfun.channel_summary(epoch_dat, channel_selection = eegfun.channels([:A]))
        @test nrow(result_subset) == n_epochs  # One channel per epoch
        @test all(result_subset.channel .== :A)
        
        # Test that statistics vary across epochs (should be different)
        epochs_A = result[result.channel .== :A, :]
        if nrow(epochs_A) > 1
            # Should have some variation across epochs (unless data is identical)
            @test length(unique(epochs_A.var)) > 1 || all(epochs_A.var .== 0)
        end
    end
    
    @testset "channel_summary MultiDataFrameEeg input validation" begin
        # Test empty epoch data
        empty_layout = eegfun.Layout(DataFrame(label = Symbol[], inc = Float64[], azi = Float64[]), nothing, nothing)
        empty_epochs = eegfun.EpochData(DataFrame[], empty_layout, 100.0, eegfun.AnalysisInfo())
        @test_throws Exception eegfun.channel_summary(empty_epochs)
    end
    
    @testset "channel_summary statistical properties" begin
        dat = create_continuous_test_data()
        result = eegfun.channel_summary(dat)
        
        # Test that constant channel has zero variance
        constant_row = result[result.channel .== :C, :]
        @test nrow(constant_row) == 1
        @test constant_row.var[1] ≈ 0.0 atol=1e-10
        @test constant_row.std[1] ≈ 0.0 atol=1e-10
        @test constant_row.range[1] ≈ 0.0 atol=1e-10
        
        # Test that varying channels have non-zero variance
        varying_rows = result[result.channel .!= :C, :]
        @test all(varying_rows.var .> 0)
        @test all(varying_rows.std .> 0)
        @test all(varying_rows.range .> 0)
        
        # Test z-scored variance properties
        if length(unique(result.var)) > 1  # If there's variance in variances
            @test abs(mean(result.zvar)) < 1e-10  # Should be approximately zero mean
            @test abs(std(result.zvar) - 1.0) < 1e-10  # Should have unit std
        end
    end
    
    @testset "channel_summary consistency" begin
        dat = create_continuous_test_data()
        
        # Test that results are consistent between calls
        result1 = eegfun.channel_summary(dat)
        result2 = eegfun.channel_summary(dat)
        @test result1 == result2
        
        # Test that subsetting gives expected results
        full_result = eegfun.channel_summary(dat)
        subset_result = eegfun.channel_summary(dat, channel_selection = eegfun.channels([:A, :B]))
        
        # The subset should match the corresponding rows from the full result
        full_subset = full_result[in([:A, :B]).(full_result.channel), :]
        @test subset_result.channel == full_subset.channel
        @test subset_result.min ≈ full_subset.min
        @test subset_result.max ≈ full_subset.max
        @test subset_result.std ≈ full_subset.std
    end
    
    @testset "channel_summary integration with selection functions" begin
        dat = create_continuous_test_data()
        
        # Test with different selection functions
        # Note: These functions are defined in the eegfun package
        result_all = eegfun.channel_summary(dat, channel_selection = eegfun.channels())
        @test nrow(result_all) == 3
        
        # Test channel selection by name
        result_specific = eegfun.channel_summary(dat, channel_selection = eegfun.channels([:A]))
        @test nrow(result_specific) == 1
        @test result_specific.channel[1] == :A
        
        # Test sample selection by range
        result_half = eegfun.channel_summary(dat, sample_selection = x -> 1:nrow(x) .<= div(nrow(x), 2))
        @test result_half isa DataFrame
        @test nrow(result_half) == 3  # Same number of channels
    end

    end # channel_summary testset
end # eegfun testset
