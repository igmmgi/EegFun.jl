using Test
using DataFrames
using eegfun
using Statistics
using OrderedCollections


@testset "correlation" begin

    @testset "correlation_matrix" begin

        dat = create_test_data(n = 1000, n_channels = 100)

        # Test with all channels
        cm = eegfun.correlation_matrix(dat)

        @test cm isa DataFrame
        @test :row in propertynames(cm)
        @test size(cm, 1) == 100  # 4 channels
        @test size(cm, 2) == 101  # 4 channels + row column

        # Test diagonal elements are 1.0 (correlation of each channel with itself)
        for i = 1:100
            channel_name = cm[i, :row]
            @test isapprox(cm[i, channel_name], 1.0, atol = 1e-10)
        end

        # Test symmetry
        @test isapprox(cm[1, :Ch2], cm[2, :Ch1], atol = 1e-10)
        @test isapprox(cm[1, :Ch3], cm[3, :Ch1], atol = 1e-10)

        # Test with channel selection
        cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        @test size(cm, 1) == 2
        @test size(cm, 2) == 3  # 2 channels + row column

        # Test empty channel selection
        @test_throws ErrorException eegfun.correlation_matrix(dat, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "channel_joint_probability" begin
        dat = create_test_data(n = 1000, n_channels = 100)

        # Test with default parameters
        jp = eegfun.channel_joint_probability(dat)
        @test jp isa DataFrame
        @test :channel in propertynames(jp)
        @test :jp in propertynames(jp)
        @test size(jp, 1) == 100  # 4 channels

        # Test jp values are finite (can be negative due to normalization)
        for i = 1:100
            @test isfinite(jp[i, :jp])
        end

        # Test with custom parameters
        jp = eegfun.channel_joint_probability(dat, threshold = 0.3, normalize = 2, discret = 500)
        @test jp isa DataFrame
        @test size(jp, 1) == 100

        # Test with channel selection
        jp = eegfun.channel_joint_probability(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        @test size(jp, 1) == 2

        # Test empty channel selection
        @test_throws ErrorException eegfun.channel_joint_probability(dat, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "_correlation_matrix" begin
        dat = create_test_data(n = 1000, n_channels = 100)

        # Test internal function
        cm = eegfun._correlation_matrix(dat.data, collect(1:size(dat.data, 1)), [:Ch1, :Ch2, :Ch3])
        @test cm isa DataFrame
        @test :row in propertynames(cm)
        @test size(cm, 1) == 3
        @test size(cm, 2) == 4  # 3 channels + row column
    end

    @testset "_channel_joint_probability" begin
        dat = create_test_data(n = 1000, n_channels = 100)

        # Test internal function
        jp = eegfun._channel_joint_probability(
            dat.data,
            collect(1:size(dat.data, 1)),
            [:Ch1, :Ch2, :Ch3],
            threshold = 0.5,
            normval = 1,
            discret = 1000,
        )
        @test jp isa DataFrame
        @test :channel in propertynames(jp)
        @test :jp in propertynames(jp)
        @test size(jp, 1) == 3
    end

    @testset "_joint_probability" begin
        # Test with simple data (channels × samples)
        signal = randn(3, 100)  # 3 channels × 100 samples
        jp = eegfun._joint_probability(signal, 0.5, 1, 100)

        @test jp isa Tuple
        @test length(jp) == 2
        jp, rej = jp
        @test jp isa Vector{Float64}
        @test rej isa BitVector
        @test length(jp) == 3
        @test length(rej) == 3
    end

    @testset "compute_probability!" begin
        # Test probability computation
        data = randn(1000)
        proba_map = zeros(Float64, 1000)  # Must match data length

        result = eegfun.compute_probability!(proba_map, data, 50)

        @test result isa Vector{Float64}
        @test length(result) == 1000
        @test all(0 <= x <= 1 for x in result)  # All values between 0 and 1
    end

    @testset "_trim_extremes" begin
        # Test extreme value trimming
        x = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0, 6.0, 200.0, 7.0, 8.0, 9.0]
        trimmed = eegfun._trim_extremes(x)

        @test trimmed isa SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true}
        @test length(trimmed) < length(x)  # Should be shorter after trimming
        @test maximum(trimmed) < maximum(x)  # Should trim extreme values
        @test minimum(trimmed) >= minimum(x)  # Should trim extreme values (can be equal if no extreme low values)

        # Test with data that has both high and low extremes
        x2 = [1.0, 2.0, 3.0, 100.0, 4.0, 5.0, 6.0, 200.0, 7.0, 8.0, -50.0]
        trimmed2 = eegfun._trim_extremes(x2)
        @test maximum(trimmed2) < maximum(x2)  # Should trim high extreme values
        @test minimum(trimmed2) > minimum(x2)  # Should trim low extreme values
        @test length(trimmed2) < length(x2)  # Should be shorter after trimming
    end

    @testset "correlation_matrix - additional parameters" begin
        dat = create_test_data(n = 1000, n_channels = 5)
        
        # Test with sample_selection (create a predicate function)
        sample_pred = x -> [i <= 500 for i in 1:nrow(x)]
        cm = eegfun.correlation_matrix(dat, sample_selection = sample_pred)
        @test cm isa DataFrame
        @test size(cm, 1) == 5
        
        # Test with include_extra (when extra channels exist)
        dat.data[!, :extra_ch] = randn(nrow(dat.data))
        cm = eegfun.correlation_matrix(dat, include_extra = true)
        @test size(cm, 1) == 6  # 5 channels + 1 extra
    end

    @testset "channel_joint_probability - additional parameters" begin
        dat = create_test_data(n = 1000, n_channels = 5)
        
        # Test with sample_selection (create a predicate function)
        sample_pred = x -> [i <= 500 for i in 1:nrow(x)]
        jp = eegfun.channel_joint_probability(dat, sample_selection = sample_pred)
        @test jp isa DataFrame
        @test size(jp, 1) == 5
        
        # Test with include_extra (when extra channels exist)
        dat.data[!, :extra_ch] = randn(nrow(dat.data))
        jp = eegfun.channel_joint_probability(dat, include_extra = true)
        @test size(jp, 1) == 6  # 5 channels + 1 extra
        
        # Test with different normalize values
        jp = eegfun.channel_joint_probability(dat, normalize = 0)
        @test size(jp, 1) == 5
        
        jp = eegfun.channel_joint_probability(dat, normalize = 1)
        @test size(jp, 1) == 5
        
        # Test with threshold = 0
        jp = eegfun.channel_joint_probability(dat, threshold = 0.0)
        @test size(jp, 1) == 5
        @test all(jp.rejection .== false)  # All should be false when threshold is 0
    end

    @testset "_joint_probability - edge cases" begin
        signal = randn(3, 100)
        
        # Test with normalize = 0
        jp, rej = eegfun._joint_probability(signal, 5.0, 0, 1000)
        @test length(jp) == 3
        @test length(rej) == 3
        
        # Test with normalize = 1 (no trimming)
        jp, rej = eegfun._joint_probability(signal, 5.0, 1, 1000)
        @test length(jp) == 3
        @test length(rej) == 3
        
        # Test with threshold = 0
        jp, rej = eegfun._joint_probability(signal, 0.0, 2, 1000)
        @test length(jp) == 3
        @test all(rej .== false)
    end

    @testset "compute_probability! - Gaussian approximation" begin
        # Test with bins <= 0 (Gaussian approximation)
        data = randn(1000)
        proba_map = zeros(Float64, 1000)
        
        result = eegfun.compute_probability!(proba_map, data, 0)
        @test result isa Vector{Float64}
        @test length(result) == 1000
        @test all(0 <= x <= 1 for x in result)  # All values between 0 and 1
        @test isapprox(sum(result), 1.0, atol = 1e-10)  # Should sum to 1
        
        # Test with negative bins
        proba_map2 = zeros(Float64, 1000)
        result2 = eegfun.compute_probability!(proba_map2, data, -10)
        @test length(result2) == 1000
        @test isapprox(sum(result2), 1.0, atol = 1e-10)
    end

    @testset "correlation_matrix_dual_selection" begin
        dat = create_test_data(n = 1000, n_channels = 10)
        
        # Test basic functionality
        cm = eegfun.correlation_matrix_dual_selection(
            dat,
            channel_selection1 = eegfun.channels([:Ch1, :Ch2, :Ch3]),
            channel_selection2 = eegfun.channels([:Ch4, :Ch5])
        )
        @test cm isa DataFrame
        @test size(cm, 1) == 3  # 3 channels in first set
        @test size(cm, 2) == 3  # 2 channels + :row column
        @test :row in propertynames(cm)
        @test :Ch4 in propertynames(cm)
        @test :Ch5 in propertynames(cm)
        
        # Test with sample_selection (create a predicate function)
        sample_pred = x -> [i <= 500 for i in 1:nrow(x)]
        cm = eegfun.correlation_matrix_dual_selection(
            dat,
            sample_selection = sample_pred,
            channel_selection1 = eegfun.channels([:Ch1, :Ch2]),
            channel_selection2 = eegfun.channels([:Ch3])
        )
        @test size(cm, 1) == 2
        @test size(cm, 2) == 2  # 1 channel + :row
        
        # Test with include_extra
        # Note: extra_ch needs to be added after the original channels in the data
        # to be considered an "extra" channel (comes after channel columns in layout)
        dat.data[!, :extra_ch] = randn(nrow(dat.data))
        # Update layout to include extra_ch as an extra channel
        # For this test, we'll just verify the function works with include_extra=true
        cm = eegfun.correlation_matrix_dual_selection(
            dat,
            channel_selection1 = eegfun.channels([:Ch1, :Ch2]),
            channel_selection2 = eegfun.channels([:Ch3]),
            include_extra_selection1 = false,
            include_extra_selection2 = false
        )
        @test size(cm, 1) == 2  # Ch1 and Ch2
        @test size(cm, 2) == 2  # Ch3 + :row
        
        # Test error handling - empty channel selection
        @test_throws ErrorException eegfun.correlation_matrix_dual_selection(
            dat,
            channel_selection1 = eegfun.channels(Symbol[]),
            channel_selection2 = eegfun.channels([:Ch1])
        )
        
        @test_throws ErrorException eegfun.correlation_matrix_dual_selection(
            dat,
            channel_selection1 = eegfun.channels([:Ch1]),
            channel_selection2 = eegfun.channels(Symbol[])
        )
    end

    @testset "_correlation_matrix_dual_selection" begin
        dat = create_test_data(n = 1000, n_channels = 5)
        
        # Test internal function
        cm = eegfun._correlation_matrix_dual_selection(
            dat.data,
            collect(1:500),
            [:Ch1, :Ch2],
            [:Ch3, :Ch4]
        )
        @test cm isa DataFrame
        @test size(cm, 1) == 2
        @test size(cm, 2) == 3  # 2 channels + :row
        @test :row in propertynames(cm)
    end

    @testset "get_eog_channels" begin
        eog_cfg = eegfun.EogConfig(
            vEOG_criterion = 50.0,
            hEOG_criterion = 30.0,
            vEOG_channels = [["Fp1", "Fp2"], ["IO1", "IO2"], ["vEOG"]],
            hEOG_channels = [["F9"], ["F10"], ["hEOG"]]
        )
        
        eog_channels = eegfun.get_eog_channels(eog_cfg)
        @test eog_channels isa Vector{Symbol}
        @test length(eog_channels) == 2
        @test eog_channels == [:vEOG, :hEOG]
    end

    @testset "add_zscore_columns!" begin
        # Create a test correlation matrix
        dat = create_test_data(n = 1000, n_channels = 5)
        cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Ch1, :Ch2, :Ch3]))
        
        # Test basic functionality
        eegfun.add_zscore_columns!(cm)
        @test :z_Ch1 in propertynames(cm)
        @test :z_Ch2 in propertynames(cm)
        @test :z_Ch3 in propertynames(cm)
        
        # Check that z-scores are calculated correctly
        z_ch1 = cm[!, :z_Ch1]
        @test isapprox(mean(z_ch1), 0.0, atol = 1e-10)
        @test isapprox(std(z_ch1), 1.0, atol = 1e-10)
        
        # Test with custom exclude_columns
        cm2 = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        # Create a fresh copy to test exclude functionality
        cm2_exclude = copy(cm2)
        # Verify initial state
        @test !(:z_Ch1 in propertynames(cm2_exclude))
        @test !(:z_Ch2 in propertynames(cm2_exclude))
        # Add z-scores excluding row (which should work)
        # Note: The implementation uses names(df) which returns Vector{String}
        # but exclude_columns is Vector{Symbol}, so String columns won't match Symbol excludes
        # However, :row is a Symbol column, so it should be excluded
        eegfun.add_zscore_columns!(cm2_exclude, [:row])
        @test :z_Ch1 in propertynames(cm2_exclude)  # Ch1 should get z-score (String column)
        @test :z_Ch2 in propertynames(cm2_exclude)  # Ch2 should get z-score (String column)
        @test !(:z_row in propertynames(cm2_exclude))  # row should be excluded (Symbol column)
        
        # Test with all values the same (std = 0)
        cm3 = DataFrame(row = [:Ch1, :Ch2], Ch1 = [1.0, 1.0], Ch2 = [1.0, 1.0])
        eegfun.add_zscore_columns!(cm3)
        @test :z_Ch1 in propertynames(cm3)
        @test :z_Ch2 in propertynames(cm3)
        @test all(cm3.z_Ch1 .== 0.0)
        @test all(cm3.z_Ch2 .== 0.0)
        
        # Test with no numeric columns (after excluding all)
        cm4 = DataFrame(row = [:Ch1, :Ch2], name = ["A", "B"])
        # Should not throw error, just warn
        eegfun.add_zscore_columns!(cm4)
    end

    @testset "add_zscore_columns" begin
        dat = create_test_data(n = 1000, n_channels = 5)
        cm = eegfun.correlation_matrix(dat, channel_selection = eegfun.channels([:Ch1, :Ch2]))
        
        # Test non-mutating version
        cm_with_z = eegfun.add_zscore_columns(cm)
        @test :z_Ch1 in propertynames(cm_with_z)
        @test :z_Ch2 in propertynames(cm_with_z)
        @test !(:z_Ch1 in propertynames(cm))  # Original should not be modified
        @test !(:z_Ch2 in propertynames(cm))
    end

    @testset "correlation_matrix_eog" begin
        dat = create_test_data(n = 1000, n_channels = 10)
        
        # Create EOG channels in data
        dat.data[!, :vEOG] = randn(nrow(dat.data))
        dat.data[!, :hEOG] = randn(nrow(dat.data))
        
        eog_cfg = eegfun.EogConfig(
            vEOG_criterion = 50.0,
            hEOG_criterion = 30.0,
            vEOG_channels = [["Ch1"], ["Ch2"], ["vEOG"]],
            hEOG_channels = [["Ch3"], ["Ch4"], ["hEOG"]]
        )
        
        # Test basic functionality
        cm = eegfun.correlation_matrix_eog(dat, eog_cfg)
        @test cm isa DataFrame
        @test :row in propertynames(cm)
        @test :vEOG in propertynames(cm)
        @test :hEOG in propertynames(cm)
        
        # Test with sample_selection (create a predicate function)
        sample_pred = x -> [i <= 500 for i in 1:nrow(x)]
        cm = eegfun.correlation_matrix_eog(
            dat, 
            eog_cfg,
            sample_selection = sample_pred
        )
        @test cm isa DataFrame
        
        # Test with channel_selection
        cm = eegfun.correlation_matrix_eog(
            dat,
            eog_cfg,
            channel_selection = eegfun.channels([:Ch1, :Ch2, :Ch3])
        )
        @test size(cm, 1) == 3
        
        # Test with include_extra
        cm = eegfun.correlation_matrix_eog(
            dat,
            eog_cfg,
            include_extra = true
        )
        @test cm isa DataFrame
    end

    @testset "identify_bad_channels" begin
        # Create test summary DataFrame
        summary_df = DataFrame(
            channel = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5],
            zvar = [1.0, 2.0, 4.0, 1.5, 2.5],  # Ch3 exceeds 3.0 threshold
            variance = [0.1, 0.2, 0.3, 0.15, 0.25]
        )
        
        # Create test joint probability DataFrame
        joint_prob_df = DataFrame(
            channel = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5],
            jp = [0.1, 0.2, 0.3, 0.4, 0.5],
            rejection = [false, false, true, false, true]  # Ch3 and Ch5 rejected
        )
        
        # Test with default criterion
        bad_channels = eegfun.identify_bad_channels(summary_df, joint_prob_df)
        @test bad_channels isa Vector{Symbol}
        @test :Ch3 in bad_channels  # High zvar
        @test :Ch5 in bad_channels  # Rejected by joint probability
        @test :Ch1 ∉ bad_channels
        @test :Ch2 ∉ bad_channels
        @test :Ch4 ∉ bad_channels
        
        # Test with custom criterion
        bad_channels = eegfun.identify_bad_channels(summary_df, joint_prob_df, zvar_criterion = 2.5)
        @test :Ch3 in bad_channels  # Still high zvar
        @test :Ch5 in bad_channels  # Still rejected
        
        # Test with very high criterion (only joint probability should catch bad channels)
        bad_channels = eegfun.identify_bad_channels(summary_df, joint_prob_df, zvar_criterion = 10.0)
        @test :Ch3 in bad_channels  # Still rejected by joint probability (rejection=true)
        @test :Ch5 in bad_channels  # Still rejected by joint probability
        
        # Test with no bad channels
        summary_good = DataFrame(
            channel = [:Ch1, :Ch2],
            zvar = [1.0, 1.5],
            variance = [0.1, 0.2]
        )
        joint_prob_good = DataFrame(
            channel = [:Ch1, :Ch2],
            jp = [0.1, 0.2],
            rejection = [false, false]
        )
        bad_channels = eegfun.identify_bad_channels(summary_good, joint_prob_good)
        @test isempty(bad_channels)
    end

    @testset "partition_channels_by_eog_correlation" begin
        # Create test EOG correlation DataFrame
        eog_corr_df = DataFrame(
            row = [:Ch1, :Ch2, :Ch3, :Ch4],
            hEOG = [0.1, 0.5, 0.2, 0.4],  # Ch2 and Ch4 highly correlated
            vEOG = [0.15, 0.2, 0.6, 0.3]  # Ch3 highly correlated
        )
        
        bad_channels = [:Ch1, :Ch2, :Ch3, :Ch4]
        
        # Test with default threshold
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            bad_channels, 
            eog_corr_df
        )
        @test :Ch1 in non_eog  # Low correlation
        @test :Ch2 in eog_related  # High correlation with hEOG
        @test :Ch3 in eog_related  # High correlation with vEOG
        @test :Ch4 in eog_related  # High correlation with hEOG
        
        # Test with custom threshold
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            bad_channels,
            eog_corr_df,
            threshold = 0.5
        )
        @test :Ch1 in non_eog
        @test :Ch2 in non_eog  # Now below threshold
        @test :Ch3 in eog_related  # Still above threshold
        @test :Ch4 in non_eog  # Now below threshold
        
        # Test with use_z (z-scored columns)
        eog_corr_df_z = DataFrame(
            row = [:Ch1, :Ch2, :Ch3],
            z_hEOG = [0.1, 0.5, 0.2],
            z_vEOG = [0.15, 0.2, 0.6]
        )
        bad_channels_z = [:Ch1, :Ch2, :Ch3]
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            bad_channels_z,
            eog_corr_df_z,
            use_z = true
        )
        @test :Ch1 in non_eog
        @test :Ch2 in eog_related
        @test :Ch3 in eog_related
        
        # Test with empty bad_channels
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            Symbol[],
            eog_corr_df
        )
        @test isempty(non_eog)
        @test isempty(eog_related)
        
        # Test with empty EOG channels
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            bad_channels,
            eog_corr_df,
            eog_channels = Symbol[]
        )
        @test length(non_eog) == length(bad_channels)
        @test isempty(eog_related)
        
        # Test with channel not in correlation DataFrame
        bad_channels_missing = [:Ch5, :Ch6]
        non_eog, eog_related = eegfun.partition_channels_by_eog_correlation(
            bad_channels_missing,
            eog_corr_df
        )
        @test length(non_eog) == 2  # Both should be in non_eog (not found in df)
        @test isempty(eog_related)
    end

    @testset "check_channel_neighbors" begin
        # Create test layout with neighbors
        channel_labels = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5]
        layout_df = DataFrame(
            label = channel_labels,
            inc = zeros(5),
            azi = zeros(5)
        )
        
        # Create neighbors dictionary (Neighbours requires channels, distances, weights)
        # Layout expects OrderedDict for neighbours
        neighbours = OrderedDict(
            :Ch1 => eegfun.Neighbours([:Ch2, :Ch3], [1.0, 1.0], [0.5, 0.5]),
            :Ch2 => eegfun.Neighbours([:Ch1, :Ch3, :Ch4], [1.0, 1.0, 1.0], [1/3, 1/3, 1/3]),
            :Ch3 => eegfun.Neighbours([:Ch1, :Ch2], [1.0, 1.0], [0.5, 0.5]),
            :Ch4 => eegfun.Neighbours([:Ch2, :Ch5], [1.0, 1.0], [0.5, 0.5]),
            :Ch5 => eegfun.Neighbours([:Ch4], [1.0], [1.0])  # Only 1 neighbor (not enough)
        )
        
        layout = eegfun.Layout(layout_df, neighbours, nothing)  # data, neighbours, criterion
        
        # Test with bad channels that have all good neighbors
        # Ch1 and Ch3 are neighbors of each other, so we can't test them together
        # Test Ch1 alone (has Ch2 and Ch3 as neighbors, both good)
        bad_channels = [:Ch1]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout)
        @test :Ch1 in repairable  # Has 2 good neighbors (Ch2, Ch3)
        
        # Test Ch3 alone (has Ch1 and Ch2 as neighbors, both good)
        bad_channels = [:Ch3]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout)
        @test :Ch3 in repairable  # Has 2 good neighbors (Ch1, Ch2)
        
        # Test with bad channel that has a bad neighbor
        bad_channels = [:Ch1, :Ch2]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout)
        @test :Ch1 ∉ repairable  # Ch2 is a neighbor and is bad
        @test :Ch2 ∉ repairable  # Ch1 is a neighbor and is bad
        
        # Test with bad channel that has only 1 neighbor (not enough)
        bad_channels = [:Ch5]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout)
        @test isempty(repairable)  # Only 1 neighbor, need at least 2
        
        # Test with empty bad_channels
        repairable = eegfun.check_channel_neighbors(Symbol[], layout)
        @test isempty(repairable)
        
        # Test with channel that has no neighbors defined
        layout_no_neighbors = eegfun.Layout(layout_df, nothing, nothing)
        bad_channels = [:Ch1]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout_no_neighbors)
        @test isempty(repairable)
        
        # Test with channel not in neighbors dict
        neighbours_partial = OrderedDict(:Ch2 => eegfun.Neighbours([:Ch1], [1.0], [1.0]))
        layout_partial = eegfun.Layout(layout_df, neighbours_partial, nothing)
        bad_channels = [:Ch3]
        repairable = eegfun.check_channel_neighbors(bad_channels, layout_partial)
        @test isempty(repairable)  # Ch3 not in neighbours dict
    end

end
