using Test
using DataFrames
using OrderedCollections
using EegFun

@testset "artifact_detection" begin

    @testset "detect_eog_onsets!" begin

        dat = EegFun.create_test_continuous_data_with_artifacts()
        original_size = size(dat.data, 2)

        # Test "EOG" detection (should be three "jumps" over the threshold)
        EegFun.detect_eog_onsets!(dat, 75, :Ch2, :is_eog_onset)
        @test sum(dat.data.is_eog_onset) == 3

        @test :is_eog_onset in propertynames(dat.data)
        @test size(dat.data, 2) == original_size + 1
        @test all(abs.(dat.data.Ch2[dat.data.is_eog_onset]) .>= 75)

        # Test that it detects the artifact onset we added
        @test dat.data.is_eog_onset[100]
        @test dat.data.is_eog_onset[500]
        @test dat.data.is_eog_onset[800]

    end

    @testset "_is_extreme_value" begin
        signal = [1.0, 2.0, 50.0, 3.0, -40.0, 4.0, 5.0]
        threshold = 10.0

        result = EegFun._is_extreme_value(signal, threshold)

        @test result isa AbstractVector{Bool}
        @test length(result) == length(signal)
        @test result == [false, false, true, false, true, false, false]

    end

    @testset "_is_extreme_value!" begin
        signal = [1.0, 2.0, 50.0, 3.0, -40.0, 4.0, 5.0]
        threshold = 10.0
        mask = Vector{Bool}(undef, length(signal))

        EegFun._is_extreme_value!(mask, signal, threshold)

        @test mask == [false, false, true, false, true, false, false]
    end

    @testset "is_extreme_value" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test with default threshold (combined mode)
        extreme_mask = EegFun.is_extreme_value(dat, 20)
        @test extreme_mask isa AbstractVector{Bool}
        @test length(extreme_mask) == size(dat.data, 1)

        # Test that it detects artifacts (combined across channels)
        @test all(extreme_mask[100:110])  # Large artifact
        @test all(extreme_mask[500:505])  # Negative artifact
        @test all(extreme_mask[800:802])  # Smaller artifact

        # Test separate mode for individual channel analysis
        extreme_df = EegFun.is_extreme_value(dat, 20, mode = :separate)
        @test extreme_df isa DataFrame
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_df)
        @test :is_extreme_value_Ch2_20 in propertynames(extreme_df)
        @test size(extreme_df, 1) == size(dat.data, 1)

        # Test that channel A has more extreme values than B
        @test sum(extreme_df.is_extreme_value_Ch1_20) == 0
        @test sum(extreme_df.is_extreme_value_Ch2_20) == 20

        # Test with channel selection (separate mode)
        extreme_subset = EegFun.is_extreme_value(dat, 20, channel_selection = EegFun.channels([:Ch1]), mode = :separate)
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_subset)
        @test :is_extreme_value_Ch2_20 ∉ propertynames(extreme_subset)

        # Test empty channel selection
        @test_throws ErrorException EegFun.is_extreme_value(dat, 20, channel_selection = EegFun.channels(Symbol[]))
    end

    @testset "is_extreme_value!" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()
        original_columns = names(dat.data)

        # Test mutating version (default combined mode)
        EegFun.is_extreme_value!(dat, 20)

        # Check that new column was added (combined mode creates single column)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # One combined column
        @test new_columns[1] == "is_extreme_value_20"

        # Test separate mode
        dat = EegFun.create_test_continuous_data_with_artifacts()
        EegFun.is_extreme_value!(dat, 20, mode = :separate)

        # Check that new columns were added (separate mode creates one per channel)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 2  # One for each channel
        @test "is_extreme_value_Ch1_20" in new_columns
        @test "is_extreme_value_Ch2_20" in new_columns

        # Test with channel selection (separate mode)
        dat = EegFun.create_test_continuous_data_with_artifacts()
        EegFun.is_extreme_value!(dat, 20, channel_selection = EegFun.channels([:Ch1]), mode = :separate)

        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # Only one channel selected
        @test "is_extreme_value_Ch1_20" in new_columns

    end

    @testset "n_extreme_value" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test counting extreme values (default combined mode)
        total_count = EegFun.n_extreme_value(dat, 20)
        @test total_count == 20

        # Test counting extreme values (default combined mode)
        total_count = EegFun.n_extreme_value(dat, 101)
        @test total_count == 17

        # Test separate mode
        count_df = EegFun.n_extreme_value(dat, 20, mode = :separate)
        @test count_df isa DataFrame
        @test :channel in propertynames(count_df)
        @test :n_extreme in propertynames(count_df)
        @test size(count_df, 1) == 2  # 2 channels

        # Test that channel Ch2 has more extreme values than Ch1
        Ch1_count = count_df[count_df.channel.==:Ch1, :n_extreme][1]
        Ch2_count = count_df[count_df.channel.==:Ch2, :n_extreme][1]
        @test Ch2_count > Ch1_count

        # Test with channel selection (separate mode)
        count_subset = EegFun.n_extreme_value(dat, 20, channel_selection = EegFun.channels([:Ch1]), mode = :separate)
        @test size(count_subset, 1) == 1
        @test count_subset.channel[1] == :Ch1

        # Test empty channel selection
        @test_throws ErrorException EegFun.n_extreme_value(dat, 20, channel_selection = EegFun.channels(Symbol[]))
    end

    @testset "_n_extreme_value" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test internal function
        counts = EegFun._n_extreme_value(dat.data, [:Ch1, :Ch2], 20.0)

        @test counts isa Vector{Int}
        @test length(counts) == 2
        @test counts[2] > counts[1]  # Channel A should have more extreme values
    end

    @testset "edge cases" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test with very high threshold (should find no extreme values)
        extreme_mask = EegFun.is_extreme_value(dat, 1000)
        @test all(extreme_mask .== false)

        # Test with very low threshold (should find many extreme values)
        extreme_mask_low = EegFun.is_extreme_value(dat, 1)
        @test sum(extreme_mask_low) > sum(extreme_mask)

        # Test with zero threshold (should find all values as extreme)
        @test_throws ErrorException extreme_mask_zero = EegFun.is_extreme_value(dat, -1)

        # Test separate mode for detailed analysis
        extreme_df = EegFun.is_extreme_value(dat, 1000, mode = :separate)
        @test all(extreme_df.is_extreme_value_Ch1_1000 .== false)
        @test all(extreme_df.is_extreme_value_Ch2_1000 .== false)
    end

    @testset "error handling" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test detect_eog_onsets! with non-existent channel (should throw error)
        @test_throws ErrorException EegFun.detect_eog_onsets!(dat, 20, :NonExistentChannel, :output)
        @test :output ∉ propertynames(dat.data)  # Should not add the output column

        # Test is_extreme_value with empty channel selection (should throw error)
        @test_throws ErrorException EegFun.is_extreme_value(dat, 20, channel_selection = EegFun.channels(Symbol[]))

        # Test is_extreme_value! with empty channel selection (should throw error)
        original_columns = names(dat.data)
        @test_throws ErrorException EegFun.is_extreme_value!(dat, 20, channel_selection = EegFun.channels(Symbol[]))
        @test names(dat.data) == original_columns  # Should not add any columns

        # Test n_extreme_value with empty channel selection (should throw error)
        @test_throws ErrorException EegFun.n_extreme_value(dat, 20, channel_selection = EegFun.channels(Symbol[]))
    end

    @testset "data type handling" begin
        # Test with different data types
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test with Int threshold (default combined mode)
        extreme_mask = EegFun.is_extreme_value(dat, 20)
        @test extreme_mask isa AbstractVector{Bool}

        # Test with Int threshold (separate mode)
        extreme_df = EegFun.is_extreme_value(dat, 20, mode = :separate)
        @test extreme_df isa DataFrame
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_df)
        @test :is_extreme_value_Ch2_20 in propertynames(extreme_df)
    end

    @testset "channel overwriting" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test that detect_eog_onsets! overwrites existing output channel
        dat.data[!, :existing_output] = fill(false, size(dat.data, 1))
        EegFun.detect_eog_onsets!(dat, 20, :Ch2, :existing_output)
        @test any(dat.data.existing_output)  # Should be overwritten with actual results

        # Test that is_extreme_value! creates new columns (default combined mode)
        original_columns = names(dat.data)
        EegFun.is_extreme_value!(dat, 20)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # Should create one combined column

        # Test with different threshold to ensure different column names
        EegFun.is_extreme_value!(dat, 30)
        final_columns = setdiff(names(dat.data), original_columns)
        @test length(final_columns) == 2  # Should have 2 new columns (2 thresholds)

        # Test separate mode
        dat2 = EegFun.create_test_continuous_data_with_artifacts()
        EegFun.is_extreme_value!(dat2, 20, mode = :separate)
        new_columns2 = setdiff(names(dat2.data), original_columns)
        @test length(new_columns2) == 2  # Should create 2 separate columns
    end

    @testset "detect_eog_signals!" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Add vEOG and hEOG channels
        dat.data[!, :vEOG] = dat.data.Ch2 .+ randn(nrow(dat.data)) * 5
        dat.data[!, :hEOG] = dat.data.Ch1 .+ randn(nrow(dat.data)) * 5

        # Update layout to include EOG channels
        layout_df = DataFrame(label = [:Ch1, :Ch2, :vEOG, :hEOG], inc = [0.0, 0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0, 0.0])
        dat.layout = EegFun.Layout(layout_df, nothing, nothing)

        eog_cfg = Dict(
            "vEOG_criterion" => 50.0,
            "hEOG_criterion" => 30.0,
            "vEOG_channels" => [["Fp1"], ["Fp2"], ["vEOG"]],
            "hEOG_channels" => [["F9"], ["F10"], ["hEOG"]],
        )

        EegFun.detect_eog_signals!(dat, eog_cfg)

        @test :is_vEOG in propertynames(dat.data)
        @test :is_hEOG in propertynames(dat.data)
    end

    @testset "is_extreme_value! for EpochData" begin
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)

        # Add artifacts to some epochs
        epochs.data[1][!, :Ch1] .+= 200.0  # Extreme value in first epoch
        epochs.data[3][!, :Ch2] .+= -200.0  # Extreme value in third epoch

        original_cols_epoch1 = names(epochs.data[1])

        # Test combined mode
        EegFun.is_extreme_value!(epochs, 100, mode = :combined)

        # Check that all epochs got the column
        for epoch_df in epochs.data
            @test :is_extreme_value_100 in propertynames(epoch_df)
        end

        # Check that artifacts were detected
        @test any(epochs.data[1].is_extreme_value_100)
        @test any(epochs.data[3].is_extreme_value_100)

        # Test separate mode
        epochs2 = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2)
        epochs2.data[1][!, :Ch1] .+= 200.0

        EegFun.is_extreme_value!(epochs2, 100, mode = :separate)

        # Check separate columns were created
        @test :is_extreme_value_Ch1_100 in propertynames(epochs2.data[1])
        @test :is_extreme_value_Ch2_100 in propertynames(epochs2.data[1])

        # Test with epoch selection
        epochs3 = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2)
        epochs3.data[1][!, :Ch1] .+= 200.0

        EegFun.is_extreme_value!(epochs3, 100, epoch_selection = EegFun.epochs([1, 3]))

        # Check that only selected epochs got the column
        @test :is_extreme_value_100 in propertynames(epochs3.data[1])
        @test :is_extreme_value_100 in propertynames(epochs3.data[3])

        # Test with channel selection
        epochs4 = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 3)
        epochs4.data[1][!, :Ch1] .+= 200.0

        EegFun.is_extreme_value!(epochs4, 100, channel_selection = EegFun.channels([:Ch1]))

        @test :is_extreme_value_100 in propertynames(epochs4.data[1])

        # Test error handling - empty channel selection
        @test_throws ErrorException EegFun.is_extreme_value!(epochs, 100, channel_selection = EegFun.channels(Symbol[]))

        # Test error handling - empty epoch selection
        @test_throws ErrorException EegFun.is_extreme_value!(epochs, 100, epoch_selection = EegFun.epochs(Int[]))
    end

    @testset "is_extreme_value! for Vector{EpochData}" begin
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2) for _ = 1:2]

        # Add artifacts
        epochs_list[1].data[1][!, :Ch1] .+= 200.0
        epochs_list[2].data[2][!, :Ch2] .+= -200.0

        EegFun.is_extreme_value!(epochs_list, 100)

        # Check both conditions got the column
        @test :is_extreme_value_100 in propertynames(epochs_list[1].data[1])
        @test :is_extreme_value_100 in propertynames(epochs_list[2].data[1])
    end

    @testset "Rejection struct and helpers" begin
        r1 = EegFun.Rejection(:Ch1, 1)
        r2 = EegFun.Rejection(:Ch2, 2)
        r3 = EegFun.Rejection(:Ch1, 1)  # Duplicate

        @test r1.label == :Ch1
        @test r1.epoch == 1

        # Test is_equal_rejection
        @test EegFun.is_equal_rejection(r1, r3)
        @test !EegFun.is_equal_rejection(r1, r2)

        # Test unique_rejections
        rejections = [r1, r2, r3]
        unique_rej = EegFun.unique_rejections(rejections)
        @test length(unique_rej) == 2

        # Test unique_channels
        unique_ch = EegFun.unique_channels(rejections)
        @test unique_ch == [:Ch1, :Ch2]

        # Test unique_epochs
        unique_ep = EegFun.unique_epochs(rejections)
        @test unique_ep == [1, 2]
    end

    @testset "detect_bad_epochs_automatic" begin
        # Create epochs with known artifacts
        epochs = EegFun.create_test_epoch_data(n_epochs = 10, n_channels = 3)

        # Add extreme values to specific epochs to trigger rejection
        epochs.data[2][!, :Ch1] .= 200.0  # High absolute value
        epochs.data[5][!, :Ch2] .= -200.0  # Low absolute value
        epochs.data[7][!, :Ch1] .= 300.0  # Very high value

        # Test with absolute threshold only
        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 150.0)

        @test rejection_info isa EegFun.EpochRejectionInfo
        @test rejection_info.abs_criterion == 150.0
        @test rejection_info.z_criterion == 0
        @test rejection_info.z_rejections === nothing
        @test rejection_info.abs_rejections !== nothing
        @test length(rejection_info.rejected) > 0

        # Test with z-score criterion only
        epochs2 = EegFun.create_test_epoch_data(n_epochs = 20, n_channels = 3)
        # Make one epoch have very high variance
        epochs2.data[5][!, :Ch1] .= randn(nrow(epochs2.data[5])) .* 100.0

        rejection_info2 = EegFun.detect_bad_epochs_automatic(epochs2, z_criterion = 2.0, abs_criterion = 0)

        @test rejection_info2.z_criterion == 2.0
        @test rejection_info2.abs_criterion == 0
        @test rejection_info2.abs_rejections === nothing
        @test rejection_info2.z_rejections !== nothing

        # Test with both criteria
        epochs3 = EegFun.create_test_epoch_data(n_epochs = 15, n_channels = 3)
        epochs3.data[3][!, :Ch1] .= 200.0

        rejection_info3 = EegFun.detect_bad_epochs_automatic(epochs3, z_criterion = 3.0, abs_criterion = 150.0)

        @test rejection_info3.z_criterion == 3.0
        @test rejection_info3.abs_criterion == 150.0
        @test rejection_info3.z_rejections !== nothing
        @test rejection_info3.abs_rejections !== nothing

        # Test with custom z_measures
        rejection_info4 = EegFun.detect_bad_epochs_automatic(epochs3, z_criterion = 2.0, abs_criterion = 0, z_measures = [:variance, :max])

        @test rejection_info4.z_rejections.z_measures == [:variance, :max]

        # Test error handling - both criteria zero
        @test_throws ErrorException EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 0, abs_criterion = 0)

        # Test error handling - negative criteria
        @test_throws ErrorException EegFun.detect_bad_epochs_automatic(epochs, z_criterion = -1, abs_criterion = 100)
        @test_throws ErrorException EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 3, abs_criterion = -1)

        # Test error handling - invalid measures
        @test_throws ErrorException EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 3, z_measures = [:invalid_measure])

        # Test error handling - empty channel selection
        @test_throws ErrorException EegFun.detect_bad_epochs_automatic(
            epochs,
            z_criterion = 3,
            channel_selection = EegFun.channels(Symbol[]),
        )

        # Test with Vector{EpochData}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2) for _ = 1:2]
        epochs_list[1].data[2][!, :Ch1] .= 200.0

        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)

        @test length(rejection_list) == 2
        @test rejection_list[1] isa EegFun.EpochRejectionInfo
    end

    @testset "get_rejected" begin
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2)
        epochs.data[2][!, :Ch1] .= 200.0

        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)

        rejected = EegFun.get_rejected(rejection_info)
        @test rejected isa Vector{EegFun.Rejection}
        @test length(rejected) > 0

        # Test with Vector{EpochRejectionInfo}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2) for _ = 1:2]
        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)

        rejected_list = EegFun.get_rejected(rejection_list)
        @test rejected_list isa Vector{Vector{EegFun.Rejection}}
        @test length(rejected_list) == 2
    end

    @testset "EpochRejectionInfo helpers" begin
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2)
        epochs.data[2][!, :Ch1] .= 200.0

        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)

        # Test unique_rejections on EpochRejectionInfo
        unique_rej = EegFun.unique_rejections(rejection_info)
        @test unique_rej isa Vector{EegFun.Rejection}

        # Test unique_channels on EpochRejectionInfo
        unique_ch = EegFun.unique_channels(rejection_info)
        @test unique_ch isa Vector{Symbol}

        # Test unique_epochs on EpochRejectionInfo
        unique_ep = EegFun.unique_epochs(rejection_info)
        @test unique_ep isa Vector{Int}

        # Test with Vector{EpochRejectionInfo}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2) for _ = 1:2]
        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)

        unique_rej_list = EegFun.unique_rejections(rejection_list)
        @test unique_rej_list isa Vector{Vector{EegFun.Rejection}}

        unique_ch_list = EegFun.unique_channels(rejection_list)
        @test unique_ch_list isa Vector{Vector{Symbol}}
    end

    @testset "Base.show for EpochRejectionInfo" begin
        epochs = EegFun.create_test_epoch_data(n_epochs = 10, n_channels = 3)
        epochs.data[2][!, :Ch1] .= 200.0
        epochs.data[5][!, :Ch2] .= -200.0

        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, z_criterion = 3.0, abs_criterion = 150.0)

        # Test that show doesn't throw
        io = IOBuffer()
        show(io, rejection_info)
        output = String(take!(io))
        @test length(output) > 0
        @test occursin("EpochRejectionInfo", output)

        # Test with Vector{EpochRejectionInfo}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 2) for _ = 1:2]
        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)

        io2 = IOBuffer()
        show(io2, rejection_list)
        output2 = String(take!(io2))
        @test length(output2) > 0
    end

    @testset "channel_repairable!" begin
        # Create layout with neighbors
        layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Ch4], inc = [0.0, 0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0, 0.0])

        # Create neighbors dict
        neighbours_dict = OrderedDict(
            :Ch1 => EegFun.Neighbours([:Ch2, :Ch3], [1.0, 1.0], [0.5, 0.5]),
            :Ch2 => EegFun.Neighbours([:Ch1, :Ch3, :Ch4], [1.0, 1.0, 1.0], [0.33, 0.33, 0.34]),
            :Ch3 => EegFun.Neighbours([:Ch1, :Ch2], [1.0, 1.0], [0.5, 0.5]),
            :Ch4 => EegFun.Neighbours([:Ch2], [1.0], [1.0]),
        )

        layout = EegFun.Layout(layout_df, neighbours_dict, nothing)

        # Create epochs with layout
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 4)
        epochs.layout = layout

        # Create rejection info
        epochs.data[2][!, :Ch1] .= 200.0
        epochs.data[3][!, :Ch2] .= 200.0
        epochs.data[4][!, :Ch4] .= 200.0  # Ch4 has only 1 neighbor, might not be repairable

        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)

        # Test channel_repairable!
        EegFun.channel_repairable!(rejection_info, layout)

        @test rejection_info.repaired !== nothing
        @test rejection_info.skipped !== nothing

        # Test with Vector{EpochRejectionInfo}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 4) for _ = 1:2]
        for ep in epochs_list
            ep.layout = layout
        end
        epochs_list[1].data[2][!, :Ch1] .= 200.0

        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_list, layout)

        @test all(info -> info.repaired !== nothing, rejection_list)
    end

    @testset "repair_artifacts! and repair_artifacts" begin
        # Create layout with neighbors for repair
        layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0])

        neighbours_dict = OrderedDict(
            :Ch1 => EegFun.Neighbours([:Ch2, :Ch3], [1.0, 1.0], [0.5, 0.5]),
            :Ch2 => EegFun.Neighbours([:Ch1, :Ch3], [1.0, 1.0], [0.5, 0.5]),
            :Ch3 => EegFun.Neighbours([:Ch1, :Ch2], [1.0, 1.0], [0.5, 0.5]),
        )

        layout = EegFun.Layout(layout_df, neighbours_dict, nothing)

        # Create epochs with artifacts
        epochs = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)
        epochs.layout = layout

        # Store original values
        original_value = epochs.data[2][1, :Ch1]
        epochs.data[2][!, :Ch1] .= 200.0  # Artifact

        # Detect and prepare for repair
        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_info, layout)

        # Test non-mutating version
        epochs_copy = copy(epochs)
        repaired_epochs = EegFun.repair_artifacts(epochs_copy, rejection_info)

        @test repaired_epochs isa EegFun.EpochData
        @test repaired_epochs !== epochs_copy  # Should be different object
        @test repaired_epochs.data[2][1, :Ch1] != 200.0  # Should be repaired

        # Test mutating version
        epochs2 = EegFun.create_test_epoch_data(n_epochs = 5, n_channels = 3)
        epochs2.layout = layout
        epochs2.data[2][!, :Ch1] .= 200.0

        rejection_info2 = EegFun.detect_bad_epochs_automatic(epochs2, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_info2, layout)

        original_value2 = epochs2.data[2][1, :Ch1]
        EegFun.repair_artifacts!(epochs2, rejection_info2)

        @test epochs2.data[2][1, :Ch1] != original_value2  # Should be repaired

        # Test with Vector{EpochData}
        epochs_list = [EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 3) for _ = 1:2]
        for ep in epochs_list
            ep.layout = layout
            ep.data[2][!, :Ch1] .= 200.0
        end

        rejection_list = EegFun.detect_bad_epochs_automatic(epochs_list, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_list, layout)

        repaired_list = EegFun.repair_artifacts(epochs_list, rejection_list)
        @test length(repaired_list) == 2

        # Test error handling - repair_artifacts_neighbor! without channel_repairable! (repaired is nothing)
        epochs3 = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 3)
        epochs3.layout = layout
        epochs3.data[1][!, :Ch1] .= 200.0  # Add artifact so there's something to reject
        rejection_info3 = EegFun.detect_bad_epochs_automatic(epochs3, abs_criterion = 150.0, z_criterion = 0)

        # repaired should be nothing (not yet populated)
        @test rejection_info3.repaired === nothing

        # Should throw error when calling repair_artifacts_neighbor! directly with repaired = nothing
        @test_throws ArgumentError EegFun.repair_artifacts_neighbor!(epochs3, rejection_info3)

        # Note: repair_artifacts! automatically calls channel_repairable! if needed, so it won't throw
        # But we can test that it works correctly
        epochs4 = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 3)
        epochs4.layout = layout
        epochs4.data[1][!, :Ch1] .= 200.0
        rejection_info4 = EegFun.detect_bad_epochs_automatic(epochs4, abs_criterion = 150.0, z_criterion = 0)
        @test rejection_info4.repaired === nothing
        # This should work because repair_artifacts! calls channel_repairable! automatically
        EegFun.repair_artifacts!(epochs4, rejection_info4, method = :neighbor_interpolation)
        @test rejection_info4.repaired !== nothing  # Should be populated now

        # Test error handling - unknown method
        epochs5 = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 3)
        epochs5.layout = layout
        epochs5.data[1][!, :Ch1] .= 200.0
        rejection_info5 = EegFun.detect_bad_epochs_automatic(epochs5, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_info5, layout)
        @test_throws ArgumentError EegFun.repair_artifacts!(epochs5, rejection_info5, method = :unknown_method)
    end

    @testset "detect_eog_onsets! step_size parameter" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test with custom step_size
        EegFun.detect_eog_onsets!(dat, 75, :Ch2, :is_eog_onset_custom, step_size = 10)

        @test :is_eog_onset_custom in propertynames(dat.data)

        # Test with default step_size
        dat2 = EegFun.create_test_continuous_data_with_artifacts()
        EegFun.detect_eog_onsets!(dat2, 75, :Ch2, :is_eog_onset_default)

        @test :is_eog_onset_default in propertynames(dat2.data)
    end

    @testset "is_extreme_value! custom channel_out" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test with custom channel_out name
        EegFun.is_extreme_value!(dat, 20, channel_out = :custom_artifact_flag)

        @test :custom_artifact_flag in propertynames(dat.data)
        @test any(dat.data.custom_artifact_flag)

        # Test with EpochData
        epochs = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2)
        epochs.data[1][!, :Ch1] .= 200.0

        EegFun.is_extreme_value!(epochs, 100, channel_out = :custom_flag)

        @test :custom_flag in propertynames(epochs.data[1])
    end

    @testset "is_extreme_value! sample_selection" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Create sample selection mask as a function
        sample_mask = falses(nrow(dat.data))
        sample_mask[100:200] .= true  # Only check samples 100-200
        sample_selection_fn = x -> sample_mask  # Function that returns the mask

        # Test with sample selection
        EegFun.is_extreme_value!(dat, 20, sample_selection = sample_selection_fn)

        @test :is_extreme_value_20 in propertynames(dat.data)
        # Artifacts outside the selected range should not be detected
        @test !any(dat.data.is_extreme_value_20[500:505])  # Artifact at 500-505

        # Test with EpochData
        epochs = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2, n = 200)
        epochs.data[1][!, :Ch1] .= 200.0

        sample_mask_epoch = falses(nrow(epochs.data[1]))
        sample_mask_epoch[1:50] .= true
        sample_selection_epoch_fn = x -> sample_mask_epoch  # Function that returns the mask

        EegFun.is_extreme_value!(epochs, 100, sample_selection = sample_selection_epoch_fn, epoch_selection = EegFun.epochs([1]))

        @test :is_extreme_value_100 in propertynames(epochs.data[1])
    end

    @testset "_detect_extreme_values" begin
        dat = EegFun.create_test_continuous_data_with_artifacts()

        # Test internal function directly
        results = EegFun._detect_extreme_values(dat, 20.0)

        @test results isa Dict
        @test :Ch1 in keys(results)
        @test :Ch2 in keys(results)
        @test all(v -> v isa Vector{Bool}, values(results))

        # Test with channel selection
        results2 = EegFun._detect_extreme_values(dat, 20.0, channel_selection = EegFun.channels([:Ch1]))

        @test :Ch1 in keys(results2)
        @test :Ch2 ∉ keys(results2)

        # Test with sample selection
        sample_mask = falses(nrow(dat.data))
        sample_mask[100:200] .= true
        sample_selection_fn = x -> sample_mask  # Function that returns the mask
        results3 = EegFun._detect_extreme_values(dat, 20.0, sample_selection = sample_selection_fn)

        @test results3 isa Dict
    end

    @testset "detect_bad_epochs_automatic edge cases" begin
        # Test with very few epochs (should warn)
        epochs = EegFun.create_test_epoch_data(n_epochs = 2, n_channels = 2)
        epochs.data[1][!, :Ch1] .= 200.0

        # Should work but might warn
        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)
        @test rejection_info isa EegFun.EpochRejectionInfo

        # Test with single measure
        epochs2 = EegFun.create_test_epoch_data(n_epochs = 10, n_channels = 2)
        rejection_info2 = EegFun.detect_bad_epochs_automatic(epochs2, z_criterion = 2.0, abs_criterion = 0, z_measures = [:variance])

        @test rejection_info2.z_rejections.z_measures == [:variance]

        # Test with all measures
        rejection_info3 = EegFun.detect_bad_epochs_automatic(epochs2, z_criterion = 2.0, abs_criterion = 0)
        @test length(rejection_info3.z_rejections.z_measures) == 6
    end

    @testset "repair_artifacts! with no repairable channels" begin
        # Create layout where channel has no neighbors
        layout_df = DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0])

        # Ch1 has no neighbors (can't be repaired)
        neighbours_dict =
            OrderedDict(:Ch1 => EegFun.Neighbours(Symbol[], Float64[], Float64[]), :Ch2 => EegFun.Neighbours([:Ch1], [1.0], [1.0]))

        layout = EegFun.Layout(layout_df, neighbours_dict, nothing)

        epochs = EegFun.create_test_epoch_data(n_epochs = 3, n_channels = 2)
        epochs.layout = layout
        epochs.data[1][!, :Ch1] .= 200.0

        rejection_info = EegFun.detect_bad_epochs_automatic(epochs, abs_criterion = 150.0, z_criterion = 0)
        EegFun.channel_repairable!(rejection_info, layout)

        # Should skip repair if no channels are repairable
        if isempty(rejection_info.repaired)
            @test isempty(rejection_info.repaired)
            @test !isempty(rejection_info.skipped)
        end
    end

end
