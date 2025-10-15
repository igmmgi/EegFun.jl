using Test
using DataFrames
using eegfun

# Helper function for artifact detection testing
function create_test_data_with_artifacts(; n::Int = 1000, fs::Int = 1000)

    t = collect(0:(n-1)) ./ fs

    # Create clean signal with some artifacts
    clean_signal = sin.(2π .* 10 .* t) .* 20
    # Add some extreme values (artifacts)
    artifact_signal = copy(clean_signal)
    artifact_signal[100:110] .= 200.0  # Large positive artifact
    artifact_signal[500:505] .= -200.0  # Large negative artifact
    artifact_signal[800:802] .= 100.0  # Smaller positive artifact

    df = DataFrame(:time => t, :triggers => zeros(Int, n), :Ch1 => clean_signal, :Ch2 => artifact_signal)
    layout = eegfun.Layout(DataFrame(label = [:Ch1, :Ch2], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

    dat = eegfun.ContinuousData(df, layout, fs, eegfun.AnalysisInfo())

    return dat

end


@testset "artifact_detection" begin
    @testset "detect_eog_onsets!" begin

        dat = create_test_data_with_artifacts()
        original_size = size(dat.data, 2)

        # Test "EOG" detection (should be three "jumps" over the threshold)
        eegfun.detect_eog_onsets!(dat, 75, :Ch2, :is_eog_onset)
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

        result = eegfun._is_extreme_value(signal, threshold)

        @test result isa AbstractVector{Bool}
        @test length(result) == length(signal)
        @test result == [false, false, true, false, true, false, false]

    end

    @testset "_is_extreme_value!" begin
        signal = [1.0, 2.0, 50.0, 3.0, -40.0, 4.0, 5.0]
        threshold = 10.0
        mask = Vector{Bool}(undef, length(signal))

        eegfun._is_extreme_value!(mask, signal, threshold)

        @test mask == [false, false, true, false, true, false, false]
    end

    @testset "is_extreme_value" begin
        dat = create_test_data_with_artifacts()

        # Test with default threshold (combined mode)
        extreme_mask = eegfun.is_extreme_value(dat, 20)
        @test extreme_mask isa AbstractVector{Bool}
        @test length(extreme_mask) == size(dat.data, 1)

        # Test that it detects artifacts (combined across channels)
        @test all(extreme_mask[100:110])  # Large artifact
        @test all(extreme_mask[500:505])  # Negative artifact
        @test all(extreme_mask[800:802])  # Smaller artifact

        # Test separate mode for individual channel analysis
        extreme_df = eegfun.is_extreme_value(dat, 20, mode = :separate)
        @test extreme_df isa DataFrame
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_df)
        @test :is_extreme_value_Ch2_20 in propertynames(extreme_df)
        @test size(extreme_df, 1) == size(dat.data, 1)

        # Test that channel A has more extreme values than B
        @test sum(extreme_df.is_extreme_value_Ch1_20) == 0
        @test sum(extreme_df.is_extreme_value_Ch2_20) == 20

        # Test with channel selection (separate mode)
        extreme_subset = eegfun.is_extreme_value(dat, 20, channel_selection = eegfun.channels([:Ch1]), mode = :separate)
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_subset)
        @test :is_extreme_value_Ch2_20 ∉ propertynames(extreme_subset)

        # Test empty channel selection
        @test_throws ErrorException eegfun.is_extreme_value(dat, 20, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "is_extreme_value!" begin
        dat = create_test_data_with_artifacts()
        original_columns = names(dat.data)

        # Test mutating version (default combined mode)
        eegfun.is_extreme_value!(dat, 20)

        # Check that new column was added (combined mode creates single column)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # One combined column
        @test new_columns[1] == "is_extreme_value_20"

        # Test separate mode
        dat = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat, 20, mode = :separate)

        # Check that new columns were added (separate mode creates one per channel)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 2  # One for each channel
        @test "is_extreme_value_Ch1_20" in new_columns
        @test "is_extreme_value_Ch2_20" in new_columns

        # Test with channel selection (separate mode)
        dat = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat, 20, channel_selection = eegfun.channels([:Ch1]), mode = :separate)

        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # Only one channel selected
        @test "is_extreme_value_Ch1_20" in new_columns

    end

    @testset "n_extreme_value" begin
        dat = create_test_data_with_artifacts()

        # Test counting extreme values (default combined mode)
        total_count = eegfun.n_extreme_value(dat, 20)
        @test total_count == 20

        # Test counting extreme values (default combined mode)
        total_count = eegfun.n_extreme_value(dat, 101)
        @test total_count == 17 

        # Test separate mode
        count_df = eegfun.n_extreme_value(dat, 20, mode = :separate)
        @test count_df isa DataFrame
        @test :channel in propertynames(count_df)
        @test :n_extreme in propertynames(count_df)
        @test size(count_df, 1) == 2  # 2 channels

        # Test that channel Ch2 has more extreme values than Ch1
        Ch1_count = count_df[count_df.channel .== :Ch1, :n_extreme][1]
        Ch2_count = count_df[count_df.channel .== :Ch2, :n_extreme][1]
        @test Ch2_count > Ch1_count

        # Test with channel selection (separate mode)
        count_subset = eegfun.n_extreme_value(dat, 20, channel_selection = eegfun.channels([:Ch1]), mode = :separate)
        @test size(count_subset, 1) == 1
        @test count_subset.channel[1] == :Ch1

        # Test empty channel selection
        @test_throws ErrorException eegfun.n_extreme_value(dat, 20, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "_n_extreme_value" begin
        dat = create_test_data_with_artifacts()

        # Test internal function
        counts = eegfun._n_extreme_value(dat.data, [:Ch1, :Ch2], 20.0)

        @test counts isa Vector{Int}
        @test length(counts) == 2
        @test counts[2] > counts[1]  # Channel A should have more extreme values
    end

    @testset "edge cases" begin
        dat = create_test_data_with_artifacts()

        # Test with very high threshold (should find no extreme values)
        extreme_mask = eegfun.is_extreme_value(dat, 1000)
        @test all(extreme_mask .== false)

        # Test with very low threshold (should find many extreme values)
        extreme_mask_low = eegfun.is_extreme_value(dat, 1)
        @test sum(extreme_mask_low) > sum(extreme_mask)

        # Test with zero threshold (should find all values as extreme)
        @test_throws ErrorException extreme_mask_zero = eegfun.is_extreme_value(dat, -1)

        # Test separate mode for detailed analysis
        extreme_df = eegfun.is_extreme_value(dat, 1000, mode = :separate)
        @test all(extreme_df.is_extreme_value_Ch1_1000 .== false)
        @test all(extreme_df.is_extreme_value_Ch2_1000 .== false)
    end

    @testset "error handling" begin
        dat = create_test_data_with_artifacts()

        # Test detect_eog_onsets! with non-existent channel (should throw error)
        @test_throws ErrorException eegfun.detect_eog_onsets!(dat, 20, :NonExistentChannel, :output)
        @test :output ∉ propertynames(dat.data)  # Should not add the output column

        # Test is_extreme_value with empty channel selection (should throw error)
        @test_throws ErrorException eegfun.is_extreme_value(dat, 20, channel_selection = eegfun.channels(Symbol[]))

        # Test is_extreme_value! with empty channel selection (should throw error)
        original_columns = names(dat.data)
        @test_throws ErrorException eegfun.is_extreme_value!(dat, 20, channel_selection = eegfun.channels(Symbol[]))
        @test names(dat.data) == original_columns  # Should not add any columns

        # Test n_extreme_value with empty channel selection (should throw error)
        @test_throws ErrorException eegfun.n_extreme_value(dat, 20, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "data type handling" begin
        # Test with different data types
        dat = create_test_data_with_artifacts()

        # Test with Int threshold (default combined mode)
        extreme_mask = eegfun.is_extreme_value(dat, 20)
        @test extreme_mask isa AbstractVector{Bool}

        # Test with Int threshold (separate mode)
        extreme_df = eegfun.is_extreme_value(dat, 20, mode = :separate)
        @test extreme_df isa DataFrame
        @test :is_extreme_value_Ch1_20 in propertynames(extreme_df)
        @test :is_extreme_value_Ch2_20 in propertynames(extreme_df)
    end

    @testset "channel overwriting" begin
        dat = create_test_data_with_artifacts()

        # Test that detect_eog_onsets! overwrites existing output channel
        dat.data[!, :existing_output] = fill(false, size(dat.data, 1))
        eegfun.detect_eog_onsets!(dat, 20, :Ch2, :existing_output)
        @test any(dat.data.existing_output)  # Should be overwritten with actual results

        # Test that is_extreme_value! creates new columns (default combined mode)
        original_columns = names(dat.data)
        eegfun.is_extreme_value!(dat, 20)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 1  # Should create one combined column

        # Test with different threshold to ensure different column names
        eegfun.is_extreme_value!(dat, 30)
        final_columns = setdiff(names(dat.data), original_columns)
        @test length(final_columns) == 2  # Should have 2 new columns (2 thresholds)

        # Test separate mode
        dat2 = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat2, 20, mode = :separate)
        new_columns2 = setdiff(names(dat2.data), original_columns)
        @test length(new_columns2) == 2  # Should create 2 separate columns
    end

end
