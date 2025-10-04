using Test
using DataFrames
using eegfun

@testset "artifact_detection" begin

    # Helper to create test data with artifacts
    function create_test_data_with_artifacts(; n::Int = 1000, fs::Int = 1000)
        t = collect(0:(n-1)) ./ fs
        # Create clean signal with some artifacts
        clean_signal = sin.(2π .* 10 .* t) .+ 0.1 .* randn(length(t))

        # Add some extreme values (artifacts)
        signal_with_artifacts = copy(clean_signal)
        signal_with_artifacts[100:110] .= 50.0  # Large artifact
        signal_with_artifacts[500:505] .= -40.0  # Negative artifact
        signal_with_artifacts[800:802] .= 30.0  # Smaller artifact

        df = DataFrame()
        df[!, :time] = t
        df[!, :triggers] = zeros(Int, n)
        df[!, :A] = signal_with_artifacts
        df[!, :B] = clean_signal .+ 0.1 .* randn(length(t))
        layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        dat = eegfun.ContinuousData(df, layout, fs, eegfun.AnalysisInfo())
        return dat
    end

    @testset "detect_eog_onsets!" begin
        dat = create_test_data_with_artifacts()
        original_size = size(dat.data, 2)


        # Test EOG detection
        eegfun.detect_eog_onsets!(dat, 20, :A, :is_eog_onset)

        @test :is_eog_onset in propertynames(dat.data)
        @test size(dat.data, 2) == original_size + 1
        @test all(dat.data.is_eog_onset .== (abs.(dat.data.A) .> 20.0))

        # Test that it detects the artifacts we added
        @test any(dat.data.is_eog_onset[100:110])  # Large artifact
        @test any(dat.data.is_eog_onset[500:505])  # Negative artifact
        @test any(dat.data.is_eog_onset[800:802])  # Smaller artifact
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
        @test any(extreme_mask[100:110])  # Large artifact
        @test any(extreme_mask[500:505])  # Negative artifact
        @test any(extreme_mask[800:802])  # Smaller artifact

        # Test separate mode for individual channel analysis
        extreme_df = eegfun.is_extreme_value(dat, 20, mode = :separate)
        @test extreme_df isa DataFrame
        @test :is_extreme_value_A_20 in propertynames(extreme_df)
        @test :is_extreme_value_B_20 in propertynames(extreme_df)
        @test size(extreme_df, 1) == size(dat.data, 1)

        # Test that channel A has more extreme values than B
        @test sum(extreme_df.is_extreme_value_A_20) > sum(extreme_df.is_extreme_value_B_20)

        # Test with channel selection (separate mode)
        extreme_subset = eegfun.is_extreme_value(dat, 20, channel_selection = eegfun.channels([:A]), mode = :separate)
        @test :is_extreme_value_A_20 in propertynames(extreme_subset)
        @test :is_extreme_value_B_20 ∉ propertynames(extreme_subset)

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
        @test any(occursin("is_extreme_value_20", string(col)) for col in new_columns)

        # Test separate mode
        dat2 = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat2, 20, mode = :separate)

        # Check that new columns were added (separate mode creates one per channel)
        new_columns2 = setdiff(names(dat2.data), original_columns)
        @test length(new_columns2) == 2  # One for each channel
        @test any(occursin("is_extreme_value_A_20", string(col)) for col in new_columns2)
        @test any(occursin("is_extreme_value_B_20", string(col)) for col in new_columns2)

        # Test with channel selection (separate mode)
        dat3 = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat3, 20, channel_selection = eegfun.channels([:A]), mode = :separate)

        new_columns3 = setdiff(names(dat3.data), original_columns)
        @test length(new_columns3) == 1  # Only one channel selected
    end

    @testset "n_extreme_value" begin
        dat = create_test_data_with_artifacts()

        # Test counting extreme values (default combined mode)
        total_count = eegfun.n_extreme_value(dat, 20)
        @test total_count isa Int
        @test total_count > 0

        # Test separate mode
        count_df = eegfun.n_extreme_value(dat, 20, mode = :separate)
        @test count_df isa DataFrame
        @test :channel in propertynames(count_df)
        @test :n_extreme in propertynames(count_df)
        @test :threshold in propertynames(count_df)
        @test size(count_df, 1) == 2  # 2 channels

        # Test that channel A has more extreme values than B
        a_count = count_df[count_df.channel .== :A, :n_extreme][1]
        b_count = count_df[count_df.channel .== :B, :n_extreme][1]
        @test a_count > b_count

        # Test threshold is correct
        @test all(count_df.threshold .== 20)

        # Test with channel selection (separate mode)
        count_subset = eegfun.n_extreme_value(dat, 20, channel_selection = eegfun.channels([:A]), mode = :separate)
        @test size(count_subset, 1) == 1
        @test count_subset.channel[1] == :A

        # Test empty channel selection
        @test_throws ErrorException eegfun.n_extreme_value(dat, 20, channel_selection = eegfun.channels(Symbol[]))
    end

    @testset "_n_extreme_value" begin
        dat = create_test_data_with_artifacts()
        channels = [:A, :B]

        # Test internal function
        counts = eegfun._n_extreme_value(dat.data, channels, 20.0)

        @test counts isa Vector{Int}
        @test length(counts) == 2
        @test counts[1] > counts[2]  # Channel A should have more extreme values
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
        extreme_mask_zero = eegfun.is_extreme_value(dat, 0)
        @test all(extreme_mask_zero .== true)

        # Test separate mode for detailed analysis
        extreme_df = eegfun.is_extreme_value(dat, 1000, mode = :separate)
        @test all(extreme_df.is_extreme_value_A_1000 .== false)
        @test all(extreme_df.is_extreme_value_B_1000 .== false)
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
        @test :is_extreme_value_A_20 in propertynames(extreme_df)
        @test :is_extreme_value_B_20 in propertynames(extreme_df)
    end

    @testset "boundary conditions" begin
        # Test with single sample data
        df_single = DataFrame()
        df_single[!, :time] = [0.0]
        df_single[!, :triggers] = [0]
        df_single[!, :A] = [100.0]  # Extreme value
        df_single[!, :B] = [1.0]    # Normal value
        layout_single = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        dat_single = eegfun.ContinuousData(df_single, layout_single, 1000, eegfun.AnalysisInfo())

        # Test detect_eog_onsets! with single sample
        eegfun.detect_eog_onsets!(dat_single, 50, :A, :is_eog)
        @test dat_single.data.is_eog[1] == true

        # Test is_extreme_value with single sample (combined mode)
        extreme_single_mask = eegfun.is_extreme_value(dat_single, 50)
        @test extreme_single_mask[1] == true

        # Test is_extreme_value with single sample (separate mode)
        extreme_single = eegfun.is_extreme_value(dat_single, 50, mode = :separate)
        @test extreme_single.is_extreme_value_A_50[1] == true
        @test extreme_single.is_extreme_value_B_50[1] == false

        # Test n_extreme_value with single sample (separate mode)
        count_single = eegfun.n_extreme_value(dat_single, 50, mode = :separate)
        @test count_single.n_extreme[count_single.channel .== :A][1] == 1
        @test count_single.n_extreme[count_single.channel .== :B][1] == 0
    end

    @testset "empty data" begin
        # Test with empty DataFrame
        df_empty = DataFrame()
        df_empty[!, :time] = Float64[]
        df_empty[!, :triggers] = Int[]
        df_empty[!, :A] = Float64[]
        df_empty[!, :B] = Float64[]
        layout_empty = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        dat_empty = eegfun.ContinuousData(df_empty, layout_empty, 1000, eegfun.AnalysisInfo())

        # Test is_extreme_value with empty data (combined mode)
        extreme_empty_mask = eegfun.is_extreme_value(dat_empty, 20)
        @test length(extreme_empty_mask) == 0

        # Test is_extreme_value with empty data (separate mode)
        extreme_empty = eegfun.is_extreme_value(dat_empty, 20, mode = :separate)
        @test size(extreme_empty, 1) == 0
        @test :is_extreme_value_A_20 in propertynames(extreme_empty)
        @test :is_extreme_value_B_20 in propertynames(extreme_empty)

        # Test n_extreme_value with empty data (separate mode)
        count_empty = eegfun.n_extreme_value(dat_empty, 20, mode = :separate)
        @test all(count_empty.n_extreme .== 0)
    end

    @testset "channel overwriting" begin
        dat = create_test_data_with_artifacts()

        # Test that detect_eog_onsets! overwrites existing output channel
        dat.data[!, :existing_output] = fill(false, size(dat.data, 1))
        eegfun.detect_eog_onsets!(dat, 20, :A, :existing_output)
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

    @testset "performance and memory" begin
        # Test with larger dataset
        dat_large = create_test_data_with_artifacts(n = 10000)

        # Test that functions work with larger data (combined mode)
        extreme_large_mask = eegfun.is_extreme_value(dat_large, 20)
        @test length(extreme_large_mask) == 10000

        # Test separate mode with larger data
        extreme_large = eegfun.is_extreme_value(dat_large, 20, mode = :separate)
        @test size(extreme_large, 1) == 10000
        @test :is_extreme_value_A_20 in propertynames(extreme_large)
        @test :is_extreme_value_B_20 in propertynames(extreme_large)

        # Test mutating version with larger data (default combined mode)
        eegfun.is_extreme_value!(dat_large, 20)
        @test any(occursin("is_extreme_value_20", string(col)) for col in names(dat_large.data))

        # Test mutating version with larger data (separate mode)
        dat_large2 = create_test_data_with_artifacts(n = 10000)
        eegfun.is_extreme_value!(dat_large2, 20, mode = :separate)
        @test any(occursin("is_extreme_value_A_20", string(col)) for col in names(dat_large2.data))
        @test any(occursin("is_extreme_value_B_20", string(col)) for col in names(dat_large2.data))
    end

end
