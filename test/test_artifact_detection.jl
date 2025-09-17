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
        eegfun.detect_eog_onsets!(dat, 20.0, :A, :is_eog_onset)
        
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
        
        # Test with default threshold
        extreme_df = eegfun.is_extreme_value(dat, 20.0)
        @test extreme_df isa DataFrame
        @test :sample in propertynames(extreme_df)
        @test :A in propertynames(extreme_df)
        @test :B in propertynames(extreme_df)
        @test size(extreme_df, 1) == size(dat.data, 1)
        
        # Test that it detects artifacts in channel A
        @test any(extreme_df.A[100:110])  # Large artifact
        @test any(extreme_df.A[500:505])  # Negative artifact
        @test any(extreme_df.A[800:802])  # Smaller artifact
        
        # Test that channel B (clean) has fewer extreme values
        @test sum(extreme_df.A) > sum(extreme_df.B)
        
        # Test with channel selection
        extreme_subset = eegfun.is_extreme_value(dat, 20.0, channel_selection = eegfun.channels([:A]))
        @test :A in propertynames(extreme_subset)
        @test :B ∉ propertynames(extreme_subset)
        
        # Test empty channel selection
        @test eegfun.is_extreme_value(dat, 20.0, channel_selection = eegfun.channels(Symbol[])) === nothing
    end

    @testset "is_extreme_value!" begin
        dat = create_test_data_with_artifacts()
        original_columns = names(dat.data)
        
        # Test mutating version
        eegfun.is_extreme_value!(dat, 20.0)
        
        # Check that new columns were added
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 2  # One for each channel
        @test any(occursin("is_extreme_value_A_20", string(col)) for col in new_columns)
        @test any(occursin("is_extreme_value_B_20", string(col)) for col in new_columns)
        
        # Test with channel selection
        dat2 = create_test_data_with_artifacts()
        eegfun.is_extreme_value!(dat2, 20.0, channel_selection = eegfun.channels([:A]))
        
        new_columns2 = setdiff(names(dat2.data), original_columns)
        @test length(new_columns2) == 1  # Only one channel selected
    end

    @testset "n_extreme_value" begin
        dat = create_test_data_with_artifacts()
        
        # Test counting extreme values
        count_df = eegfun.n_extreme_value(dat, 20.0)
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
        @test all(count_df.threshold .== 20.0)
        
        # Test with channel selection
        count_subset = eegfun.n_extreme_value(dat, 20.0, channel_selection = eegfun.channels([:A]))
        @test size(count_subset, 1) == 1
        @test count_subset.channel[1] == :A
        
        # Test empty channel selection
        @test eegfun.n_extreme_value(dat, 20.0, channel_selection = eegfun.channels(Symbol[])) === nothing
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
        extreme_df = eegfun.is_extreme_value(dat, 1000.0)
        @test all(extreme_df.A .== false)
        @test all(extreme_df.B .== false)
        
        # Test with very low threshold (should find many extreme values)
        extreme_df_low = eegfun.is_extreme_value(dat, 0.1)
        @test sum(extreme_df_low.A) > sum(extreme_df.A)
        @test sum(extreme_df_low.B) > sum(extreme_df.B)
        
        # Test with zero threshold (should find all values as extreme)
        extreme_df_zero = eegfun.is_extreme_value(dat, 0.0)
        @test all(extreme_df_zero.A .== true)
        @test all(extreme_df_zero.B .== true)
        
        # Test with negative threshold (should find all values as extreme)
        extreme_df_neg = eegfun.is_extreme_value(dat, -10.0)
        @test all(extreme_df_neg.A .== true)
        @test all(extreme_df_neg.B .== true)
    end

    @testset "error handling" begin
        dat = create_test_data_with_artifacts()
        
        # Test detect_eog_onsets! with non-existent channel (should return nothing and log error)
        result = eegfun.detect_eog_onsets!(dat, 20.0, :NonExistentChannel, :output)
        @test result === nothing
        @test :output ∉ propertynames(dat.data)  # Should not add the output column
        
        # Test is_extreme_value with empty channel selection (should return nothing and log error)
        result = eegfun.is_extreme_value(dat, 20.0, channel_selection = eegfun.channels(Symbol[]))
        @test result === nothing
        
        # Test is_extreme_value! with empty channel selection (should return nothing and log error)
        original_columns = names(dat.data)
        result = eegfun.is_extreme_value!(dat, 20.0, channel_selection = eegfun.channels(Symbol[]))
        @test result === nothing
        @test names(dat.data) == original_columns  # Should not add any columns
        
        # Test n_extreme_value with empty channel selection (should return nothing and log error)
        result = eegfun.n_extreme_value(dat, 20.0, channel_selection = eegfun.channels(Symbol[]))
        @test result === nothing
    end

    @testset "data type handling" begin
        # Test with different data types
        dat = create_test_data_with_artifacts()
        
        # Test with Float32 threshold (should be converted to Float64)
        extreme_df = eegfun.is_extreme_value(dat, Float32(20.0))
        @test extreme_df isa DataFrame
        @test :A in propertynames(extreme_df)
        @test :B in propertynames(extreme_df)
        
        # Test with Int threshold (should be converted to Float64)
        extreme_df_int = eegfun.is_extreme_value(dat, 20)
        @test extreme_df_int isa DataFrame
        @test extreme_df_int.A == extreme_df.A
        @test extreme_df_int.B == extreme_df.B
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
        eegfun.detect_eog_onsets!(dat_single, 50.0, :A, :is_eog)
        @test dat_single.data.is_eog[1] == true
        
        # Test is_extreme_value with single sample
        extreme_single = eegfun.is_extreme_value(dat_single, 50.0)
        @test extreme_single.A[1] == true
        @test extreme_single.B[1] == false
        
        # Test n_extreme_value with single sample
        count_single = eegfun.n_extreme_value(dat_single, 50.0)
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
        
        # Test is_extreme_value with empty data
        extreme_empty = eegfun.is_extreme_value(dat_empty, 20.0)
        @test size(extreme_empty, 1) == 0
        @test :A in propertynames(extreme_empty)
        @test :B in propertynames(extreme_empty)
        
        # Test n_extreme_value with empty data
        count_empty = eegfun.n_extreme_value(dat_empty, 20.0)
        @test all(count_empty.n_extreme .== 0)
    end

    @testset "channel overwriting" begin
        dat = create_test_data_with_artifacts()
        
        # Test that detect_eog_onsets! overwrites existing output channel
        dat.data[!, :existing_output] = fill(false, size(dat.data, 1))
        eegfun.detect_eog_onsets!(dat, 20.0, :A, :existing_output)
        @test any(dat.data.existing_output)  # Should be overwritten with actual results
        
        # Test that is_extreme_value! creates new columns even if similar ones exist
        original_columns = names(dat.data)
        eegfun.is_extreme_value!(dat, 20.0)
        new_columns = setdiff(names(dat.data), original_columns)
        @test length(new_columns) == 2  # Should create new columns
        
        # Test with different threshold to ensure different column names
        eegfun.is_extreme_value!(dat, 30.0)
        final_columns = setdiff(names(dat.data), original_columns)
        @test length(final_columns) == 4  # Should have 4 new columns (2 thresholds × 2 channels)
    end

    @testset "performance and memory" begin
        # Test with larger dataset
        dat_large = create_test_data_with_artifacts(n=10000)
        
        # Test that functions work with larger data
        extreme_large = eegfun.is_extreme_value(dat_large, 20.0)
        @test size(extreme_large, 1) == 10000
        @test :A in propertynames(extreme_large)
        @test :B in propertynames(extreme_large)
        
        # Test mutating version with larger data
        eegfun.is_extreme_value!(dat_large, 20.0)
        @test any(occursin("is_extreme_value_A_20", string(col)) for col in names(dat_large.data))
        @test any(occursin("is_extreme_value_B_20", string(col)) for col in names(dat_large.data))
    end

end
