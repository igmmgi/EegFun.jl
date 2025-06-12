using Test
using DataFrames
using eegfun

@testset "Data Utilities" begin

    # Test head and tail
    @testset "head and tail" begin

        df = DataFrame(time = (0:9) ./ 1000, a=1:10, b=11:20)
        layout = DataFrame(label=[:a, :b])
        eeg = eegfun.ContinuousData(df, layout, 1000, eegfun.AnalysisInfo())

        # Test head
        @test eegfun.head(eeg) === nothing  # Default n=5
        @test eegfun.head(eeg, n=3) === nothing

        # Test tail
        @test eegfun.tail(eeg) === nothing  # Default n=5
        @test eegfun.tail(eeg, n=3) === nothing
    end

    # Test datarange
    @testset "datarange" begin
        @test eegfun.datarange([1.0, 2.0, 3.0]) == 2.0
        @test eegfun.datarange([-1.0, 0.0, 1.0]) == 2.0
        @test eegfun.datarange([1.0]) == 0.0
    end

    # Test colmeans
    @testset "colmeans" begin
        # Test DataFrame
        df = DataFrame(a=[1.0, 2.0, 3.0], b=[4.0, 5.0, 6.0])
        @test eegfun.colmeans(df, [:a]) ≈ [1.0, 2.0, 3.0]  # Each row's mean of column a
        @test eegfun.colmeans(df, [:b]) ≈ [4.0, 5.0, 6.0]  # Each row's mean of column b
        @test eegfun.colmeans(df, [:a, :b]) ≈ [2.5, 3.5, 4.5]  # Each row's mean across both columns

        # Test Matrix
        mat = [1.0 2.0; 3.0 4.0; 5.0 6.0]  # 3×2 matrix
        @test eegfun.colmeans(mat) ≈ [1.5, 3.5, 5.5]  # Each row's mean across both columns
        @test eegfun.colmeans(mat, [1]) ≈ [1.0, 3.0, 5.0]  # Each row's mean of first column
        @test eegfun.colmeans(mat, [2]) ≈ [2.0, 4.0, 6.0]  # Each row's mean of second column
    end

    # Test data_limits
    @testset "data_limits" begin
        df = DataFrame(time=[1.0, 2.0, 3.0], value=[4.0, 5.0, 6.0])
        
        # Test data_limits_x
        @test eegfun.data_limits_x(df) == (1.0, 3.0)
        @test eegfun.data_limits_x(df, col=:value) == (4.0, 6.0)

        # Test data_limits_y
        @test eegfun.data_limits_y(df, [:value]) == [4.0, 6.0]
        @test eegfun.data_limits_y(df, [:time, :value]) == [1.0, 6.0]
        @test eegfun.data_limits_y(df, [:time]) == [1.0, 3.0]  # Single column

        # Test empty data
        empty_df = DataFrame(time=Float64[], value=Float64[])
        @test eegfun.data_limits_x(empty_df) === nothing
        @test eegfun.data_limits_y(empty_df, [:value]) === nothing
    end

    # Test to_data_frame
    @testset "to_data_frame" begin
        # Create test epoch data
        epoch1 = DataFrame(time=[1.0, 2.0], value=[3.0, 4.0])
        epoch2 = DataFrame(time=[5.0, 6.0], value=[7.0, 8.0])
        layout = DataFrame(label=[:time, :value])
        epoch_data = eegfun.EpochData([epoch1, epoch2], layout, 1000, eegfun.AnalysisInfo())
        
        # Test single EpochData
        result = eegfun.to_data_frame(epoch_data)
        @test size(result) == (4, 2)
        @test result.time == [1.0, 2.0, 5.0, 6.0]
        @test result.value == [3.0, 4.0, 7.0, 8.0]

        # Test Vector of EpochData
        epoch_data_vec = [epoch_data, epoch_data]
        result = eegfun.to_data_frame(epoch_data_vec)
        @test size(result) == (8, 2)
        @test result.time == [1.0, 2.0, 5.0, 6.0, 1.0, 2.0, 5.0, 6.0]
        @test result.value == [3.0, 4.0, 7.0, 8.0, 3.0, 4.0, 7.0, 8.0]

        # Test empty EpochData
        empty_epoch = eegfun.EpochData(DataFrame[], layout, 1000, eegfun.AnalysisInfo())
        @test size(eegfun.to_data_frame(empty_epoch)) == (0, 0)
        @test size(eegfun.to_data_frame([empty_epoch])) == (0, 0)
    end
end 