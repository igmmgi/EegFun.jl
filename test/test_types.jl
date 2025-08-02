using Test
using DataFrames
using eegfun

@testset "Type System Tests" begin
    @testset "Abstract Types" begin
        @test eegfun.EegData <: Any
        @test eegfun.SingleDataFrameEeg <: eegfun.EegData
        @test eegfun.MultiDataFrameEeg <: eegfun.EegData
    end

    @testset "AnalysisInfo" begin
        # Test default constructor
        info = eegfun.AnalysisInfo()
        @test info.reference == :none
        @test info.hp_filter == 0.0
        @test info.lp_filter == 0.0

        # Test custom constructor
        info = eegfun.AnalysisInfo(reference=:avg, hp_filter=0.1, lp_filter=30.0)
        @test info.reference == :avg
        @test info.hp_filter == 0.1
        @test info.lp_filter == 30.0

        # Test reference function
        @test eegfun.reference(info) == :avg

        # Test filter_info function
        @test eegfun.filter_info(info) == [0.1, 30.0]
    end

    @testset "ContinuousData" begin
        # Create test data
        df = DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001))
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        
        # Test constructor
        data = eegfun.ContinuousData(df, layout, 1000, info)
        @test data isa eegfun.ContinuousData
        @test data isa eegfun.SingleDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info

        # Test basic functions
        @test eegfun.channels(data) == [:ch1, :ch2]
        @test eegfun.times(data) == df.time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 1
        @test eegfun.duration(data) ≈ 1.0
        @test eegfun.has_channels(data, [:ch1, :ch2])
    end

    @testset "ErpData" begin
        # Create test data
        df = DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001))
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        
        # Test constructor
        data = eegfun.ErpData(df, layout, 1000, info, 10)
        @test data isa eegfun.ErpData
        @test data isa eegfun.SingleDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info
        @test data.n_epochs == 10

        # Test basic functions
        @test eegfun.channels(data) == [:ch1, :ch2]
        @test eegfun.times(data) == df.time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 1
        @test eegfun.duration(data) ≈ 1.0
    end

    @testset "EpochData" begin
        # Create test data
        epochs = [DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001)) for _ in 1:3]
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        
        # Test constructor
        data = eegfun.EpochData(epochs, layout, 1000, info)
        @test data isa eegfun.EpochData
        @test data isa eegfun.MultiDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info

        # Test basic functions
        @test eegfun.channels(data) == [:ch1, :ch2]
        @test eegfun.times(data) == epochs[1].time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 3
        @test eegfun.duration(data) ≈ 1.0
    end

    @testset "Interval Types" begin
        # Test IntervalIdx
        idx = eegfun.IntervalIdx(1, 10)
        @test idx.interval_start == 1
        @test idx.interval_end == 10

        # Test IntervalTime
        time = eegfun.IntervalTime(0.0, 1.0)
        @test time.interval_start == 0.0
        @test time.interval_end == 1.0
    end

    @testset "ICA Types" begin
        # Test IcaPrms
        prms = eegfun.IcaPrms(0.1, 100, 0.001, 0.9, 0.98, 1e12, 0.5, 1e12, 0.9, 0.0, 1e-7)
        @test prms.l_rate == 0.1
        @test prms.max_iter == 100
        @test prms.w_change == 0.001

        # Test InfoIca
        unmixing = rand(3, 3)
        mixing = rand(3, 3)
        sphere = rand(3, 3)
        variance = rand(3)
        info = eegfun.InfoIca(unmixing, mixing, sphere, variance, 1.0, zeros(3), [:ic1, :ic2, :ic3], [:ch1, :ch2, :ch3], Dict{Int, Matrix{Float64}}())
        @test size(info.unmixing) == (3, 3)
        @test size(info.mixing) == (3, 3)
        @test length(info.ica_label) == 3
        @test length(info.data_label) == 3
    end

    @testset "Neighbours" begin
        # Test Neighbours struct
        neighbours = eegfun.Neighbours([:Fp1, :Fp2], [1.0, 2.0], [0.5, 0.5])
        @test length(neighbours.electrodes) == 2
        @test neighbours.distances == [1.0, 2.0]
        @test neighbours.weights == [0.5, 0.5]
    end

    @testset "Common Channel Functions" begin
        # Create two datasets with different channels
        df1 = DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001))
        df2 = DataFrame(time=0:0.001:1, ch2=rand(1001), ch3=rand(1001))
        layout1 = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        layout2 = DataFrame(label=[:ch2, :ch3], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        
        data1 = eegfun.ContinuousData(df1, layout1, 1000, info)
        data2 = eegfun.ContinuousData(df2, layout2, 1000, info)

        # Test common_channels
        @test eegfun.common_channels(data1, data2) == [:ch2]
    end

    @testset "DataFrame Functions" begin
        # Test sample_rate for DataFrame with different sampling rates
        df = DataFrame(time=0:0.001:1)  # 1000 Hz
        @test eegfun.sample_rate(df) == 1000

        df = DataFrame(time=0:0.002:1)  # 500 Hz
        @test eegfun.sample_rate(df) == 500

        df = DataFrame(time=0:0.0005:1)  # 2000 Hz
        @test eegfun.sample_rate(df) == 2000
    end

    @testset "Channel Functions" begin
        # Create test data with extra columns
        df = DataFrame(
            time=0:0.001:1,
            ch1=rand(1001),
            ch2=rand(1001),
            sample=1:1001,
            triggers=zeros(1001),
            extra1=rand(1001),
            extra2=rand(1001)
        )
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        data = eegfun.ContinuousData(df, layout, 1000, info)

        # Test extra_channels
        extra = eegfun.extra_channels(data)
        @test :extra1 in extra
        @test :extra2 in extra
        @test !(:time in extra)  # time should be excluded
        @test !(:sample in extra)  # sample should be excluded
        @test !(:triggers in extra)  # triggers should be excluded
        @test !(:ch1 in extra)  # channel should be excluded
        @test !(:ch2 in extra)  # channel should be excluded

        # Test has_channels with empty vector
        @test eegfun.has_channels(data, Symbol[])

        # Test has_channels with non-existent channel
        @test !eegfun.has_channels(data, [:nonexistent])

        # Test common_channels with no overlap
        df2 = DataFrame(time=0:0.001:1, ch3=rand(1001), ch4=rand(1001))
        layout2 = DataFrame(label=[:ch3, :ch4], x=[0.0, 1.0], y=[0.0, 0.0])
        data2 = eegfun.ContinuousData(df2, layout2, 1000, info)
        @test isempty(eegfun.common_channels(data, data2))
    end

    @testset "MultiDataFrameEeg Functions" begin
        # Create test data with multiple epochs
        epochs = [
            DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001)),
            DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001)),
            DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001))
        ]
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        data = eegfun.EpochData(epochs, layout, 1000, info)

        # Test data function for MultiDataFrameEeg
        combined_data = eegfun.data(data)
        @test combined_data isa DataFrame
        @test nrow(combined_data) == 3003  # 3 epochs * 1001 samples
        @test :ch1 in propertynames(combined_data)
        @test :ch2 in propertynames(combined_data)
        @test :time in propertynames(combined_data)
    end

    @testset "Display Functions" begin
        # Create test data
        df = DataFrame(time=0:0.001:1, ch1=rand(1001), ch2=rand(1001))
        layout = DataFrame(label=[:ch1, :ch2], x=[0.0, 1.0], y=[0.0, 0.0])
        info = eegfun.AnalysisInfo()
        data = eegfun.ContinuousData(df, layout, 1000, info)

        # Test show for EegData
        io = IOBuffer()
        show(io, data)
        output = String(take!(io))
        @test occursin("Type: eegfun.ContinuousData", output)
        @test occursin("Size: 1 (epoch) x 1001 (rows) x 2 (columns)", output)
        @test occursin("Labels: ch1, ch2", output)
        @test occursin("Duration: 1.0 S", output)
        @test occursin("Sample Rate: 1000", output)

        # Test show for AnalysisInfo
        io = IOBuffer()
        show(io, info)
        output = String(take!(io))
        @test occursin("Reference: none", output)
        @test occursin("HP Filter: 0.0", output)
        @test occursin("LP Filter: 0.0", output)
    end
end 