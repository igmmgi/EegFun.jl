using Test
using DataFrames
using OrderedCollections
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
        info = eegfun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)
        @test info.reference == :avg
        @test info.hp_filter == 0.1
        @test info.lp_filter == 30.0

        # Test reference function
        @test eegfun.reference(info) == :avg

        # Test filter_info function
        @test eegfun.filter_info(info) == [0.1, 30.0]

        # Test edge cases
        @test eegfun.filter_info(eegfun.AnalysisInfo()) == [0.0, 0.0]
        @test eegfun.reference(eegfun.AnalysisInfo(reference = :mastoid)) == :mastoid
    end

    @testset "Layout" begin
        # Test basic layout creation
        layout_df = DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0])
        layout = eegfun.Layout(layout_df, nothing, nothing)
        @test layout.data == layout_df
        @test layout.neighbours === nothing
        @test layout.criterion === nothing

        # Test layout with neighbours
        neighbours_dict = OrderedDict{Symbol, eegfun.Neighbours}()
        neighbours_dict[:ch1] = eegfun.Neighbours([:ch2], [1.0], [0.5])
        layout_with_neigh = eegfun.Layout(layout_df, neighbours_dict, 10.0)
        @test layout_with_neigh.neighbours == neighbours_dict
        @test layout_with_neigh.criterion == 10.0

        # Test layout copy
        layout_copy = copy(layout)
        @test layout_copy.data == layout.data
        @test layout_copy.neighbours === layout.neighbours
        @test layout_copy.criterion === layout.criterion
        @test layout_copy !== layout  # Should be a different object

        # Test edge cases
        empty_layout = eegfun.Layout(DataFrame(label = Symbol[], x = Float64[], y = Float64[]), nothing, nothing)
        @test size(empty_layout.data) == (0, 3)  # 0 rows, 3 columns (label, x, y)
    end

    @testset "Neighbours" begin
        # Test basic neighbours creation
        neighbours = eegfun.Neighbours([:Fp1, :Fp2], [1.0, 2.0], [0.5, 0.5])
        @test length(neighbours.electrodes) == 2
        @test neighbours.distances == [1.0, 2.0]
        @test neighbours.weights == [0.5, 0.5]

        # Test edge cases
        empty_neighbours = eegfun.Neighbours(Symbol[], Float64[], Float64[])
        @test isempty(empty_neighbours.electrodes)
        @test isempty(empty_neighbours.distances)
        @test isempty(empty_neighbours.weights)
    end

    @testset "ContinuousData" begin
        # Create test data
        df = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        # Test constructor
        data = eegfun.ContinuousData(df, layout, 1000, info)
        @test data isa eegfun.ContinuousData
        @test data isa eegfun.SingleDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info

        # Test basic functions
        @test eegfun.channel_labels(data) == [:ch1, :ch2]
        @test eegfun.all_data(data).time == df.time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 1
        @test eegfun.duration(data) ≈ 1.0
        @test eegfun.has_channels(data, [:ch1, :ch2])

        # Test edge cases
        empty_df = DataFrame(time = Float64[], ch1 = Float64[], ch2 = Float64[])
        empty_data = eegfun.ContinuousData(empty_df, layout, 1000, info)
        @test eegfun.n_samples(empty_data) == 0
        @test eegfun.duration(empty_data) == 0.0
        @test eegfun.n_channels(empty_data) == 2  # Layout still has 2 channels
    end

    @testset "ErpData" begin
        # Create test data
        df = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        # Test constructor
        data = eegfun.ErpData(df, layout, 1000, info, 10)
        @test data isa eegfun.ErpData
        @test data isa eegfun.SingleDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info
        @test data.n_epochs == 10

        # Test basic functions
        @test eegfun.channel_labels(data) == [:ch1, :ch2]
        @test eegfun.all_data(data).time == df.time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 10
        @test eegfun.duration(data) ≈ 1.0

        # Test edge cases
        empty_df = DataFrame(time = Float64[], ch1 = Float64[], ch2 = Float64[])
        empty_erp = eegfun.ErpData(empty_df, layout, 1000, info, 0)
        @test eegfun.n_samples(empty_erp) == 0
        @test eegfun.duration(empty_erp) == 0.0
        @test eegfun.n_epochs(empty_erp) == 0
    end

    @testset "EpochData" begin
        # Create test data
        epochs = [DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001)) for _ = 1:3]
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        # Test constructor
        data = eegfun.EpochData(epochs, layout, 1000, info)
        @test data isa eegfun.EpochData
        @test data isa eegfun.MultiDataFrameEeg
        @test data.sample_rate == 1000
        @test data.analysis_info == info

        # Test basic functions
        @test eegfun.channel_labels(data) == [:ch1, :ch2]
        @test data.data[1].time == epochs[1].time
        @test eegfun.sample_rate(data) == 1000
        @test eegfun.n_samples(data) == 1001
        @test eegfun.n_channels(data) == 2
        @test eegfun.n_epochs(data) == 3
        @test eegfun.duration(data) ≈ 1.0

        # Test edge cases
        empty_epochs = eegfun.EpochData(DataFrame[], layout, 1000, info)
        @test eegfun.n_epochs(empty_epochs) == 0
        @test isempty(empty_epochs.data)

        # Test with epochs of different lengths
        mixed_epochs = [
            DataFrame(time = 0:0.001:0.5, ch1 = rand(501), ch2 = rand(501)),
            DataFrame(time = 0:0.001:1.0, ch1 = rand(1001), ch2 = rand(1001))
        ]
        mixed_data = eegfun.EpochData(mixed_epochs, layout, 1000, info)
        @test eegfun.n_epochs(mixed_data) == 2
        @test eegfun.n_samples(mixed_data) == 501  # Should return length of first epoch
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

        # Test edge cases
        zero_interval = eegfun.IntervalIdx(5, 5)
        @test zero_interval.interval_start == 5
        @test zero_interval.interval_end == 5

        negative_time = eegfun.IntervalTime(-1.0, 1.0)
        @test negative_time.interval_start == -1.0
        @test negative_time.interval_end == 1.0
    end

    @testset "EpochCondition" begin
        # Test basic EpochCondition
        condition = eegfun.EpochCondition(
            name = "test_condition",
            trigger_sequences = [[1, 2, 3]]
        )
        @test condition.name == "test_condition"
        @test condition.trigger_sequences == [[1, 2, 3]]
        @test condition.reference_index == 1
        @test condition.timing_pairs === nothing
        @test condition.min_interval === nothing
        @test condition.max_interval === nothing
        @test condition.after === nothing
        @test condition.before === nothing

        # Test complex EpochCondition
        complex_condition = eegfun.EpochCondition(
            name = "complex_condition",
            trigger_sequences = [[1, :any, 3], [1:5]],
            reference_index = 2,
            timing_pairs = [(1, 3)],
            min_interval = 0.1,
            max_interval = 1.0,
            after = 5,
            before = 10
        )
        @test complex_condition.name == "complex_condition"
        @test length(complex_condition.trigger_sequences) == 2
        @test complex_condition.reference_index == 2
        @test complex_condition.timing_pairs == [(1, 3)]
        @test complex_condition.min_interval == 0.1
        @test complex_condition.max_interval == 1.0
        @test complex_condition.after == 5
        @test complex_condition.before == 10

        # Test edge cases
        empty_condition = eegfun.EpochCondition(
            name = "empty",
            trigger_sequences = Vector{Vector{Union{Int,Symbol,UnitRange{Int}}}}()
        )
        @test isempty(empty_condition.trigger_sequences)
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
        # Create a simple layout for testing
        test_layout = eegfun.Layout(
            DataFrame(label = [:ch1, :ch2, :ch3], x = [0.0, 1.0, 2.0], y = [0.0, 0.0, 0.0]),
            nothing,
            nothing,
        )
        info = eegfun.InfoIca(
            unmixing,
            mixing,
            sphere,
            variance,
            1.0,
            zeros(3),
            [:ic1, :ic2, :ic3],
            eegfun.OrderedDict{Int,Matrix{Float64}}(),
            test_layout,
        )
        @test size(info.unmixing) == (3, 3)
        @test size(info.mixing) == (3, 3)
        @test length(info.ica_label) == 3
        @test length(info.layout.data.label) == 3

        # Test edge cases
        empty_ica = eegfun.InfoIca(
            Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0),
            Float64[],
            1.0,
            Float64[],
            Symbol[],
            OrderedDict{Int,Matrix{Float64}}(),
            eegfun.Layout(DataFrame(label = Symbol[], x = Float64[], y = Float64[]), nothing, nothing)
        )
        @test size(empty_ica.unmixing) == (0, 0)
        @test isempty(empty_ica.ica_label)
        @test isempty(empty_ica.variance)

        # Test InfoIca copy
        info_copy = copy(info)
        @test info_copy.unmixing == info.unmixing
        @test info_copy.mixing == info.mixing
        @test info_copy !== info  # Should be a different object
    end

    @testset "Common Channel Functions" begin
        # Create two datasets with different channels
        df1 = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        df2 = DataFrame(time = 0:0.001:1, ch2 = rand(1001), ch3 = rand(1001))
        layout1 = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        layout2 = eegfun.Layout(DataFrame(label = [:ch2, :ch3], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        data1 = eegfun.ContinuousData(df1, layout1, 1000, info)
        data2 = eegfun.ContinuousData(df2, layout2, 1000, info)

        # Test common_channels
        @test eegfun.common_channels(data1, data2) == [:ch2]

        # Test edge cases
        empty_layout = eegfun.Layout(DataFrame(label = Symbol[], x = Float64[], y = Float64[]), nothing, nothing)
        empty_data = eegfun.ContinuousData(DataFrame(), empty_layout, 1000, info)
        @test isempty(eegfun.common_channels(empty_data, data1))
        @test isempty(eegfun.common_channels(data1, empty_data))
    end

    @testset "DataFrame Functions" begin
        # Test sample_rate for DataFrame with different sampling rates
        df = DataFrame(time = 0:0.001:1)  # 1000 Hz
        @test eegfun.sample_rate(df) == 1000

        df = DataFrame(time = 0:0.002:1)  # 500 Hz
        @test eegfun.sample_rate(df) == 500

        df = DataFrame(time = 0:0.0005:1)  # 2000 Hz
        @test eegfun.sample_rate(df) == 2000

        # Test edge cases
        empty_df = DataFrame(time = Float64[])
        @test_throws InexactError eegfun.sample_rate(empty_df)  # mean(diff([])) returns NaN, Int(NaN) throws error

        single_point_df = DataFrame(time = [0.0])
        @test_throws InexactError eegfun.sample_rate(single_point_df)  # diff([0.0]) returns empty array, mean([]) returns NaN
    end

    @testset "Channel Functions" begin
        # Create test data with extra columns
        df = DataFrame(
            time = 0:0.001:1,
            ch1 = rand(1001),
            ch2 = rand(1001),
            sample = 1:1001,
            triggers = zeros(1001),
            extra1 = rand(1001),
            extra2 = rand(1001),
        )
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()
        data = eegfun.ContinuousData(df, layout, 1000, info)

        # Test extra_labels
        extra = eegfun.extra_labels(data)
        @test :extra1 in extra
        @test :extra2 in extra
        @test :sample in extra  # sample is after channels, so it's extra
        @test :triggers in extra  # triggers is after channels, so it's extra
        @test !(:time in extra)  # time should be excluded (it's metadata)
        @test !(:ch1 in extra)  # channel should be excluded
        @test !(:ch2 in extra)  # channel should be excluded

        # Test has_channels with empty vector
        @test eegfun.has_channels(data, Symbol[])

        # Test has_channels with non-existent channel
        @test !eegfun.has_channels(data, [:nonexistent])

        # Test common_channels with no overlap
        df2 = DataFrame(time = 0:0.001:1, ch3 = rand(1001), ch4 = rand(1001))
        layout2 = eegfun.Layout(DataFrame(label = [:ch3, :ch4], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        data2 = eegfun.ContinuousData(df2, layout2, 1000, info)
        @test isempty(eegfun.common_channels(data, data2))

        # Test edge cases
        # Create data with only time column (no channels, no extra columns)
        time_only_df = DataFrame(time = Float64[])
        time_only_layout = eegfun.Layout(DataFrame(label = Symbol[], x = Float64[], y = Float64[]), nothing, nothing)
        time_only_data = eegfun.ContinuousData(time_only_df, time_only_layout, 1000, info)
        @test isempty(eegfun.extra_labels(time_only_data))  # Returns empty vector when no extra columns exist
    end

    @testset "MultiDataFrameEeg Functions" begin
        # Create test data with multiple epochs
        epochs = [
            DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001)),
            DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001)),
            DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001)),
        ]
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()
        data = eegfun.EpochData(epochs, layout, 1000, info)

        # Test all_data function for MultiDataFrameEeg
        combined_data = eegfun.all_data(data)
        @test combined_data isa DataFrame
        @test nrow(combined_data) == 3003  # 3 epochs * 1001 samples
        @test :ch1 in propertynames(combined_data)
        @test :ch2 in propertynames(combined_data)
        @test :time in propertynames(combined_data)

        # Test edge cases
        empty_epoch_data = eegfun.EpochData(DataFrame[], layout, 1000, info)
        empty_combined = eegfun.all_data(empty_epoch_data)
        @test empty_combined isa DataFrame
        @test nrow(empty_combined) == 0
    end

    @testset "Display Functions" begin
        # Create test data
        df = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
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

        # Test show for Layout
        io = IOBuffer()
        show(io, layout)
        output = String(take!(io))
        @test occursin("Layout (2 channels)", output)

        # Test show for InfoIca
        unmixing = rand(2, 2)
        mixing = rand(2, 2)
        sphere = rand(2, 2)
        variance = [0.6, 0.4]
        test_layout = eegfun.Layout(
            DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]),
            nothing,
            nothing,
        )
        ica_info = eegfun.InfoIca(
            unmixing,
            mixing,
            sphere,
            variance,
            1.0,
            zeros(2),
            [:ic1, :ic2],
            OrderedDict{Int,Matrix{Float64}}(),
            test_layout,
        )
        io = IOBuffer()
        show(io, ica_info)
        output = String(take!(io))
        @test occursin("InfoIca Result", output)
        @test occursin("Components: 2", output)
        @test occursin("Channels: 2", output)
    end

    @testset "Type Conversion and Validation" begin
        # Test that types are properly constructed and validated
        df = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        # Test that all constructors work without errors
        @test eegfun.ContinuousData(df, layout, 1000, info) isa eegfun.ContinuousData
        @test eegfun.ErpData(df, layout, 1000, info, 5) isa eegfun.ErpData
        @test eegfun.EpochData([df, df], layout, 1000, info) isa eegfun.EpochData

        # Test that abstract types are properly implemented
        cont_data = eegfun.ContinuousData(df, layout, 1000, info)
        erp_data = eegfun.ErpData(df, layout, 1000, info, 5)
        epoch_data = eegfun.EpochData([df, df], layout, 1000, info)

        @test cont_data isa eegfun.EegData
        @test erp_data isa eegfun.EegData
        @test epoch_data isa eegfun.EegData

        @test cont_data isa eegfun.SingleDataFrameEeg
        @test erp_data isa eegfun.SingleDataFrameEeg
        @test !(epoch_data isa eegfun.SingleDataFrameEeg)

        @test !(cont_data isa eegfun.MultiDataFrameEeg)
        @test !(erp_data isa eegfun.MultiDataFrameEeg)
        @test epoch_data isa eegfun.MultiDataFrameEeg
    end

    @testset "Error Handling" begin
        # Test that invalid inputs are handled appropriately
        df = DataFrame(time = 0:0.001:1, ch1 = rand(1001), ch2 = rand(1001))
        layout = eegfun.Layout(DataFrame(label = [:ch1, :ch2], x = [0.0, 1.0], y = [0.0, 0.0]), nothing, nothing)
        info = eegfun.AnalysisInfo()

        # Test with empty DataFrame vector (this should work)
        @test eegfun.EpochData(DataFrame[], layout, 1000, info) isa eegfun.EpochData

        # Note: The constructors don't currently validate negative values
        # These tests would need to be added to the constructors if validation is desired
        # For now, we test that the constructors accept these values without error
        @test eegfun.ContinuousData(df, layout, -1000, info) isa eegfun.ContinuousData
        @test eegfun.ErpData(df, layout, 1000, info, -5) isa eegfun.ErpData
    end
end
