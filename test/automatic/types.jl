using Test
using DataFrames
using OrderedCollections
using EegFun

@testset "Type Definitions and Display Functions" begin
    @testset "Abstract Types" begin
        # Test that abstract types are properly defined
        @test EegFun.EegData isa Type
        @test EegFun.SingleDataFrameEeg <: EegFun.EegData
        @test EegFun.MultiDataFrameEeg <: EegFun.EegData
    end

    @testset "AnalysisInfo" begin
        # Test default constructor
        info = EegFun.AnalysisInfo()
        @test info.reference == :none
        @test info.hp_filter == 0.0
        @test info.lp_filter == 0.0

        # Test keyword constructor
        info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)
        @test info.reference == :avg
        @test info.hp_filter == 0.1
        @test info.lp_filter == 30.0

        # Test mutability
        info.reference = :mastoid
        @test info.reference == :mastoid
    end

    @testset "Neighbours" begin
        # Test constructor
        neighbours = EegFun.Neighbours([:Fz, :Cz], [10.0, 15.0], [0.5, 0.5])
        @test neighbours.channels == [:Fz, :Cz]
        @test neighbours.distances == [10.0, 15.0]
        @test neighbours.weights == [0.5, 0.5]

        # Test immutability - trying to set a non-existent field throws FieldError
        @test_throws FieldError neighbours.electrodes = [:Pz]
    end

    @testset "Layout" begin
        # Create test data
        df = DataFrame(label = [:Fz, :Cz, :Pz], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0])

        # Test constructor
        layout = EegFun.Layout(df, nothing, nothing)
        @test layout.data == df
        @test layout.neighbours === nothing
        @test layout.criterion === nothing

        # Test mutability
        layout.criterion = 50.0
        @test layout.criterion == 50.0

        # Test with neighbours
        neighbours_dict = OrderedDict{Symbol,EegFun.Neighbours}()
        neighbours_dict[:Fz] = EegFun.Neighbours([:Cz], [10.0], [1.0])
        layout_with_neighbours = EegFun.Layout(df, neighbours_dict, 50.0)
        @test layout_with_neighbours.neighbours == neighbours_dict
        @test layout_with_neighbours.criterion == 50.0
    end

    @testset "Layout copy function" begin
        # Create test layout
        df = DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0])
        neighbours_dict = OrderedDict{Symbol,EegFun.Neighbours}()
        neighbours_dict[:Fz] = EegFun.Neighbours([:Cz], [10.0], [1.0])
        layout = EegFun.Layout(df, neighbours_dict, 50.0)

        # Test copy
        copied = copy(layout)
        @test copied.data == layout.data
        @test copied.neighbours == layout.neighbours
        @test copied.criterion == layout.criterion

        # Test independence
        layout.data[1, :label] = :Pz
        @test copied.data[1, :label] == :Fz  # Copy should be unchanged
    end

    @testset "ContinuousData" begin
        # Create test data
        df = DataFrame(time = [0.1, 0.2, 0.3], Fz = [1.0, 2.0, 3.0], Cz = [4.0, 5.0, 6.0])
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        analysis_info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)

        # Test constructor
        continuous_data = EegFun.ContinuousData("test_data", df, layout, 250, analysis_info)
        @test continuous_data.data == df
        @test continuous_data.layout == layout
        @test continuous_data.sample_rate == 250
        @test continuous_data.analysis_info == analysis_info

        # Test inheritance
        @test continuous_data isa EegFun.SingleDataFrameEeg
        @test continuous_data isa EegFun.EegData
    end

    @testset "ErpData" begin
        # Create test data
        df = DataFrame(time = [0.1, 0.2, 0.3], Fz = [1.0, 2.0, 3.0], Cz = [4.0, 5.0, 6.0])
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        analysis_info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)

        # Test constructor
        erp_data = EegFun.ErpData("test_data", 1, "condition_1", df, layout, 250, analysis_info, 10)
        @test erp_data.data == df
        @test erp_data.layout == layout
        @test erp_data.sample_rate == 250
        @test erp_data.analysis_info == analysis_info
        @test erp_data.n_epochs == 10

        # Test inheritance
        @test erp_data isa EegFun.SingleDataFrameEeg
        @test erp_data isa EegFun.EegData
    end

    @testset "EpochData" begin
        # Create test data
        epoch1 = DataFrame(time = [0.1, 0.2], Fz = [1.0, 2.0], Cz = [3.0, 4.0])
        epoch2 = DataFrame(time = [0.1, 0.2], Fz = [5.0, 6.0], Cz = [7.0, 8.0])
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        analysis_info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)

        # Test constructor
        epoch_data = EegFun.EpochData("test_data", 1, "condition_1", [epoch1, epoch2], layout, 250, analysis_info)
        @test epoch_data.data == [epoch1, epoch2]
        @test epoch_data.layout == layout
        @test epoch_data.sample_rate == 250
        @test epoch_data.analysis_info == analysis_info

        # Test inheritance
        @test epoch_data isa EegFun.MultiDataFrameEeg
        @test epoch_data isa EegFun.EegData
    end

    @testset "IntervalTime" begin
        # Test constructor
        interval = EegFun.IntervalTime(start = 0.1, stop = 0.5)
        @test interval.start == 0.1
        @test interval.stop == 0.5

        # Test immutability
        @test_throws ErrorException interval.start = 0.2
    end

    @testset "EpochCondition" begin
        # Test default constructor
        condition = EegFun.EpochCondition(name = "test_condition", trigger_sequences = [[1, 2, 3]])
        @test condition.name == "test_condition"
        @test condition.trigger_sequences == [[1, 2, 3]]
        @test condition.reference_index == 1
        @test condition.timing_pairs === nothing
        @test condition.min_interval === nothing
        @test condition.max_interval === nothing
        @test condition.after === nothing
        @test condition.before === nothing

        # Test with all parameters
        condition = EegFun.EpochCondition(
            name = "complex_condition",
            trigger_sequences = [[1, :any, 3], [1:5]],
            reference_index = 2,
            timing_pairs = [(1, 2)],
            min_interval = 0.1,
            max_interval = 1.0,
            after = 10,
            before = 20,
        )
        @test condition.name == "complex_condition"
        @test condition.trigger_sequences == [[1, :any, 3], [1:5]]
        @test condition.reference_index == 2
        @test condition.timing_pairs == [(1, 2)]
        @test condition.min_interval == 0.1
        @test condition.max_interval == 1.0
        @test condition.after == 10
        @test condition.before == 20
    end

    @testset "IntervalIndex" begin
        # Test constructor
        interval = EegFun.IntervalIndex(start = 10, stop = 50)
        @test interval.start == 10
        @test interval.stop == 50

        # Test immutability
        @test_throws ErrorException interval.start = 20
    end

    @testset "IcaPrms" begin
        # Test constructor with all parameters
        ica_params = EegFun.IcaPrms(
            0.001,  # l_rate
            1000,   # max_iter
            1e-6,   # w_change
            0.9,    # anneal_deg
            0.8,    # anneal_step
            1e12,   # blowup
            0.5,    # blowup_fac
            1000.0, # max_weight
            0.9,    # restart_factor
            0.0,    # degconst
            1e-6,    # default_stop
        )

        @test ica_params.l_rate == 0.001
        @test ica_params.max_iter == 1000
        @test ica_params.w_change == 1e-6
        @test ica_params.anneal_deg == 0.9
        @test ica_params.anneal_step == 0.8
        @test ica_params.blowup == 1e12
        @test ica_params.blowup_fac == 0.5
        @test ica_params.max_weight == 1000.0
        @test ica_params.restart_factor == 0.9
        @test ica_params.degconst == 0.0
        @test ica_params.default_stop == 1e-6

        # Test mutability
        ica_params.l_rate = 0.002
        @test ica_params.l_rate == 0.002
    end

    @testset "InfoIca" begin
        # Create test data
        unmixing = [1.0 0.0; 0.0 1.0]
        mixing = [1.0 0.0; 0.0 1.0]
        sphere = [1.0 0.0; 0.0 1.0]
        variance = [0.5, 0.3]
        scale = 1.0
        mean = [0.0, 0.0]
        ica_label = [:IC1, :IC2]
        removed_activations = OrderedDict{Int,Matrix{Float64}}()
        removed_activations[1] = [1.0 2.0; 3.0 4.0]
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

        # Test constructor
        is_sub_gaussian = falses(2)  # All super-Gaussian for test (all false)
        ica_info = EegFun.InfoIca(
            "test_file.bdf",
            unmixing,
            mixing,
            sphere,
            variance,
            scale,
            mean,
            ica_label,
            removed_activations,
            layout,
            is_sub_gaussian,
        )

        @test ica_info.unmixing == unmixing
        @test ica_info.mixing == mixing
        @test ica_info.sphere == sphere
        @test ica_info.variance == variance
        @test ica_info.scale == scale
        @test ica_info.mean == mean
        @test ica_info.ica_label == ica_label
        @test ica_info.removed_activations == removed_activations
        @test ica_info.layout == layout
        @test ica_info.is_sub_gaussian == is_sub_gaussian
        @test ica_info.filename == "test_file.bdf"

        # Test immutability
        @test_throws ErrorException ica_info.scale = 2.0
    end

    @testset "Display functions" begin
        # Test Layout display
        df = DataFrame(label = [:Fz, :Cz, :Pz], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0])
        layout = EegFun.Layout(df, nothing, nothing)

        # Test that show doesn't throw errors
        @test_nowarn show(stdout, layout)
        @test_nowarn show(stdout, MIME"text/plain"(), layout)

        # Test AnalysisInfo display
        info = EegFun.AnalysisInfo(reference = :avg, hp_filter = 0.1, lp_filter = 30.0)
        @test_nowarn show(stdout, info)

        # Test EegData display (skip due to dependency issues)
        # continuous_data = EegFun.ContinuousData(df, layout, 250, info)
        # @test_nowarn show(stdout, continuous_data)

        # Test InfoIca display
        unmixing = [1.0 0.0; 0.0 1.0]
        mixing = [1.0 0.0; 0.0 1.0]
        sphere = [1.0 0.0; 0.0 1.0]
        variance = [0.5, 0.3]
        ica_label = [:IC1, :IC2]
        removed_activations = OrderedDict{Int,Matrix{Float64}}()
        ica_info = EegFun.InfoIca(
            "test_file.bdf",
            unmixing,
            mixing,
            sphere,
            variance,
            1.0,
            [0.0, 0.0],
            ica_label,
            removed_activations,
            layout,
            falses(2),
        )
        @test_nowarn show(stdout, ica_info)
        @test_nowarn show(stdout, MIME"text/plain"(), ica_info)

        # Test Neighbours OrderedDict display
        neighbours_dict = OrderedDict{Symbol,EegFun.Neighbours}()
        neighbours_dict[:Fz] = EegFun.Neighbours([:Cz], [10.0], [1.0])
        @test_nowarn show(stdout, neighbours_dict)
        @test_nowarn show(stdout, MIME"text/plain"(), neighbours_dict)
    end

    @testset "Copy functions" begin
        # Test InfoIca copy
        unmixing = [1.0 0.0; 0.0 1.0]
        mixing = [1.0 0.0; 0.0 1.0]
        sphere = [1.0 0.0; 0.0 1.0]
        variance = [0.5, 0.3]
        ica_label = [:IC1, :IC2]
        removed_activations = OrderedDict{Int,Matrix{Float64}}()
        removed_activations[1] = [1.0 2.0; 3.0 4.0]
        layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        ica_info = EegFun.InfoIca(
            "test_file.bdf",
            unmixing,
            mixing,
            sphere,
            variance,
            1.0,
            [0.0, 0.0],
            ica_label,
            removed_activations,
            layout,
            falses(2),
        )

        copied = copy(ica_info)
        @test copied.unmixing == ica_info.unmixing
        @test copied.mixing == ica_info.mixing
        @test copied.sphere == ica_info.sphere
        @test copied.variance == ica_info.variance
        @test copied.scale == ica_info.scale
        @test copied.mean == ica_info.mean
        @test copied.ica_label == ica_info.ica_label
        @test copied.removed_activations == ica_info.removed_activations
        @test copied.layout == ica_info.layout

        # Test independence
        ica_info.unmixing[1, 1] = 2.0
        @test copied.unmixing[1, 1] == 1.0  # Copy should be unchanged
    end

    @testset "Filename functions" begin
        # Test filename function for SingleDataFrameEeg
        df = DataFrame(time = [0.1, 0.2], Fz = [1.0, 2.0])
        layout = EegFun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)
        analysis_info = EegFun.AnalysisInfo()
        continuous_data = EegFun.ContinuousData("test_file.jld2", df, layout, 250, analysis_info)

        @test EegFun.filename(continuous_data) == "test_file.jld2"

        # Test filename function for MultiDataFrameEeg
        epoch1 = DataFrame(time = [0.1], Fz = [1.0])
        epoch2 = DataFrame(time = [0.2], Fz = [2.0])
        epoch_data = EegFun.EpochData("test_file.jld2", 1, "condition_1", [epoch1, epoch2], layout, 250, analysis_info)

        @test EegFun.filename(epoch_data) == "test_file.jld2"
    end
end
