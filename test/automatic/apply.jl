using Test
using DataFrames
using Random
using Statistics
using OrderedCollections
using EegFun

@testset "apply" begin
    Random.seed!(1234)

    @testset "apply_analysis_settings! - ContinuousData" begin
        # Create test data
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        original_data = copy(dat.data)

        # Test 1: Empty settings (should do nothing)
        settings = EegFun.AnalysisSettings()
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.data == original_data

        # Test 2: High-pass filter only
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.1, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 0.0

        # Test 3: Low-pass filter only
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 40.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == 0.0
        @test dat.analysis_info.lp_filter == 40.0

        # Test 4: Both filters
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.1, 40.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 40.0

        # Test 5: Rereference (average)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        original_mean = mean(Matrix(dat.data[:, [:Ch1, :Ch2, :Ch3]]))
        settings = EegFun.AnalysisSettings(0.0, 0.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.reference == :avg
        # After average reference, mean should be approximately zero
        new_mean = mean(Matrix(dat.data[:, [:Ch1, :Ch2, :Ch3]]))
        @test abs(new_mean) < 1e-10

        # Test 6: Channel repair (neighbor interpolation)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 5)
        # Create a layout with neighbors for repair
        layout_df = DataFrame(label = [:Ch1, :Ch2, :Ch3, :Ch4, :Ch5], inc = zeros(5), azi = zeros(5))
        # Create neighbors dict manually for testing
        neighbours_dict = OrderedCollections.OrderedDict(
            :Ch1 => EegFun.Neighbours([:Ch2, :Ch3], [1.0, 1.0], [1.0, 1.0]),
            :Ch2 => EegFun.Neighbours([:Ch1, :Ch3, :Ch4], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
            :Ch3 => EegFun.Neighbours([:Ch1, :Ch2, :Ch4], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
            :Ch4 => EegFun.Neighbours([:Ch2, :Ch3, :Ch5], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]),
            :Ch5 => EegFun.Neighbours([:Ch3, :Ch4], [1.0, 1.0], [1.0, 1.0]),
        )
        layout = EegFun.Layout(layout_df, neighbours_dict, nothing)
        dat.layout = layout

        # Set a channel to NaN to test repair
        dat.data[!, :Ch3] .= NaN
        original_ch3 = copy(dat.data.Ch3)

        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, [:Ch3], :neighbor_interpolation, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        # Ch3 should be repaired (no longer all NaN)
        @test !all(isequal(NaN), dat.data.Ch3)

        # Test 7: Selected regions
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, [(0.1, 0.3), (0.5, 0.7)], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test hasproperty(dat.data, :selected_region)
        @test dat.data.selected_region isa AbstractVector{Bool}
        @test length(dat.data.selected_region) == nrow(dat.data)

        # Test 8: Combined settings
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.1, 40.0, :avg, Symbol[], :none, [(0.1, 0.3)], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 40.0
        @test dat.analysis_info.reference == :avg
        @test hasproperty(dat.data, :selected_region)
    end

    @testset "apply_analysis_settings! - ErpData" begin
        # Create test ERP data
        erp = create_test_erp_data(participant = 1, condition = 1, fs = 1000, n_channels = 3)
        original_data = copy(erp.data)

        # Test with high-pass filter
        settings = EegFun.AnalysisSettings(0.1, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(erp, settings)
        @test erp.analysis_info.hp_filter == 0.1

        # Test with rereference
        erp = create_test_erp_data(participant = 1, condition = 1, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(erp, settings)
        @test erp.analysis_info.reference == :avg

        # Test with selected regions
        erp = create_test_erp_data(participant = 1, condition = 1, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, [(-0.1, 0.1)], Int[])
        EegFun.apply_analysis_settings!(erp, settings)
        @test hasproperty(erp.data, :selected_region)
    end

    @testset "apply_analysis_settings! - EpochData" begin
        # Create test epoch data
        epochs = create_test_epoch_data(n = 500, fs = 1000, n_channels = 3, n_epochs = 3)
        original_data = copy(epochs.data)

        # Test with high-pass filter
        settings = EegFun.AnalysisSettings(0.1, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(epochs, settings)
        @test epochs.analysis_info.hp_filter == 0.1

        # Test with rereference
        epochs = create_test_epoch_data(n = 500, fs = 1000, n_channels = 3, n_epochs = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(epochs, settings)
        @test epochs.analysis_info.reference == :avg

        # Test with selected regions
        epochs = create_test_epoch_data(n = 500, fs = 1000, n_channels = 3, n_epochs = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, [(-0.1, 0.1)], Int[])
        EegFun.apply_analysis_settings!(epochs, settings)
        # Each epoch DataFrame should have selected_region column
        for epoch_df in epochs.data
            @test hasproperty(epoch_df, :selected_region)
        end
    end

    @testset "apply_analysis_settings! - with ICA" begin
        # Create test data
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)

        # Create mock ICA result
        ica_result = EegFun.run_ica(dat, n_components = 2)

        # Test without ICA component removal
        settings = EegFun.AnalysisSettings(0.1, 40.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, ica_result, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 40.0
        @test dat.analysis_info.reference == :avg

        # Test with ICA component removal
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        ica_result = EegFun.run_ica(dat, n_components = 2)
        original_n_channels = EegFun.n_channels(dat)

        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], [1])
        EegFun.apply_analysis_settings!(dat, ica_result, settings)
        # After removing 1 component, should still have same number of channels
        # (ICA removal affects data values, not channel count)
        @test EegFun.n_channels(dat) == original_n_channels

        # Test with all settings including ICA
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        ica_result = EegFun.run_ica(dat, n_components = 2)
        settings = EegFun.AnalysisSettings(0.1, 40.0, :avg, Symbol[], :none, [(0.1, 0.3)], [1])
        EegFun.apply_analysis_settings!(dat, ica_result, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 40.0
        @test dat.analysis_info.reference == :avg
        @test hasproperty(dat.data, :selected_region)
    end

    @testset "apply_analysis_settings! - Observable support" begin
        # Test with Observable (Makie.jl)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.Observable(EegFun.AnalysisSettings(0.1, 40.0, :avg, Symbol[], :none, Tuple{Float64,Float64}[], Int[]))

        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == 0.1
        @test dat.analysis_info.lp_filter == 40.0
        @test dat.analysis_info.reference == :avg

        # Test with Observable and ICA
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        ica_result = EegFun.run_ica(dat, n_components = 2)
        settings = EegFun.Observable(EegFun.AnalysisSettings(0.1, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], [1]))

        EegFun.apply_analysis_settings!(dat, ica_result, settings)
        @test dat.analysis_info.hp_filter == 0.1
    end

    @testset "apply_analysis_settings - non-mutating" begin
        # Test non-mutating version
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        original_hp = dat.analysis_info.hp_filter
        original_lp = dat.analysis_info.lp_filter

        settings = EegFun.AnalysisSettings(0.1, 40.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        dat_new = EegFun.apply_analysis_settings(dat, settings)

        # Original should be unchanged
        @test dat.analysis_info.hp_filter == original_hp
        @test dat.analysis_info.lp_filter == original_lp

        # New data should be modified
        @test dat_new.analysis_info.hp_filter == 0.1
        @test dat_new.analysis_info.lp_filter == 40.0

        # Test with ICA
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        ica_result = EegFun.run_ica(dat, n_components = 2)
        original_hp = dat.analysis_info.hp_filter

        settings = EegFun.AnalysisSettings(0.1, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], [1])
        dat_new = EegFun.apply_analysis_settings(dat, ica_result, settings)

        # Original should be unchanged
        @test dat.analysis_info.hp_filter == original_hp

        # New data should be modified
        @test dat_new.analysis_info.hp_filter == 0.1
    end

    @testset "Edge cases" begin
        # Test with zero filters (should not apply)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        original_hp = dat.analysis_info.hp_filter
        original_lp = dat.analysis_info.lp_filter
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.hp_filter == original_hp
        @test dat.analysis_info.lp_filter == original_lp

        # Test with :none reference (should not apply)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        original_ref = dat.analysis_info.reference
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test dat.analysis_info.reference == original_ref

        # Test with empty repaired channels (should not apply)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        # Should complete without error

        # Test with empty selected regions (should not apply)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, settings)
        @test !hasproperty(dat.data, :selected_region)

        # Test with empty ICA components (should not apply)
        dat = create_test_continuous_data(n = 1000, fs = 1000, n_channels = 3)
        ica_result = EegFun.run_ica(dat, n_components = 2)
        settings = EegFun.AnalysisSettings(0.0, 0.0, :none, Symbol[], :none, Tuple{Float64,Float64}[], Int[])
        EegFun.apply_analysis_settings!(dat, ica_result, settings)
        # Should complete without error
    end
end
