"""
Test suite for src/analysis/gfp.jl and src/plots/plot_gfp.jl
"""

using Test
using DataFrames
using Statistics
using GLMakie

# Use generic create_test_erp_data from test_utils.jl
# create_test_erp_data(participant, condition, n_timepoints, n_channels)

@testset "Global Field Power (GFP)" begin
    @testset "Basic GFP calculation" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate GFP
        gfp_result = eegfun.gfp(erp_data)

        # Check structure
        @test gfp_result isa DataFrame
        @test hasproperty(gfp_result, :time)
        @test hasproperty(gfp_result, :gfp)
        @test hasproperty(gfp_result, :condition)
        @test hasproperty(gfp_result, :condition_name)

        # Check dimensions
        @test nrow(gfp_result) == 100

        # Check that GFP values are non-negative (std deviation can't be negative)
        @test all(gfp_result.gfp .>= 0)

        # Check that GFP is not constant (we created varying patterns)
        @test std(gfp_result.gfp) > 0
    end

    @testset "GFP normalization" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate normalized GFP
        gfp_result = eegfun.gfp(erp_data, normalize = true)

        # Check that values are in 0-100 range
        @test all(0 .<= gfp_result.gfp .<= 100)

        # Check that min is 0 and max is 100 (normalization should give this)
        @test minimum(gfp_result.gfp) ≈ 0.0 atol = 1e-10
        @test maximum(gfp_result.gfp) ≈ 100.0 atol = 1e-10
    end

    @testset "GFP with channel selection" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate GFP with subset of channels
        gfp_result = eegfun.gfp(erp_data, channel_selection = eegfun.channels([:Ch1, :Ch2, :Ch3]))

        @test gfp_result isa DataFrame
        @test nrow(gfp_result) == 100

        # Should produce different result than using all channels
        gfp_all = eegfun.gfp(erp_data)
        @test gfp_result.gfp != gfp_all.gfp
    end

    @testset "GFP for multiple datasets" begin
        erps = [create_test_erp_data(1, 1, 100, 5) for _ = 1:3]

        # Calculate GFP for all
        gfp_results = eegfun.gfp(erps)

        @test gfp_results isa Vector{DataFrame}
        @test length(gfp_results) == 3

        for gfp_result in gfp_results
            @test hasproperty(gfp_result, :gfp)
            @test nrow(gfp_result) == 100
        end
    end

    @testset "GFP error handling" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Should error with no channels selected
        @test_throws Exception eegfun.gfp(erp_data, channel_selection = eegfun.channels([]))
    end

    @testset "GFP calculation correctness" begin
        # Create simple data where GFP can be verified manually
        df = DataFrame()
        df.time = [0.0, 0.1, 0.2]
        df.condition = [1, 1, 1]
        df.condition_name = ["test", "test", "test"]

        # Three channels with known values
        df.Ch1 = [1.0, 2.0, 3.0]
        df.Ch2 = [2.0, 3.0, 4.0]
        df.Ch3 = [3.0, 4.0, 5.0]

        layout = eegfun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
            nothing,
            nothing,
        )
        erp_data = eegfun.ErpData(df, layout, 1000.0, eegfun.AnalysisInfo(), 10)

        # Calculate GFP
        gfp_result = eegfun.gfp(erp_data, normalize = false)

        # Verify GFP values (std of [1,2,3], [2,3,4], [3,4,5])
        expected_gfp = [std([1.0, 2.0, 3.0]), std([2.0, 3.0, 4.0]), std([3.0, 4.0, 5.0])]
        @test all(abs.(gfp_result.gfp .- expected_gfp) .< 1e-10)
    end
end

@testset "Global Dissimilarity" begin
    @testset "Basic dissimilarity calculation" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate dissimilarity
        gd_result = eegfun.global_dissimilarity(erp_data)

        # Check structure
        @test gd_result isa DataFrame
        @test hasproperty(gd_result, :time)
        @test hasproperty(gd_result, :dissimilarity)

        # Check dimensions
        @test nrow(gd_result) == 100

        # Check that dissimilarity values are non-negative
        @test all(gd_result.dissimilarity .>= 0)
    end

    @testset "Dissimilarity normalization" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate normalized dissimilarity
        gd_result = eegfun.global_dissimilarity(erp_data, normalize = true)

        # Check that values are in 0-100 range
        @test all(0 .<= gd_result.dissimilarity .<= 100)
        @test minimum(gd_result.dissimilarity) ≈ 0.0 atol = 1e-10
        @test maximum(gd_result.dissimilarity) ≈ 100.0 atol = 1e-10
    end

    @testset "Dissimilarity for multiple datasets" begin
        erps = [create_test_erp_data(1, 1, 100, 5) for _ = 1:3]

        # Calculate dissimilarity for all
        gd_results = eegfun.global_dissimilarity(erps)

        @test gd_results isa Vector{DataFrame}
        @test length(gd_results) == 3

        for gd_result in gd_results
            @test hasproperty(gd_result, :dissimilarity)
            @test nrow(gd_result) == 100
        end
    end
end

@testset "GFP and Dissimilarity combined" begin
    @testset "Basic combined calculation" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate both
        result = eegfun.gfp_and_dissimilarity(erp_data)

        # Check structure
        @test result isa DataFrame
        @test hasproperty(result, :time)
        @test hasproperty(result, :gfp)
        @test hasproperty(result, :dissimilarity)
        @test nrow(result) == 100
    end

    @testset "Combined calculation with normalization" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate both with normalization
        result = eegfun.gfp_and_dissimilarity(erp_data, normalize = true)

        # Check normalization
        @test all(0 .<= result.gfp .<= 100)
        @test all(0 .<= result.dissimilarity .<= 100)
    end

    @testset "Combined vs separate calculations" begin
        erp_data = create_test_erp_data(1, 1, 100, 5)

        # Calculate combined
        combined = eegfun.gfp_and_dissimilarity(erp_data, normalize = true)

        # Calculate separately
        gfp_only = eegfun.gfp(erp_data, normalize = true)
        gd_only = eegfun.global_dissimilarity(erp_data, normalize = true)

        # Results should match
        @test all(abs.(combined.gfp .- gfp_only.gfp) .< 1e-10)
        @test all(abs.(combined.dissimilarity .- gd_only.dissimilarity) .< 1e-10)
    end

    @testset "Combined for multiple datasets" begin
        erps = [create_test_erp_data(1, 1, 100, 5) for _ = 1:3]

        # Calculate for all
        results = eegfun.gfp_and_dissimilarity(erps)

        @test results isa Vector{DataFrame}
        @test length(results) == 3

        for result in results
            @test hasproperty(result, :gfp)
            @test hasproperty(result, :dissimilarity)
        end
    end
end

@testset "Plot GFP" begin
    @testset "Basic plot creation" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)

        # Create plot (without displaying)
        fig = eegfun.plot_gfp(erp_data, display_plot = false)

        @test fig isa Figure
    end

    @testset "Plot with ERP traces" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)

        fig = eegfun.plot_gfp(erp_data, display_plot = false, show_erp_traces = true)

        @test fig isa Figure
    end

    @testset "Plot with dissimilarity" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)

        fig = eegfun.plot_gfp(erp_data, display_plot = false, show_dissimilarity = true)

        @test fig isa Figure
    end

    @testset "Plot all panels" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)

        fig = eegfun.plot_gfp(erp_data, display_plot = false, show_erp_traces = true, show_dissimilarity = true)

        @test fig isa Figure
    end

    @testset "Plot multiple datasets" begin
        erps = [create_test_erp_data(1, 1, 50, 5) for _ = 1:3]

        fig = eegfun.plot_gfp(erps, display_plot = false)

        @test fig isa Figure
    end

    @testset "Plot from pre-computed GFP" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)
        gfp_result = eegfun.gfp(erp_data, normalize = true)

        fig = eegfun.plot_gfp(gfp_result, display_plot = false)

        @test fig isa Figure
    end

    @testset "Plot from pre-computed GFP with dissimilarity" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)
        result = eegfun.gfp_and_dissimilarity(erp_data, normalize = true)

        fig = eegfun.plot_gfp(result, display_plot = false, show_dissimilarity = true)

        @test fig isa Figure
    end

    @testset "Plot multiple pre-computed GFP" begin
        erps = [create_test_erp_data(1, 1, 50, 5) for _ = 1:3]
        gfp_results = eegfun.gfp.(erps, normalize = true)

        fig = eegfun.plot_gfp(gfp_results, display_plot = false)

        @test fig isa Figure
    end

    @testset "Plot with custom styling" begin
        erp_data = create_test_erp_data(1, 1, 50, 5)

        fig = eegfun.plot_gfp(
            erp_data,
            display_plot = false,
            color = :blue,
            linewidth = 3,
            xlim = (-0.1, 0.5),
            title = "Test GFP",
        )

        @test fig isa Figure
    end

    @testset "Plot error handling" begin
        # DataFrame without required columns
        bad_df = DataFrame(x = [1, 2, 3])

        @test_throws Exception eegfun.plot_gfp(bad_df, display_plot = false)
    end
end

@testset "GFP with realistic data patterns" begin
    @testset "Constant signal (zero GFP)" begin
        # Create data where all channels have the same value
        df = DataFrame()
        df.time = [0.0, 0.1, 0.2]
        df.condition = [1, 1, 1]
        df.condition_name = ["test", "test", "test"]
        df.Ch1 = [1.0, 2.0, 3.0]
        df.Ch2 = [1.0, 2.0, 3.0]
        df.Ch3 = [1.0, 2.0, 3.0]

        layout = eegfun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
            nothing,
            nothing,
        )
        erp_data = eegfun.ErpData(df, layout, 1000.0, eegfun.AnalysisInfo(), 10)

        # Calculate GFP
        gfp_result = eegfun.gfp(erp_data, normalize = false)

        # GFP should be zero (all channels identical)
        @test all(gfp_result.gfp .≈ 0.0)
    end

    @testset "High variance signal (high GFP)" begin
        # Create data with high variance across channels
        df = DataFrame()
        df.time = [0.0, 0.1, 0.2]
        df.condition = [1, 1, 1]
        df.condition_name = ["test", "test", "test"]
        df.Ch1 = [10.0, 10.0, 10.0]
        df.Ch2 = [0.0, 0.0, 0.0]
        df.Ch3 = [-10.0, -10.0, -10.0]

        layout = eegfun.Layout(
            DataFrame(label = [:Ch1, :Ch2, :Ch3], inc = [0.0, 0.0, 0.0], azi = [0.0, 0.0, 0.0]),
            nothing,
            nothing,
        )
        erp_data = eegfun.ErpData(df, layout, 1000.0, eegfun.AnalysisInfo(), 10)

        # Calculate GFP
        gfp_result = eegfun.gfp(erp_data, normalize = false)

        # GFP should be high (large variance)
        @test all(gfp_result.gfp .> 5.0)
    end
end
