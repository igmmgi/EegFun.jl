"""
Test suite for src/analysis/lrp.jl
"""

using Test
using DataFrames


# Use generic create_test_lrp_erp from test_utils.jl
# create_test_lrp_erp(participant, condition, n_timepoints, signal_scale)

@testset "LRP Calculation" begin

    @testset "Basic LRP calculation with automatic pair detection" begin
        # Create left and right response ERPs
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        # Calculate LRP
        lrp_data = EegFun.lrp(erp_left, erp_right)

        # Verify structure
        @test lrp_data isa EegFun.ErpData
        @test lrp_data.sample_rate == erp_left.sample_rate
        @test nrow(lrp_data.data) == nrow(erp_left.data)

        # Verify LRP channels were created
        lrp_channels = EegFun.channel_labels(lrp_data)
        @test :C3 in lrp_channels
        @test :C4 in lrp_channels
        @test :C1 in lrp_channels
        @test :C2 in lrp_channels
        @test :Fp1 in lrp_channels
        @test :Fp2 in lrp_channels

        # Verify midline channel is not included
        @test :Fz ∉ lrp_channels

        # Verify condition name
        @test lrp_data.condition_name == "lrp"

        # Verify metadata columns are preserved
        @test hasproperty(lrp_data.data, :time)
        @test hasproperty(lrp_data.data, :sample)
        @test hasproperty(lrp_data.data, :participant)

        # Verify n_epochs is minimum of the two inputs
        @test lrp_data.n_epochs == min(erp_left.n_epochs, erp_right.n_epochs)
    end

    @testset "Manual LRP calculation verification" begin
        # Create simple test data with known values
        n_timepoints = 50
        time = collect(range(-0.2, 0.8, length = n_timepoints))

        # Create left response ERP
        df_left = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            condition = fill(1, n_timepoints),
            condition_name = fill("left", n_timepoints),
            participant = fill(1, n_timepoints),
            file = fill("test", n_timepoints),
            C3 = fill(1.0, n_timepoints),
            C4 = fill(2.0, n_timepoints),
        )

        # Create right response ERP
        df_right = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            condition = fill(2, n_timepoints),
            condition_name = fill("right", n_timepoints),
            participant = fill(1, n_timepoints),
            file = fill("test", n_timepoints),
            C3 = fill(3.0, n_timepoints),
            C4 = fill(4.0, n_timepoints),
        )

        # Create layouts and ErpData objects
        layout = EegFun.Layout(DataFrame(label = [:C3, :C4], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

        erp_left = EegFun.ErpData("test_data", 1, "condition_1", df_left, layout, 250, EegFun.AnalysisInfo(), 10)
        erp_right = EegFun.ErpData("test_data", 2, "condition_2", df_right, layout, 250, EegFun.AnalysisInfo(), 10)

        # Calculate LRP
        lrp_data = EegFun.lrp(erp_left, erp_right)

        # Manual calculation:
        # C3_left = 1.0, C4_left = 2.0, C3_right = 3.0, C4_right = 4.0
        # LRP_C3 = 0.5 * ((C3_right - C4_right) + (C4_left - C3_left))
        #        = 0.5 * ((3.0 - 4.0) + (2.0 - 1.0))
        #        = 0.5 * (-1.0 + 1.0)
        #        = 0.0
        expected_lrp_c3 = 0.5 * ((3.0 - 4.0) + (2.0 - 1.0))
        @test all(abs.(lrp_data.data.C3 .- expected_lrp_c3) .< 1e-10)

        # LRP_C4 = 0.5 * ((C4_right - C3_right) + (C3_left - C4_left))
        #        = 0.5 * ((4.0 - 3.0) + (1.0 - 2.0))
        #        = 0.5 * (1.0 - 1.0)
        #        = 0.0
        expected_lrp_c4 = 0.5 * ((4.0 - 3.0) + (1.0 - 2.0))
        @test all(abs.(lrp_data.data.C4 .- expected_lrp_c4) .< 1e-10)
    end

    @testset "LRP with specific channel selection" begin
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        # Calculate LRP for specific left channels only (C3, C1)
        # Function will automatically pair with C4, C2
        lrp_data = EegFun.lrp(erp_left, erp_right, channel_selection = EegFun.channels([:C3, :C1]))

        # Verify only specified channels are in the result
        lrp_channels = EegFun.channel_labels(lrp_data)
        @test :C3 in lrp_channels
        @test :C4 in lrp_channels
        @test :C1 in lrp_channels
        @test :C2 in lrp_channels

        # Verify Fp1/Fp2 are not included (not specified)
        @test :Fp1 ∉ lrp_channels
        @test :Fp2 ∉ lrp_channels
    end

    @testset "LRP layout filtering" begin
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        lrp_data = EegFun.lrp(erp_left, erp_right)

        # Verify layout contains only LRP channels
        layout_labels = lrp_data.layout.data.label
        @test length(layout_labels) == length(EegFun.channel_labels(lrp_data))
        @test all(ch in layout_labels for ch in EegFun.channel_labels(lrp_data))
    end

    @testset "Error handling" begin
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        @testset "Mismatched sample rates" begin
            erp_wrong_rate = create_test_lrp_erp(1, 2, 100, 1.5)
            erp_wrong_rate.sample_rate = 500  # Different sample rate

            @test_throws Exception EegFun.lrp(erp_left, erp_wrong_rate)
        end

        @testset "Mismatched time points" begin
            erp_wrong_length = create_test_lrp_erp(1, 2, 50, 1.5)  # Different length

            @test_throws Exception EegFun.lrp(erp_left, erp_wrong_length)
        end

        @testset "Invalid channel selection" begin
            # Non-existent channels
            @test_throws Exception EegFun.lrp(
                erp_left,
                erp_right,
                channel_selection = EegFun.channels([:NonExistent1, :NonExistent2]),
            )
        end

        @testset "Even channel selection (should error)" begin
            # Selecting even channel should fail (needs odd for left hemisphere)
            @test_throws Exception EegFun.lrp(
                erp_left,
                erp_right,
                channel_selection = EegFun.channels([:C4]),  # Even number not valid
            )
        end

        @testset "No lateral pairs detected" begin
            # Create ERP with only midline channels
            n_timepoints = 100
            time = collect(range(-0.2, 0.8, length = n_timepoints))

            df = DataFrame(
                time = time,
                sample = 1:n_timepoints,
                condition = fill(1, n_timepoints),
                condition_name = fill("test", n_timepoints),
                participant = fill(1, n_timepoints),
                file = fill("test", n_timepoints),
                Fz = randn(n_timepoints),
                Cz = randn(n_timepoints),
            )

            layout = EegFun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
            erp_no_pairs = EegFun.ErpData("test_data", 1, "condition_1", df, layout, 250, EegFun.AnalysisInfo(), 10)

            @test_throws Exception EegFun.lrp(erp_no_pairs, erp_no_pairs)
        end
    end

    @testset "LRP with complex channel naming" begin
        # Test with channels that have more complex names
        n_timepoints = 100
        time = collect(range(-0.2, 0.8, length = n_timepoints))

        df = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            condition = fill(1, n_timepoints),
            condition_name = fill("test", n_timepoints),
            participant = fill(1, n_timepoints),
            file = fill("test", n_timepoints),
            CP3 = randn(n_timepoints),
            CP4 = randn(n_timepoints),
            PO7 = randn(n_timepoints),
            PO8 = randn(n_timepoints),
        )

        layout =
            EegFun.Layout(DataFrame(label = [:CP3, :CP4, :PO7, :PO8], inc = zeros(4), azi = zeros(4)), nothing, nothing)
        erp1 = EegFun.ErpData("test_data", 1, "condition_1", df, layout, 250, EegFun.AnalysisInfo(), 10)
        erp2 = EegFun.ErpData("test_data", 2, "condition_2", copy(df), layout, 250, EegFun.AnalysisInfo(), 10)

        lrp_data = EegFun.lrp(erp1, erp2)

        # Should detect CP3/CP4 but not PO7/PO8 (7 and 8 are not consecutive odd/even)
        lrp_channels = EegFun.channel_labels(lrp_data)
        @test :CP3 in lrp_channels
        @test :CP4 in lrp_channels

        # PO7/PO8 detection: 7 is odd, but 8 is even and is 7+1, so it should be detected
        @test :PO7 in lrp_channels
        @test :PO8 in lrp_channels
    end

    @testset "Analysis info preservation" begin
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        # Modify analysis info
        erp_left.analysis_info.reference = :avg
        erp_left.analysis_info.hp_filter = 0.1
        erp_left.analysis_info.lp_filter = 30.0

        lrp_data = EegFun.lrp(erp_left, erp_right)

        # Verify analysis info is preserved
        @test lrp_data.analysis_info.reference == :avg
        @test lrp_data.analysis_info.hp_filter == 0.1
        @test lrp_data.analysis_info.lp_filter == 30.0
    end

    @testset "Time vector preservation" begin
        erp_left = create_test_lrp_erp(1, 1, 100, 1.0)
        erp_right = create_test_lrp_erp(1, 2, 100, 1.5)

        lrp_data = EegFun.lrp(erp_left, erp_right)

        # Verify time vectors are identical
        @test all(lrp_data.data.time .≈ erp_left.data.time)
        @test all(lrp_data.data.time .≈ erp_right.data.time)
    end

    @testset "LRP symmetry property" begin
        # Create symmetric test data
        n_timepoints = 50
        time = collect(range(-0.2, 0.8, length = n_timepoints))

        # For symmetric data where left and right are mirror images,
        # LRP should show specific relationships

        df_left = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            condition = fill(1, n_timepoints),
            condition_name = fill("left", n_timepoints),
            participant = fill(1, n_timepoints),
            file = fill("test", n_timepoints),
            C3 = fill(1.0, n_timepoints),  # Left hemisphere: low activity
            C4 = fill(5.0, n_timepoints),  # Right hemisphere: high activity (contralateral)
        )

        df_right = DataFrame(
            time = time,
            sample = 1:n_timepoints,
            condition = fill(2, n_timepoints),
            condition_name = fill("right", n_timepoints),
            participant = fill(1, n_timepoints),
            file = fill("test", n_timepoints),
            C3 = fill(5.0, n_timepoints),  # Left hemisphere: high activity (contralateral)
            C4 = fill(1.0, n_timepoints),  # Right hemisphere: low activity
        )

        layout = EegFun.Layout(DataFrame(label = [:C3, :C4], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

        erp_left = EegFun.ErpData("test_data", 1, "condition_1", df_left, layout, 250, EegFun.AnalysisInfo(), 10)
        erp_right = EegFun.ErpData("test_data", 2, "condition_2", df_right, layout, 250, EegFun.AnalysisInfo(), 10)

        lrp_data = EegFun.lrp(erp_left, erp_right)

        # For perfectly symmetric data, LRP at C3 and C4 should have opposite signs
        # C3 should be positive, C4 should be negative (or vice versa)
        # LRP_C3 = 0.5 * ((5.0 - 1.0) + (5.0 - 1.0)) = 4.0
        # LRP_C4 = 0.5 * ((1.0 - 5.0) + (1.0 - 5.0)) = -4.0
        @test all(abs.(lrp_data.data.C3 .+ lrp_data.data.C4) .< 1e-10)  # C3 = -C4

        # C3 should be positive (contralateral > ipsilateral)
        @test all(lrp_data.data.C3 .> 0)
        @test all(abs.(lrp_data.data.C3 .- 4.0) .< 1e-10)

        # C4 should be negative
        @test all(lrp_data.data.C4 .< 0)
        @test all(abs.(lrp_data.data.C4 .+ 4.0) .< 1e-10)
    end

    @testset "Batch LRP calculation for multiple condition pairs" begin
        # Create test data with multiple conditions (16 conditions: odd=left, even=right)
        erps = EegFun.ErpData[]
        for cond = 1:16
            # Alternate between "left" and "right" patterns
            is_left = isodd(cond)
            erp = create_test_lrp_erp(1, cond, 100, is_left ? 1.0 : 1.5)
            push!(erps, erp)
        end

        # Define condition pairs (odd=left, even=right)
        condition_pairs = [(i, i+1) for i = 1:2:15]

        # Calculate LRP for all pairs
        lrp_results = EegFun.lrp(erps, condition_pairs)

        # Verify structure
        @test lrp_results isa Vector{EegFun.ErpData}
        @test length(lrp_results) == 8  # 8 pairs from 16 conditions

        # Verify each LRP result
        for (idx, lrp_data) in enumerate(lrp_results)
            @test lrp_data isa EegFun.ErpData
            @test lrp_data.condition == idx

            # Verify condition name format
            left_cond = 2*idx - 1
            right_cond = 2*idx
            expected_name = "lrp_$(left_cond)_$(right_cond)"
            @test lrp_data.condition_name == expected_name

            # Verify channels are present
            @test :C3 in EegFun.channel_labels(lrp_data)
            @test :C4 in EegFun.channel_labels(lrp_data)
        end
    end

    @testset "Batch LRP with specific pairs" begin
        # Create test data
        erps = [create_test_lrp_erp(1, i, 100, 1.0) for i = 1:8]

        # Calculate LRP for specific pairs only
        lrp_results = EegFun.lrp(erps, [(1, 2), (3, 4), (7, 8)])

        @test length(lrp_results) == 3
        @test lrp_results[1].condition_name == "lrp_1_2"
        @test lrp_results[2].condition_name == "lrp_3_4"
        @test lrp_results[3].condition_name == "lrp_7_8"
    end

    @testset "Batch LRP error handling" begin
        erps = [create_test_lrp_erp(1, i, 100, 1.0) for i = 1:4]

        @testset "Invalid condition index - too high" begin
            @test_throws Exception EegFun.lrp(erps, [(1, 2), (5, 6)])
        end

        @testset "Invalid condition index - zero" begin
            @test_throws Exception EegFun.lrp(erps, [(0, 1)])
        end

        @testset "Invalid condition index - negative" begin
            @test_throws Exception EegFun.lrp(erps, [(-1, 2)])
        end
    end

    @testset "Batch LRP with custom channel selection" begin
        erps = [create_test_lrp_erp(1, i, 100, 1.0) for i = 1:4]

        # Calculate LRP for specific left channels only (C3)
        # Function will automatically pair with C4
        lrp_results = EegFun.lrp(erps, [(1, 2), (3, 4)], channel_selection = EegFun.channels([:C3]))

        @test length(lrp_results) == 2

        # Verify only C3/C4 are in the results
        for lrp_data in lrp_results
            channels = EegFun.channel_labels(lrp_data)
            @test :C3 in channels
            @test :C4 in channels
            @test :C1 ∉ channels  # Should not be included
            @test :C2 ∉ channels
        end
    end
end
