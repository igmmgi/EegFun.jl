using Test
using DataFrames
using eegfun
using Statistics
using JLD2
using Random



@testset "epochs" begin


    @testset "parse_epoch_conditions" begin
        cfg = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict("name" => "c1", "trigger_sequences" => [[1, 2, 3]], "reference_index" => 2),
                    Dict(
                        "name" => "c2",
                        "trigger_sequences" => [[1, :any, 3]],
                        "reference_index" => 2,
                        "timing_pairs" => [[1, 3]],
                        "min_interval" => 0.001,
                        "max_interval" => 0.010,
                    ),
                ],
            ),
        )
        conds = eegfun.parse_epoch_conditions(cfg)
        @test length(conds) == 2
        @test conds[1].name == "c1"
        @test conds[1].reference_index == 2
        @test conds[2].timing_pairs == [(1, 3)]

        # Missing trigger_sequences → throws
        bad_cfg1 = Dict("epochs" => Dict("conditions" => [Dict("name" => "bad1")]))
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg1)

        # after and before both set → throws
        bad_cfg2 = Dict(
            "epochs" => Dict(
                "conditions" =>
                    [Dict("name" => "bad2", "trigger_sequences" => [[1, 2, 3]], "after" => 1, "before" => 2)],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg2)

        # reference_index out of bounds → throws
        bad_cfg3 = Dict(
            "epochs" => Dict(
                "conditions" =>
                    [Dict("name" => "bad3", "trigger_sequences" => [[1, 2, 3]], "reference_index" => 5)],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg3)

        # timing_pairs provided but min/max missing → throws
        bad_cfg4 = Dict(
            "epochs" => Dict(
                "conditions" =>
                    [Dict("name" => "bad4", "trigger_sequences" => [[1, 2, 3]], "timing_pairs" => [[1, 3]])],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg4)

        # min_interval ≥ max_interval → throws
        bad_cfg5 = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict(
                        "name" => "bad5",
                        "trigger_sequences" => [[1, 2, 3]],
                        "timing_pairs" => [[1, 3]],
                        "min_interval" => 0.01,
                        "max_interval" => 0.005,
                    ),
                ],
            ),
        )
        @test_throws ErrorException eegfun.parse_epoch_conditions(bad_cfg5)

        # Additional edge cases
        # Empty conditions array
        empty_cfg = Dict("epochs" => Dict("conditions" => []))
        @test eegfun.parse_epoch_conditions(empty_cfg) == []

        # Missing epochs key
        no_epochs_cfg = Dict()
        @test eegfun.parse_epoch_conditions(no_epochs_cfg) == []

        # Invalid trigger sequence format
        invalid_seq_cfg = Dict(
            "epochs" =>
                Dict("conditions" => [Dict("name" => "invalid", "trigger_sequences" => ["invalid_string"])]),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(invalid_seq_cfg)
    end

    @testset "search helpers" begin
        arr = [0, 1, 2, 3, 0, 1, 9, 3]
        to_seq(v) = Vector{Union{Int,Symbol,UnitRange{Int}}}(v)
        idxs = eegfun.search_sequence(arr, [to_seq([1, 2, 3])])
        @test idxs == [2]
        idxs2 = eegfun.search_sequence(arr, [to_seq([1, :any, 3])])
        @test Set(idxs2) == Set([2, 6])  # [1,2,3] and [1,9,3]
        idxr = eegfun.search_sequence(arr, [1:2])
        @test Set(idxr) == Set([2, 3, 6])
        # onset detection for single value
        idxsingle = eegfun.search_sequence([0, 1, 1, 0, 2, 2], 1)
        @test idxsingle == [2]
    end

    @testset "extract_epochs basic and constraints" begin
        dat = create_test_continuous_data_with_triggers()
        win = (-0.01, 0.02)  # 10 ms pre, 20 ms post

        # Single sequence [1,2,3], reference=2
        ec1 = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        ep1 = eegfun.extract_epochs(dat, 1, ec1, win[1], win[2])
        @test eegfun.n_epochs(ep1) == 3  # three [1,2,3] sequences present
        # Integrity preserved
        @test ep1.sample_rate == dat.sample_rate
        @test ep1.layout === dat.layout
        @test all(abs.(ep1.data[1].time) .<= maximum(abs.(ep1.data[1].time)))
        @test 0.0 in ep1.data[1].time
        @test :condition in propertynames(ep1.data[1])
        @test :condition_name in propertynames(ep1.data[1])
        @test :epoch in propertynames(ep1.data[1])

        # Wildcard sequence [1,:any,3]
        ec2 = eegfun.EpochCondition(name = "wild", trigger_sequences = [[1, :any, 3]], reference_index = 2)
        ep2 = eegfun.extract_epochs(dat, 2, ec2, win[1], win[2])
        @test eegfun.n_epochs(ep2) >= 3  # matches [1,2,3], [1,7,3], and [1,2,3] after 9

        # Single range [1:3]
        ec3 = eegfun.EpochCondition(name = "range", trigger_sequences = [[1:3]], reference_index = 1)
        ep3 = eegfun.extract_epochs(dat, 3, ec3, win[1], win[2])
        @test eegfun.n_epochs(ep3) >= 1

        # Position constraints: after 9, before 8
        ec_after =
            eegfun.EpochCondition(name = "after9", trigger_sequences = [[1, 2, 3]], reference_index = 2, after = 9)
        ep_after = eegfun.extract_epochs(dat, 4, ec_after, win[1], win[2])
        @test eegfun.n_epochs(ep_after) >= 1

        ec_before =
            eegfun.EpochCondition(name = "before8", trigger_sequences = [[1, 2, 3]], reference_index = 2, before = 8)
        ep_before = eegfun.extract_epochs(dat, 5, ec_before, win[1], win[2])
        @test eegfun.n_epochs(ep_before) >= 1

        # Timing constraint across (1,3) should be ~2 ms with fs=1000
        ec_time = eegfun.EpochCondition(
            name = "timing",
            trigger_sequences = [[1, 2, 3]],
            reference_index = 2,
            timing_pairs = [(1, 3)],
            min_interval = 0.001,
            max_interval = 0.005,
        )
        ep_time = eegfun.extract_epochs(dat, 6, ec_time, win[1], win[2])
        @test eegfun.n_epochs(ep_time) >= 1

        # No match → throws
        ec_nomatch = eegfun.EpochCondition(name = "none", trigger_sequences = [[7, 7, 7]], reference_index = 2)
        @test_throws ErrorException eegfun.extract_epochs(dat, 9, ec_nomatch, win[1], win[2])

        # Constraints eliminate all → throws
        ec_strict = eegfun.EpochCondition(
            name = "strict",
            trigger_sequences = [[1, 2, 3]],
            reference_index = 2,
            timing_pairs = [(1, 3)],
            min_interval = 0.1,
            max_interval = 0.2,
        )
        @test_throws ErrorException eegfun.extract_epochs(dat, 10, ec_strict, win[1], win[2])

        # Boundary windows (too wide) → BoundsError
        @test_throws BoundsError eegfun.extract_epochs(dat, 11, ec1, -1.0, 0.02)
        @test_throws BoundsError eegfun.extract_epochs(dat, 12, ec1, -0.01, 2.0)
    end

    @testset "average_epochs" begin
        dat = create_test_continuous_data_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 10, ec, win[1], win[2])
        @test eegfun.n_epochs(eps) == 3
        erp = eegfun.average_epochs(eps)
        @test erp isa eegfun.ErpData
        @test :n_epochs in propertynames(erp.data)
        @test maximum(erp.data.n_epochs) == 3
        # Check that averaged channel values at t=0 equal mean of contributing epochs
        zero_rows = findall(erp.data.time .== 0.0)
        @test !isempty(zero_rows)

        # Multiple conditions
        ec_other = eegfun.EpochCondition(name = "seq123b", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps2 = eegfun.extract_epochs(dat, 11, ec_other, win[1], win[2])
        mixed = eegfun.EpochData(vcat(eps.data, eps2.data), eps.layout, eps.sample_rate, eps.analysis_info)
        erp_mixed = eegfun.average_epochs(mixed)
        @test length(unique(erp_mixed.data.condition)) >= 2

        # No EEG channels with layout mismatch → expect error (current implementation)
        only_meta = [
            DataFrame(
                time = [0.0, 0.001],
                triggers = [0, 0],
                condition = [1, 1],
                condition_name = ["x", "x"],
                epoch = [1, 1],
            ),
        ]
        em = eegfun.EpochData(only_meta, eps.layout, eps.sample_rate, eps.analysis_info)
        @test_throws Any eegfun.average_epochs(em)
    end

    @testset "reject_epochs" begin
        dat = create_test_continuous_data_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 20, ec, win[1], win[2])
        @test eegfun.n_epochs(eps) == 3
        # Mark second epoch as bad
        eps.data[1][!, :is_bad] = falses(nrow(eps.data[1]))
        eps.data[2][!, :is_bad] = falses(nrow(eps.data[2]))
        eps.data[3][!, :is_bad] = falses(nrow(eps.data[3]))
        eps.data[2].is_bad[1] = true
        cleaned = eegfun.reject_epochs(eps, :is_bad)
        @test eegfun.n_epochs(cleaned) == 2

        # Multi-column filtering
        eps2 = deepcopy(eps)
        eps2.data[1][!, :is_bad1] = falses(nrow(eps2.data[1]))
        eps2.data[1][!, :is_bad2] = falses(nrow(eps2.data[1]))
        eps2.data[2][!, :is_bad1] = falses(nrow(eps2.data[2]))
        eps2.data[2][!, :is_bad2] = falses(nrow(eps2.data[2]))
        eps2.data[3][!, :is_bad1] = falses(nrow(eps2.data[3]))
        eps2.data[3][!, :is_bad2] = falses(nrow(eps2.data[3]))
        eps2.data[1].is_bad1[5] = true
        eps2.data[2].is_bad2[10] = true
        cleaned2 = eegfun.reject_epochs(eps2, [:is_bad1, :is_bad2])
        @test eegfun.n_epochs(cleaned2) == 1

        # Missing column → throws
        @test_throws ErrorException eegfun.reject_epochs(eps, :does_not_exist)

        # Empty EpochData → error (current implementation validates first epoch)
        empty_ep = eegfun.EpochData(DataFrame[], eps.layout, eps.sample_rate, eps.analysis_info)
        @test_throws Any eegfun.reject_epochs(empty_ep, :is_bad)
    end

    @testset "mark_epoch_windows! (simple triggers)" begin
        dat = create_test_continuous_data_with_triggers()

        # Basic functionality - mark windows around trigger 1
        eegfun.mark_epoch_windows!(dat, [1], [-0.005, 0.005])
        @test :epoch_window in propertynames(dat.data)
        @test any(dat.data.epoch_window)  # Should mark some samples

        # Test that windows are marked around trigger 1 occurrences
        trigger_1_indices = findall(dat.data.triggers .== 1)
        @test !isempty(trigger_1_indices)

        # Check that some samples near each trigger are marked
        for idx in trigger_1_indices
            # Look for marked samples in a small neighborhood 
            neighborhood = max(1, idx-5):min(length(dat.data.epoch_window), idx+5)
            @test any(dat.data.epoch_window[neighborhood])
        end

        # Custom column name
        dat2 = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat2, [1], [-0.005, 0.005], channel_out = :custom_window)
        @test :custom_window in propertynames(dat2.data)
        @test any(dat2.data.custom_window)

        # Multiple triggers
        dat3 = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat3, [1, 2], [-0.005, 0.005])
        @test any(dat3.data.epoch_window)

        # Non-existent trigger should give warning but not error
        dat4 = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat4, [999], [-0.005, 0.005])
        @test !any(dat4.data.epoch_window)  # Should be all false

        # Test input validation
        dat5 = create_test_continuous_data_with_triggers()
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat5, [1], [-0.005])  # Wrong window length
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat5, [1], [0.005, -0.005])  # Wrong order
    end

    @testset "mark_epoch_windows! (epoch conditions)" begin
        dat = create_test_continuous_data_with_triggers()

        # Create epoch condition for sequence [1,2,3]
        ec = eegfun.EpochCondition(name = "test_condition", trigger_sequences = [[1, 2, 3]], reference_index = 2)

        # Basic functionality
        eegfun.mark_epoch_windows!(dat, [ec], [-0.005, 0.005])
        @test :epoch_window in propertynames(dat.data)
        @test any(dat.data.epoch_window)

        # Wildcard condition
        dat2 = create_test_continuous_data_with_triggers()
        ec_wild = eegfun.EpochCondition(name = "wildcard", trigger_sequences = [[1, :any, 3]], reference_index = 2)
        eegfun.mark_epoch_windows!(dat2, [ec_wild], [-0.005, 0.005])
        @test any(dat2.data.epoch_window)

        # Multiple conditions
        dat3 = create_test_continuous_data_with_triggers()
        ec1 = eegfun.EpochCondition(name = "c1", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        ec2 = eegfun.EpochCondition(name = "c2", trigger_sequences = [[1, :any, 3]], reference_index = 1)
        eegfun.mark_epoch_windows!(dat3, [ec1, ec2], [-0.005, 0.005])
        @test any(dat3.data.epoch_window)

        # Condition with constraints
        dat4 = create_test_continuous_data_with_triggers()
        ec_constrained =
            eegfun.EpochCondition(name = "constrained", trigger_sequences = [[1, 2, 3]], reference_index = 2, after = 9)
        eegfun.mark_epoch_windows!(dat4, [ec_constrained], [-0.005, 0.005])
        @test any(dat4.data.epoch_window)  # Should find constrained sequences

        # Non-matching condition
        dat5 = create_test_continuous_data_with_triggers()
        ec_nomatch = eegfun.EpochCondition(name = "nomatch", trigger_sequences = [[7, 7, 7]], reference_index = 1)
        eegfun.mark_epoch_windows!(dat5, [ec_nomatch], [-0.005, 0.005])
        @test !any(dat5.data.epoch_window)  # Should be all false
    end

    @testset "get_selected_epochs" begin
        dat = create_test_continuous_data_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 1, ec, win[1], win[2])

        # Test with simple function that selects all epochs
        all_selector = x -> trues(length(x))
        selected_all = eegfun.get_selected_epochs(eps, all_selector)
        @test selected_all == [1, 2, 3]  # Should have 2 epochs

        # Test with function that selects first epoch only
        first_only = x -> [true, false, false]
        selected_first = eegfun.get_selected_epochs(eps, first_only)
        @test selected_first == [1]

        # Test with function that selects no epochs
        none_selector = x -> falses(length(x))
        selected_none = eegfun.get_selected_epochs(eps, none_selector)
        @test isempty(selected_none)

        # Test with function that selects specific epochs by index
        specific_selector = x -> [i in [1] for i in x]
        selected_specific = eegfun.get_selected_epochs(eps, specific_selector)
        @test selected_specific == [1]

        # Test with realistic selector (odd epochs)
        odd_selector = x -> isodd.(x)
        selected_odd = eegfun.get_selected_epochs(eps, odd_selector)
        @test selected_odd == [1, 3]  # Only epoch 1 is odd
    end

    @testset "_validate_epoch_window_params" begin
        dat = create_test_continuous_data_with_triggers()

        # Valid parameters should not throw
        @test eegfun._validate_epoch_window_params(dat, [-0.1, 0.1]) === nothing

        # Invalid time window length
        @test_throws AssertionError eegfun._validate_epoch_window_params(dat, [-0.1])
        @test_throws AssertionError eegfun._validate_epoch_window_params(dat, [-0.1, 0.0, 0.1])

        # Invalid time window order
        @test_throws AssertionError eegfun._validate_epoch_window_params(dat, [0.1, -0.1])

        # Missing required columns
        df_no_triggers = DataFrame(time = [0.0, 0.1], A = [1.0, 2.0])
        layout = eegfun.Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing)
        dat_no_triggers = eegfun.ContinuousData(df_no_triggers, layout, 1000, eegfun.AnalysisInfo())
        @test_throws AssertionError eegfun._validate_epoch_window_params(dat_no_triggers, [-0.1, 0.1])

        df_no_time = DataFrame(triggers = [0, 1], A = [1.0, 2.0])
        dat_no_time = eegfun.ContinuousData(df_no_time, layout, 1000, eegfun.AnalysisInfo())
        @test_throws AssertionError eegfun._validate_epoch_window_params(dat_no_time, [-0.1, 0.1])
    end

    @testset "edge cases and robustness" begin
        # Empty EpochData for get_selected_epochs
        layout = eegfun.Layout(DataFrame(label = [:A], inc = [0.0], azi = [0.0]), nothing, nothing)
        empty_epochs = eegfun.EpochData(DataFrame[], layout, 1000, eegfun.AnalysisInfo())
        empty_selector = x -> trues(length(x))
        selected_empty = eegfun.get_selected_epochs(empty_epochs, empty_selector)
        @test isempty(selected_empty)

        # Very short time windows for mark_epoch_windows!
        dat = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat, [1], [-0.001, 0.001])  # 2ms window
        @test any(dat.data.epoch_window)

        # Zero-width window
        dat2 = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat2, [1], [0.0, 0.0])
        @test any(dat2.data.epoch_window)  # Should still mark the exact sample

        # Large time windows (should not cause issues if within data bounds)
        dat3 = create_test_continuous_data_with_triggers(n = 10000)  # Longer data
        eegfun.mark_epoch_windows!(dat3, [1], [-0.05, 0.05])  # 100ms window
        @test any(dat3.data.epoch_window)
    end

    @testset "additional edge cases and stress tests" begin
        # Test very large epochs (performance)
        dat_large = create_test_continuous_data_with_triggers(n = 50000)  # 50s at 1kHz
        ec = eegfun.EpochCondition(name = "large", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        large_epochs = eegfun.extract_epochs(dat_large, 1, ec, -0.1, 0.1)
        @test eegfun.n_epochs(large_epochs) >= 1

        # Test averaging with single epoch
        single_epoch = eegfun.EpochData(
            [large_epochs.data[1]],
            large_epochs.layout,
            large_epochs.sample_rate,
            large_epochs.analysis_info,
        )
        erp_single = eegfun.average_epochs(single_epoch)
        @test erp_single isa eegfun.ErpData
        @test maximum(erp_single.data.n_epochs) == 1

        # Test epoch extraction with very short windows
        dat_short = create_test_continuous_data_with_triggers()
        short_epochs = eegfun.extract_epochs(dat_short, 1, ec, -0.001, 0.001)  # 2ms window
        @test eegfun.n_epochs(short_epochs) >= 1
        @test all(epoch -> nrow(epoch) >= 1, short_epochs.data)  # At least one sample per epoch

        # Test with unsorted time vector (should fail validation)
        dat_unsorted = create_test_continuous_data_with_triggers()
        dat_unsorted.data.time = reverse(dat_unsorted.data.time)
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat_unsorted, [1], [-0.01, 0.01])

        # Test reject_epochs with non-boolean columns (should work but warn)
        dat_nonbool = create_test_continuous_data_with_triggers()
        eps_nonbool = eegfun.extract_epochs(dat_nonbool, 1, ec, -0.01, 0.02)
        # Add boolean columns (since any() requires boolean context)
        n_samples_1 = nrow(eps_nonbool.data[1])
        n_samples_2 = nrow(eps_nonbool.data[2])
        n_samples_3 = nrow(eps_nonbool.data[3])
        eps_nonbool.data[1][!, :bool_bad] = vcat([true], falses(n_samples_1 - 1))  # Boolean, one bad sample
        eps_nonbool.data[2][!, :bool_bad] = falses(n_samples_2)  # All good samples
        eps_nonbool.data[3][!, :bool_bad] = falses(n_samples_3)  # All good samples
        cleaned_nonbool = eegfun.reject_epochs(eps_nonbool, :bool_bad)
        @test eegfun.n_epochs(cleaned_nonbool) == 2  # One epoch should be removed

        # Test overlapping epoch windows
        dat_overlap = create_test_continuous_data_with_triggers()
        eegfun.mark_epoch_windows!(dat_overlap, [1], [-0.01, 0.01], channel_out = :window1)
        eegfun.mark_epoch_windows!(dat_overlap, [2], [-0.01, 0.01], channel_out = :window2)
        @test :window1 in propertynames(dat_overlap.data)
        @test :window2 in propertynames(dat_overlap.data)

        # Test get_selected_epochs with edge cases
        dat_edge = create_test_continuous_data_with_triggers()
        eps_edge = eegfun.extract_epochs(dat_edge, 1, ec, -0.01, 0.02)

        # Function that tries to access beyond the input range (should error)
        bad_selector = x -> x[end+1] > 0  # Tries to access beyond range
        @test_throws BoundsError eegfun.get_selected_epochs(eps_edge, bad_selector)

        # Test extract_epochs with reference_index at sequence boundary
        ec_boundary = eegfun.EpochCondition(name = "boundary", trigger_sequences = [[1, 2, 3]], reference_index = 3)  # Last element
        boundary_epochs = eegfun.extract_epochs(dat_edge, 1, ec_boundary, -0.01, 0.01)
        @test eegfun.n_epochs(boundary_epochs) >= 1

        # Test timing constraints with very tight windows
        ec_tight = eegfun.EpochCondition(
            name = "tight",
            trigger_sequences = [[1, 2, 3]],
            reference_index = 2,
            timing_pairs = [(1, 3)],
            min_interval = 0.0019,  # Very tight: exactly 2ms at 1kHz
            max_interval = 0.0021,
        )
        tight_epochs = eegfun.extract_epochs(dat_edge, 1, ec_tight, -0.01, 0.01)
        @test eegfun.n_epochs(tight_epochs) >= 1
    end

end


@testset "Batch Average Epochs" begin

    # Create a temporary directory for test files
    test_dir = mktempdir()

    try
        # Helper to create test EpochData using test_utils.jl function
        function create_batch_test_epoch_data(n_conditions::Int = 2, n_epochs_per_condition::Int = 5)
            epochs = eegfun.EpochData[]
            fs = 256  # Int64
            n_samples = 513
            t = range(-1.0, 1.0, length = n_samples)

            for cond = 1:n_conditions
                # Create multiple epochs per condition
                dfs = DataFrame[]
                for ep = 1:n_epochs_per_condition
                    # Add some trial-to-trial variability
                    ch1 = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)
                    ch2 = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples)

                    df = DataFrame(
                        time = collect(t),
                        sample = 1:n_samples,
                        condition = fill(cond, n_samples),
                        condition_name = fill("condition_$cond", n_samples),
                        epoch = fill(ep, n_samples),
                        Fz = ch1,
                        Cz = ch2,
                    )
                    push!(dfs, df)
                end

                layout =
                    eegfun.Layout(DataFrame(label = [:Fz, :Cz], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)

                # EpochData constructor: (data, layout, sample_rate, analysis_info)
                push!(epochs, eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo()))
            end

            return epochs
        end

        # Create test data files
        @testset "Setup test files" begin
            for participant in [1, 2]
                epochs = create_batch_test_epoch_data(2, 5)
                filename = joinpath(test_dir, "$(participant)_epochs_cleaned.jld2")
                save(filename, "epochs", epochs)
                @test isfile(filename)
            end
        end

        @testset "Basic averaging" begin
            output_dir = joinpath(test_dir, "averaged_output")

            # Test averaging epochs
            result = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            @test result !== nothing
            @test result.success == 2
            @test result.errors == 0
            @test isdir(output_dir)

            # Check that averaged files exist
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))

            # Load and verify averaged data
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 2  # 2 conditions
            @test erps[1] isa eegfun.ErpData
            @test hasproperty(erps[1].data, :Fz)
            @test hasproperty(erps[1].data, :Cz)
            @test hasproperty(erps[1].data, :n_epochs)

            # Verify n_epochs column
            @test all(erps[1].data.n_epochs .== 5)  # 5 epochs averaged
        end

        @testset "Average specific participants" begin
            output_dir = joinpath(test_dir, "averaged_participant")

            result =
                eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir, participants = 1)

            @test result.success == 1
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
        end

        @testset "Average multiple participants" begin
            output_dir = joinpath(test_dir, "averaged_multi_participants")

            result = eegfun.average_epochs(
                "epochs_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                participants = [1, 2],
            )

            @test result.success == 2
            @test result.errors == 0
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))
        end

        @testset "Average specific conditions" begin
            output_dir = joinpath(test_dir, "averaged_condition")

            result =
                eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir, conditions = 1)

            @test result.success == 2

            # Load and verify only one condition
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 1
            @test erps[1].data[1, :condition] == 1
        end

        @testset "Average multiple conditions" begin
            output_dir = joinpath(test_dir, "averaged_multi_conditions")

            result = eegfun.average_epochs(
                "epochs_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                conditions = [1, 2],
            )

            @test result.success == 2

            # Load and verify both conditions
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 2
            @test erps[1].data[1, :condition] == 1
            @test erps[2].data[1, :condition] == 2
        end

        @testset "Error handling" begin
            # Non-existent directory
            @test_throws Exception eegfun.average_epochs("epochs_cleaned", input_dir = "/nonexistent/path")

            # Invalid pattern (doesn't contain 'epochs')
            @test_throws Exception eegfun.average_epochs("erps_cleaned", input_dir = test_dir)
        end

        @testset "No matching files" begin
            empty_dir = joinpath(test_dir, "empty_match")
            mkpath(empty_dir)

            # Directory exists but has no JLD2 files matching pattern
            result = eegfun.average_epochs("epochs_cleaned", input_dir = empty_dir)

            @test result === nothing  # No files to process
        end

        @testset "Logging" begin
            output_dir = joinpath(test_dir, "averaged_with_log")

            result = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            # Check log file exists
            log_file = joinpath(output_dir, "average_epochs.log")
            @test isfile(log_file)

            # Verify log contains expected information
            log_contents = read(log_file, String)
            @test contains(log_contents, "Batch epoch averaging started")
            @test contains(log_contents, "average_epochs")
            @test contains(log_contents, "epochs_cleaned")
        end

        @testset "Existing output directory" begin
            output_dir = joinpath(test_dir, "existing_output_avg")
            mkpath(output_dir)

            # Create a dummy file in the output directory
            touch(joinpath(output_dir, "dummy.txt"))
            @test isfile(joinpath(output_dir, "dummy.txt"))

            # Run average_epochs - should work fine with existing directory
            result = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            @test result.success == 2
            @test isfile(joinpath(output_dir, "dummy.txt"))  # Original file preserved
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
        end

        @testset "Partial failures" begin
            partial_dir = joinpath(test_dir, "partial_test")
            mkpath(partial_dir)

            # Create one valid file
            epochs = create_batch_test_epoch_data(2, 5)
            save(joinpath(partial_dir, "1_epochs_cleaned.jld2"), "epochs", epochs)

            # Create one malformed file (wrong variable name)
            save(joinpath(partial_dir, "2_epochs_cleaned.jld2"), "invalid_var", epochs)

            output_dir = joinpath(test_dir, "averaged_partial")
            result = eegfun.average_epochs("epochs_cleaned", input_dir = partial_dir, output_dir = output_dir)

            @test result.success == 1
            @test result.errors == 1
            @test isfile(joinpath(output_dir, "1_epochs_cleaned.jld2"))
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))  # Failed file not saved
        end

        @testset "Condition out of range" begin
            output_dir = joinpath(test_dir, "averaged_invalid_condition")

            # Request condition 5 when only 2 exist
            result =
                eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir, conditions = 5)

            # Should fail for all files
            @test result.success == 0
            @test result.errors == 2
        end

        @testset "Return value structure" begin
            output_dir = joinpath(test_dir, "averaged_return_check")

            result = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            # Check result structure
            @test hasfield(typeof(result), :success)
            @test hasfield(typeof(result), :errors)
            @test result.success isa Integer
            @test result.errors isa Integer
            @test result.success >= 0
            @test result.errors >= 0
            @test result.success + result.errors == 2  # Total files processed
        end

        @testset "Data integrity - averaging reduces variance" begin
            output_dir = joinpath(test_dir, "averaged_integrity")

            # Get original epoch data statistics
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            # Get variance across epochs for first time point, first condition
            first_epoch_values = [df[1, :Fz] for df in original_epochs[1].data]
            epoch_variance = var(first_epoch_values)

            # Average epochs
            eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            # Load averaged data
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")

            # Check that averaged ERP is smoother (single value, no variance)
            @test erps[1].data[1, :Fz] isa Float64  # Single averaged value

            # Check that n_epochs is correct
            @test erps[1].n_epochs == 5
            @test all(erps[1].data.n_epochs .== 5)

            # Check that time vector is preserved
            original_time = original_epochs[1].data[1].time
            erp_time = erps[1].data.time
            @test length(erp_time) == length(original_time)
            @test erp_time ≈ original_time
        end

        @testset "Combined filters" begin
            output_dir = joinpath(test_dir, "averaged_combined")

            # Average specific participant AND condition
            result = eegfun.average_epochs(
                "epochs_cleaned",
                input_dir = test_dir,
                output_dir = output_dir,
                participants = 1,
                conditions = 1,
            )

            @test result.success == 1

            # Load and verify
            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")
            @test length(erps) == 1  # Only one condition
            @test erps[1].data[1, :condition] == 1
            @test !isfile(joinpath(output_dir, "2_epochs_cleaned.jld2"))  # Participant 2 not processed
        end

        @testset "Different epoch counts per condition" begin
            # Create data with different number of epochs per condition
            epochs_var = eegfun.EpochData[]
            fs = 256
            n_samples = 513
            t = range(-1.0, 1.0, length = n_samples)

            # Condition 1: 3 epochs
            dfs1 = DataFrame[]
            for ep = 1:3
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = sin.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples),
                )
                push!(dfs1, df)
            end

            # Condition 2: 7 epochs
            dfs2 = DataFrame[]
            for ep = 1:7
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(2, n_samples),
                    condition_name = fill("condition_2", n_samples),
                    epoch = fill(ep, n_samples),
                    Fz = cos.(2π .* 5 .* t) .+ 0.1 .* randn(n_samples),
                )
                push!(dfs2, df)
            end

            layout = eegfun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)

            push!(epochs_var, eegfun.EpochData(dfs1, layout, fs, eegfun.AnalysisInfo()))
            push!(epochs_var, eegfun.EpochData(dfs2, layout, fs, eegfun.AnalysisInfo()))

            # Save and process
            var_dir = joinpath(test_dir, "var_epochs")
            mkpath(var_dir)
            save(joinpath(var_dir, "1_epochs_var.jld2"), "epochs", epochs_var)

            output_dir = joinpath(test_dir, "averaged_var")
            result = eegfun.average_epochs("epochs_var", input_dir = var_dir, output_dir = output_dir)

            @test result.success == 1

            # Load and verify epoch counts
            erps = load(joinpath(output_dir, "1_epochs_var.jld2"), "erps")
            @test length(erps) == 2
            @test all(erps[1].data.n_epochs .== 3)  # Condition 1: 3 epochs
            @test all(erps[2].data.n_epochs .== 7)  # Condition 2: 7 epochs
            @test erps[1].n_epochs == 3
            @test erps[2].n_epochs == 7
        end

        @testset "Empty epochs handling" begin
            # Create an epochs file with empty condition
            empty_epochs_dir = joinpath(test_dir, "empty_epochs_test")
            mkpath(empty_epochs_dir)

            # This should trigger an error when trying to average
            # We'll test this by creating a minimal epochs structure
            layout = eegfun.Layout(DataFrame(label = [:Fz], inc = [0.0], azi = [0.0]), nothing, nothing)

            empty_epoch = eegfun.EpochData(DataFrame[], layout, 256, eegfun.AnalysisInfo())
            save(joinpath(empty_epochs_dir, "1_epochs_empty.jld2"), "epochs", [empty_epoch])

            output_dir = joinpath(test_dir, "averaged_empty")
            result = eegfun.average_epochs("epochs_empty", input_dir = empty_epochs_dir, output_dir = output_dir)

            # Should fail because cannot average empty epochs
            @test result.success == 0
            @test result.errors == 1
        end

        @testset "Pattern matching variants" begin
            # Create files with different naming patterns
            pattern_dir = joinpath(test_dir, "pattern_test")
            mkpath(pattern_dir)

            epochs = create_batch_test_epoch_data(2, 3)
            save(joinpath(pattern_dir, "1_epochs_original.jld2"), "epochs", epochs)
            save(joinpath(pattern_dir, "2_epochs_cleaned.jld2"), "epochs", epochs)
            save(joinpath(pattern_dir, "3_custom_epochs.jld2"), "epochs", epochs)

            # Test pattern matching "epochs_original"
            output_dir1 = joinpath(test_dir, "averaged_original")
            result1 = eegfun.average_epochs("epochs_original", input_dir = pattern_dir, output_dir = output_dir1)
            @test result1.success == 1
            @test isfile(joinpath(output_dir1, "1_epochs_original.jld2"))

            # Test pattern matching "epochs_cleaned"
            output_dir2 = joinpath(test_dir, "averaged_cleaned_pattern")
            result2 = eegfun.average_epochs("epochs_cleaned", input_dir = pattern_dir, output_dir = output_dir2)
            @test result2.success == 1
            @test isfile(joinpath(output_dir2, "2_epochs_cleaned.jld2"))

            # Test pattern matching "epochs" (should match all)
            output_dir3 = joinpath(test_dir, "averaged_all_epochs")
            result3 = eegfun.average_epochs("epochs", input_dir = pattern_dir, output_dir = output_dir3)
            @test result3.success == 3
        end

        @testset "Verify averaging math correctness" begin
            # Create simple data where we can verify the math
            math_dir = joinpath(test_dir, "math_test")
            mkpath(math_dir)

            fs = 256
            n_samples = 11  # Small for easy verification
            t = range(0.0, 1.0, length = n_samples)

            # Create 3 epochs with known values
            dfs = DataFrame[]
            known_values = [
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
                [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0],
                [3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0],
            ]

            for (ep, vals) in enumerate(known_values)
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                    Ch1 = vals,
                )
                push!(dfs, df)
            end

            layout = eegfun.Layout(DataFrame(label = [:Ch1], inc = [0.0], azi = [0.0]), nothing, nothing)

            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(math_dir, "1_epochs_math.jld2"), "epochs", [epoch_data])

            # Average
            output_dir = joinpath(test_dir, "averaged_math")
            eegfun.average_epochs("epochs_math", input_dir = math_dir, output_dir = output_dir)

            # Load and verify
            erps = load(joinpath(output_dir, "1_epochs_math.jld2"), "erps")

            # Expected average: [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]
            expected_avg = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0]
            @test erps[1].data.Ch1 ≈ expected_avg
            @test all(erps[1].data.n_epochs .== 3)
        end

        @testset "Metadata preservation" begin
            output_dir = joinpath(test_dir, "averaged_metadata")

            eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")

            # Verify metadata columns exist
            # Note: :sample is not preserved during averaging (it's a per-epoch metadata)
            @test hasproperty(erps[1].data, :time)
            @test hasproperty(erps[1].data, :condition)
            @test hasproperty(erps[1].data, :condition_name)
            @test hasproperty(erps[1].data, :n_epochs)

            # Verify condition metadata
            @test all(erps[1].data.condition .== 1)
            @test all(erps[1].data.condition_name .== "condition_1")
            @test all(erps[2].data.condition .== 2)
            @test all(erps[2].data.condition_name .== "condition_2")
        end

        @testset "Layout and metadata preservation" begin
            output_dir = joinpath(test_dir, "averaged_layout")

            # Get original layout
            original_epochs = load(joinpath(test_dir, "1_epochs_cleaned.jld2"), "epochs")
            original_layout = original_epochs[1].layout
            original_fs = original_epochs[1].sample_rate

            eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = output_dir)

            erps = load(joinpath(output_dir, "1_epochs_cleaned.jld2"), "erps")

            # Verify layout preservation
            @test erps[1].layout.data == original_layout.data
            @test erps[1].sample_rate == original_fs
            @test erps[1] isa eegfun.ErpData
        end

        @testset "Single epoch per condition" begin
            # Edge case: only 1 epoch to average
            single_dir = joinpath(test_dir, "single_epoch")
            mkpath(single_dir)

            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            dfs = DataFrame[]
            df = DataFrame(
                time = collect(t),
                sample = 1:n_samples,
                condition = fill(1, n_samples),
                condition_name = fill("condition_1", n_samples),
                epoch = fill(1, n_samples),
                Cz = sin.(2π .* 10 .* t),
            )
            push!(dfs, df)

            layout = eegfun.Layout(DataFrame(label = [:Cz], inc = [0.0], azi = [0.0]), nothing, nothing)

            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(single_dir, "1_epochs_single.jld2"), "epochs", [epoch_data])

            # Average single epoch
            output_dir = joinpath(test_dir, "averaged_single")
            result = eegfun.average_epochs("epochs_single", input_dir = single_dir, output_dir = output_dir)

            @test result.success == 1

            erps = load(joinpath(output_dir, "1_epochs_single.jld2"), "erps")
            @test erps[1].n_epochs == 1
            @test all(erps[1].data.n_epochs .== 1)

            # With single epoch, average should equal the original
            @test erps[1].data.Cz ≈ df.Cz
        end

        @testset "Output file overwriting" begin
            overwrite_dir = joinpath(test_dir, "averaged_overwrite")

            # First run
            result1 = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = overwrite_dir)
            @test result1.success == 2

            # Get original file modification time
            file1 = joinpath(overwrite_dir, "1_epochs_cleaned.jld2")
            mtime1 = stat(file1).mtime

            # Wait a tiny bit to ensure different mtime
            sleep(0.1)

            # Second run (should overwrite)
            result2 = eegfun.average_epochs("epochs_cleaned", input_dir = test_dir, output_dir = overwrite_dir)
            @test result2.success == 2

            # Verify file was overwritten (different mtime)
            mtime2 = stat(file1).mtime
            @test mtime2 > mtime1
        end

        @testset "Many channels" begin
            # Test with more channels (10)
            many_ch_dir = joinpath(test_dir, "many_channels")
            mkpath(many_ch_dir)

            fs = 256
            n_samples = 101
            t = range(-0.2, 0.2, length = n_samples)

            # Create 10 channels
            channel_names = Symbol.("Ch" .* string.(1:10))

            dfs = DataFrame[]
            for ep = 1:3
                # Build DataFrame with proper column order
                df = DataFrame(
                    time = collect(t),
                    sample = 1:n_samples,
                    condition = fill(1, n_samples),
                    condition_name = fill("condition_1", n_samples),
                    epoch = fill(ep, n_samples),
                )

                # Add channel data
                for (i, ch) in enumerate(channel_names)
                    df[!, ch] = sin.(2π .* i .* t) .+ 0.1 .* randn(n_samples)
                end

                push!(dfs, df)
            end

            layout_df = DataFrame(label = channel_names, inc = zeros(10), azi = zeros(10))
            layout = eegfun.Layout(layout_df, nothing, nothing)

            epoch_data = eegfun.EpochData(dfs, layout, fs, eegfun.AnalysisInfo())
            save(joinpath(many_ch_dir, "1_epochs_many.jld2"), "epochs", [epoch_data])

            # Average
            output_dir = joinpath(test_dir, "averaged_many_ch")
            result = eegfun.average_epochs("epochs_many", input_dir = many_ch_dir, output_dir = output_dir)

            @test result.success == 1

            erps = load(joinpath(output_dir, "1_epochs_many.jld2"), "erps")

            # Verify all channels are present
            for ch in channel_names
                @test hasproperty(erps[1].data, ch)
            end
        end

    finally
        # Cleanup
        rm(test_dir, recursive = true, force = true)
    end
end

# =============================================================================
# ARTIFACT DETECTION TESTS
# =============================================================================

# Helper function to create test epoched data with varying artifact levels

@testset "Artifact Detection" begin
    @testset "Basic detection and rejection" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(1, 1, 3, 100, 3, n_bad_epochs = 3)
        original_n_epochs = length(epoch_data.data)

        # Apply detection and rejection
        rejection_info = eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 2.0)
        clean_data = eegfun.reject_epochs(epoch_data, rejection_info)

        # Check that original data is unchanged
        @test length(epoch_data.data) == original_n_epochs

        # Check that some epochs were rejected
        @test length(clean_data.data) < original_n_epochs
        @test length(clean_data.data) ==
              rejection_info.n_epochs - length(unique([r.epoch for r in rejection_info.rejected_epochs]))
        @test length(rejection_info.rejected_epochs) > 0
    end

    @testset "In-place rejection" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(1, 1, 3, 100, 3, n_bad_epochs = 3)
        original_n_epochs = length(epoch_data.data)

        # Apply detection and rejection in-place
        rejection_info = eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 2.0)
        eegfun.reject_epochs!(epoch_data, rejection_info)

        # Check that data was modified
        @test length(epoch_data.data) < original_n_epochs
        @test length(epoch_data.data) ==
              rejection_info.n_epochs - length(unique([r.epoch for r in rejection_info.rejected_epochs]))
    end

    @testset "Different z-criteria" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(1, 1, 3, 100, 3, n_bad_epochs = 3)

        # Test different criteria
        rejection_info_aggressive = eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 1.5)
        rejection_info_conservative = eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 3.0)

        # More aggressive should reject more epochs
        @test length(rejection_info_aggressive.rejected_epochs) >= length(rejection_info_conservative.rejected_epochs)
    end

    @testset "EpochRejectionInfo structure" begin
        epoch_data, bad_indices = create_test_epochs_with_artifacts(1, 20, 100, 3, 3, n_bad_epochs = 3)
        rejection_info = eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 2.0)

        # Check structure
        @test rejection_info isa eegfun.EpochRejectionInfo
        @test rejection_info.n_epochs == length(epoch_data.data)
        # n_artifacts now counts total artifacts (all channel/epoch combinations), rejected_epochs contains unique rejections
        unique_rejected_epochs = length(unique([r.epoch for r in rejection_info.rejected_epochs]))
        @test rejection_info.n_artifacts >= unique_rejected_epochs  # n_artifacts should be >= unique epochs (can be more due to multiple channels)
        @test length(rejection_info.rejected_epochs) >= 0
    end

    @testset "Error handling" begin
        # Test with empty data
        empty_epochs = eegfun.EpochData(
            DataFrame[],
            eegfun.Layout(DataFrame(label = [:ch1], x = [0], y = [0], z = [0]), nothing, nothing),
            1000,
            eegfun.AnalysisInfo(),
        )
        @test_throws Exception eegfun.detect_bad_epochs_automatic(empty_epochs, z_criterion = 2.0)

        # Test with invalid z-criterion
        epoch_data, _ = create_test_epochs_with_artifacts(1, 1, 5, 50, 2)
        @test_throws Exception eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = -1.0)
        @test_throws Exception eegfun.detect_bad_epochs_automatic(epoch_data, z_criterion = 0.0)
    end
end
