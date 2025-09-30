using Test
using DataFrames
using eegfun
using Statistics

@testset "epochs" begin

    # Helper: create simple continuous data with designed trigger patterns
    function create_continuous_with_triggers(; n::Int = 1000, fs::Int = 1000)
        t = collect(0:(n-1)) ./ fs
        triggers = zeros(Int, n)
        # Sequence [1,2,3] at consecutive samples starting at idx1
        idx1 = 101
        triggers[idx1] = 1; triggers[idx1+1] = 2; triggers[idx1+2] = 3
        # Wildcard sequence [1, :any, 3] at idx2
        idx2 = 301
        triggers[idx2] = 1; triggers[idx2+1] = 7; triggers[idx2+2] = 3
        # Single-range [1:3] at idx3 (value 2)
        idx3 = 501
        triggers[idx3] = 2
        # Add an 'after' marker (9) before first sequence and a 'before' marker (8) after
        triggers[idx1-10] = 9
        triggers[idx1+5] = 8

        # Second [1,2,3] to test averaging and removal
        idx1b = 701
        triggers[idx1b] = 1; triggers[idx1b+1] = 2; triggers[idx1b+2] = 3

        # Two channels
        A = sin.(2π .* 5 .* t)
        B = cos.(2π .* 7 .* t)
        df = DataFrame(time = t, triggers = triggers, A = A, B = B)
        layout = eegfun.Layout(DataFrame(label = [:A, :B], inc = [0.0, 0.0], azi = [0.0, 0.0]), nothing, nothing)
        dat = eegfun.ContinuousData(copy(df, copycols=true), layout, fs, eegfun.AnalysisInfo())
        return dat
    end

    @testset "parse_epoch_conditions" begin
        cfg = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict(
                        "name" => "c1",
                        "trigger_sequences" => [[1, 2, 3]],
                        "reference_index" => 2,
                    ),
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
                "conditions" => [
                    Dict(
                        "name" => "bad2",
                        "trigger_sequences" => [[1, 2, 3]],
                        "after" => 1,
                        "before" => 2,
                    ),
                ],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg2)

        # reference_index out of bounds → throws
        bad_cfg3 = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict(
                        "name" => "bad3",
                        "trigger_sequences" => [[1, 2, 3]],
                        "reference_index" => 5,
                    ),
                ],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(bad_cfg3)

        # timing_pairs provided but min/max missing → throws
        bad_cfg4 = Dict(
            "epochs" => Dict(
                "conditions" => [
                    Dict(
                        "name" => "bad4",
                        "trigger_sequences" => [[1, 2, 3]],
                        "timing_pairs" => [[1, 3]],
                    ),
                ],
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
            "epochs" => Dict(
                "conditions" => [
                    Dict(
                        "name" => "invalid",
                        "trigger_sequences" => ["invalid_string"],
                    ),
                ],
            ),
        )
        @test_throws Exception eegfun.parse_epoch_conditions(invalid_seq_cfg)
    end

    @testset "search helpers" begin
        arr = [0, 1, 2, 3, 0, 1, 9, 3]
        to_seq(v) = Vector{Union{Int,Symbol,UnitRange{Int}}}(v)
        idxs = eegfun.search_sequences(arr, [to_seq([1, 2, 3])])
        @test idxs == [2]
        idxs2 = eegfun.search_sequences(arr, [to_seq([1, :any, 3])])
        @test Set(idxs2) == Set([2, 6])  # [1,2,3] and [1,9,3]
        idxr = eegfun.search_trigger_ranges(arr, [1:2])
        @test Set(idxr) == Set([2, 3, 6])
        # onset detection for single value
        idxsingle = eegfun.search_sequence([0, 1, 1, 0, 2, 2], 1)
        @test idxsingle == [2]
    end

    @testset "extract_epochs basic and constraints" begin
        dat = create_continuous_with_triggers()
        win = (-0.01, 0.02)  # 10 ms pre, 20 ms post

        # Single sequence [1,2,3], reference=2
        ec1 = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        ep1 = eegfun.extract_epochs(dat, 1, ec1, win[1], win[2])
        @test eegfun.n_epochs(ep1) == 2  # two [1,2,3] sequences present
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
        @test eegfun.n_epochs(ep2) >= 2  # matches both [1,2,3] and [1,7,3]

        # Single range [1:3]
        ec3 = eegfun.EpochCondition(name = "range", trigger_sequences = [[1:3]], reference_index = 1)
        ep3 = eegfun.extract_epochs(dat, 3, ec3, win[1], win[2])
        @test eegfun.n_epochs(ep3) >= 1

        # Position constraints: after 9, before 8
        ec_after = eegfun.EpochCondition(name = "after9", trigger_sequences = [[1, 2, 3]], reference_index = 2, after = 9)
        ep_after = eegfun.extract_epochs(dat, 4, ec_after, win[1], win[2])
        @test eegfun.n_epochs(ep_after) >= 1

        ec_before = eegfun.EpochCondition(name = "before8", trigger_sequences = [[1, 2, 3]], reference_index = 2, before = 8)
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
        dat = create_continuous_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 10, ec, win[1], win[2])
        @test eegfun.n_epochs(eps) == 2
        erp = eegfun.average_epochs(eps)
        @test erp isa eegfun.ErpData
        @test :n_epochs in propertynames(erp.data)
        @test maximum(erp.data.n_epochs) == 2
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
        only_meta = [DataFrame(time = [0.0, 0.001], triggers = [0, 0], condition = [1, 1], condition_name = ["x", "x"], epoch = [1, 1])]
        em = eegfun.EpochData(only_meta, eps.layout, eps.sample_rate, eps.analysis_info)
        @test_throws Any eegfun.average_epochs(em)
    end

    @testset "remove_bad_epochs" begin
        dat = create_continuous_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 20, ec, win[1], win[2])
        @test eegfun.n_epochs(eps) == 2
        # Mark second epoch as bad
        eps.data[1][!, :is_bad] = falses(nrow(eps.data[1]))
        eps.data[2][!, :is_bad] = falses(nrow(eps.data[2]))
        eps.data[2].is_bad[1] = true
        cleaned = eegfun.remove_bad_epochs(eps, :is_bad)
        @test eegfun.n_epochs(cleaned) == 1

        # Multi-column filtering
        eps2 = deepcopy(eps)
        eps2.data[1][!, :is_bad1] = falses(nrow(eps2.data[1]))
        eps2.data[1][!, :is_bad2] = falses(nrow(eps2.data[1]))
        eps2.data[2][!, :is_bad1] = falses(nrow(eps2.data[2]))
        eps2.data[2][!, :is_bad2] = falses(nrow(eps2.data[2]))
        eps2.data[1].is_bad1[5] = true
        eps2.data[2].is_bad2[10] = true
        cleaned2 = eegfun.remove_bad_epochs(eps2, [:is_bad1, :is_bad2])
        @test eegfun.n_epochs(cleaned2) == 0

        # Missing column → throws
        @test_throws ErrorException eegfun.remove_bad_epochs(eps, :does_not_exist)

        # Empty EpochData → error (current implementation validates first epoch)
        empty_ep = eegfun.EpochData(DataFrame[], eps.layout, eps.sample_rate, eps.analysis_info)
        @test_throws Any eegfun.remove_bad_epochs(empty_ep, :is_bad)
    end

    @testset "mark_epoch_windows! (simple triggers)" begin
        dat = create_continuous_with_triggers()
        
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
        dat2 = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat2, [1], [-0.005, 0.005], channel_out = :custom_window)
        @test :custom_window in propertynames(dat2.data)
        @test any(dat2.data.custom_window)
        
        # Multiple triggers
        dat3 = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat3, [1, 2], [-0.005, 0.005])
        @test any(dat3.data.epoch_window)
        
        # Non-existent trigger should give warning but not error
        dat4 = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat4, [999], [-0.005, 0.005])
        @test !any(dat4.data.epoch_window)  # Should be all false
        
        # Test input validation
        dat5 = create_continuous_with_triggers()
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat5, [1], [-0.005])  # Wrong window length
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat5, [1], [0.005, -0.005])  # Wrong order
    end

    @testset "mark_epoch_windows! (epoch conditions)" begin
        dat = create_continuous_with_triggers()
        
        # Create epoch condition for sequence [1,2,3]
        ec = eegfun.EpochCondition(
            name = "test_condition",
            trigger_sequences = [[1, 2, 3]], 
            reference_index = 2
        )
        
        # Basic functionality
        eegfun.mark_epoch_windows!(dat, [ec], [-0.005, 0.005])
        @test :epoch_window in propertynames(dat.data)
        @test any(dat.data.epoch_window)
        
        # Wildcard condition
        dat2 = create_continuous_with_triggers()
        ec_wild = eegfun.EpochCondition(
            name = "wildcard",
            trigger_sequences = [[1, :any, 3]], 
            reference_index = 2
        )
        eegfun.mark_epoch_windows!(dat2, [ec_wild], [-0.005, 0.005])
        @test any(dat2.data.epoch_window)
        
        # Multiple conditions
        dat3 = create_continuous_with_triggers()
        ec1 = eegfun.EpochCondition(name = "c1", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        ec2 = eegfun.EpochCondition(name = "c2", trigger_sequences = [[1, :any, 3]], reference_index = 1)
        eegfun.mark_epoch_windows!(dat3, [ec1, ec2], [-0.005, 0.005])
        @test any(dat3.data.epoch_window)
        
        # Condition with constraints
        dat4 = create_continuous_with_triggers()
        ec_constrained = eegfun.EpochCondition(
            name = "constrained",
            trigger_sequences = [[1, 2, 3]], 
            reference_index = 2,
            after = 9
        )
        eegfun.mark_epoch_windows!(dat4, [ec_constrained], [-0.005, 0.005])
        @test any(dat4.data.epoch_window)  # Should find constrained sequences
        
        # Non-matching condition
        dat5 = create_continuous_with_triggers()
        ec_nomatch = eegfun.EpochCondition(
            name = "nomatch",
            trigger_sequences = [[7, 7, 7]], 
            reference_index = 1
        )
        eegfun.mark_epoch_windows!(dat5, [ec_nomatch], [-0.005, 0.005])
        @test !any(dat5.data.epoch_window)  # Should be all false
    end

    @testset "get_selected_epochs" begin
        dat = create_continuous_with_triggers()
        win = (-0.01, 0.02)
        ec = eegfun.EpochCondition(name = "seq123", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        eps = eegfun.extract_epochs(dat, 1, ec, win[1], win[2])
        
        # Test with simple function that selects all epochs
        all_selector = x -> trues(length(x))
        selected_all = eegfun.get_selected_epochs(eps, all_selector)
        @test selected_all == [1, 2]  # Should have 2 epochs
        
        # Test with function that selects first epoch only
        first_only = x -> [true, false]
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
        @test selected_odd == [1]  # Only epoch 1 is odd
    end

    @testset "_validate_epoch_window_params" begin
        dat = create_continuous_with_triggers()
        
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
        dat = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat, [1], [-0.001, 0.001])  # 2ms window
        @test any(dat.data.epoch_window)
        
        # Zero-width window
        dat2 = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat2, [1], [0.0, 0.0])
        @test any(dat2.data.epoch_window)  # Should still mark the exact sample
        
        # Large time windows (should not cause issues if within data bounds)
        dat3 = create_continuous_with_triggers(n = 10000)  # Longer data
        eegfun.mark_epoch_windows!(dat3, [1], [-0.05, 0.05])  # 100ms window
        @test any(dat3.data.epoch_window)
    end

    @testset "additional edge cases and stress tests" begin
        # Test very large epochs (performance)
        dat_large = create_continuous_with_triggers(n = 50000)  # 50s at 1kHz
        ec = eegfun.EpochCondition(name = "large", trigger_sequences = [[1, 2, 3]], reference_index = 2)
        large_epochs = eegfun.extract_epochs(dat_large, 1, ec, -0.1, 0.1)
        @test eegfun.n_epochs(large_epochs) >= 1
        
        # Test averaging with single epoch
        single_epoch = eegfun.EpochData([large_epochs.data[1]], large_epochs.layout, large_epochs.sample_rate, large_epochs.analysis_info)
        erp_single = eegfun.average_epochs(single_epoch)
        @test erp_single isa eegfun.ErpData
        @test maximum(erp_single.data.n_epochs) == 1
        
        # Test epoch extraction with very short windows
        dat_short = create_continuous_with_triggers()
        short_epochs = eegfun.extract_epochs(dat_short, 1, ec, -0.001, 0.001)  # 2ms window
        @test eegfun.n_epochs(short_epochs) >= 1
        @test all(epoch -> nrow(epoch) >= 1, short_epochs.data)  # At least one sample per epoch
        
        # Test with unsorted time vector (should fail validation)
        dat_unsorted = create_continuous_with_triggers()
        dat_unsorted.data.time = reverse(dat_unsorted.data.time)
        @test_throws AssertionError eegfun.mark_epoch_windows!(dat_unsorted, [1], [-0.01, 0.01])
        
        # Test remove_bad_epochs with non-boolean columns (should work but warn)
        dat_nonbool = create_continuous_with_triggers()
        eps_nonbool = eegfun.extract_epochs(dat_nonbool, 1, ec, -0.01, 0.02)
        # Add boolean columns (since any() requires boolean context)
        n_samples_1 = nrow(eps_nonbool.data[1])
        n_samples_2 = nrow(eps_nonbool.data[2])
        eps_nonbool.data[1][!, :bool_bad] = vcat([true], falses(n_samples_1 - 1))  # Boolean, one bad sample
        eps_nonbool.data[2][!, :bool_bad] = falses(n_samples_2)  # All good samples
        cleaned_nonbool = eegfun.remove_bad_epochs(eps_nonbool, :bool_bad)
        @test eegfun.n_epochs(cleaned_nonbool) == 1  # One epoch should be removed
        
        # Test overlapping epoch windows
        dat_overlap = create_continuous_with_triggers()
        eegfun.mark_epoch_windows!(dat_overlap, [1], [-0.01, 0.01], channel_out = :window1)
        eegfun.mark_epoch_windows!(dat_overlap, [2], [-0.01, 0.01], channel_out = :window2)
        @test :window1 in propertynames(dat_overlap.data)
        @test :window2 in propertynames(dat_overlap.data)
        
        # Test get_selected_epochs with edge cases
        dat_edge = create_continuous_with_triggers()
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
            max_interval = 0.0021
        )
        tight_epochs = eegfun.extract_epochs(dat_edge, 1, ec_tight, -0.01, 0.01)
        @test eegfun.n_epochs(tight_epochs) >= 1
    end

end


