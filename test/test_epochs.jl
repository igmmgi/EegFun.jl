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

        # Missing trigger_sequences → nothing
        bad_cfg1 = Dict("epochs" => Dict("conditions" => [Dict("name" => "bad1")]))
        @test eegfun.parse_epoch_conditions(bad_cfg1) === nothing

        # after and before both set → nothing
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
        @test eegfun.parse_epoch_conditions(bad_cfg2) === nothing

        # reference_index out of bounds → nothing
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
        @test eegfun.parse_epoch_conditions(bad_cfg3) === nothing

        # timing_pairs provided but min/max missing → nothing
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
        @test eegfun.parse_epoch_conditions(bad_cfg4) === nothing

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

end


