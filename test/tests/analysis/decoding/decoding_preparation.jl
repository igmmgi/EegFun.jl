# Tests for trial balancing, pseudotrial creation, and MVPA data preparation.

using Test
using DataFrames
using EegFun

@testset "Decoding Preparation" begin
    # Helper to generate mock participant epoch files
    function make_mock_epoch(sub::Int, cond::Int; fs = 256, n_samples = 100, n_channels = 3)
        ep = EegFun.create_test_epoch_data(
            n = n_samples,
            fs = fs,
            n_epochs = 4,
            n_channels = n_channels,
            condition = cond,
            seed = sub * 10 + cond,
        )
        return EegFun.EpochData(
            "cond_$(cond)_sub_$(sub).jld2",
            cond,
            "condition_$(cond)",
            ep.data,
            ep.layout,
            ep.sample_rate,
            ep.analysis_info,
        )
    end

    p1_c1 = make_mock_epoch(1, 1)
    p2_c1 = make_mock_epoch(2, 1)
    p1_c2 = make_mock_epoch(1, 2)
    p2_c2 = make_mock_epoch(2, 2)

    all_epochs = [p1_c1, p2_c1, p1_c2, p2_c2]

    # 1. Standard preparation
    prepared = EegFun.prepare_decoding(all_epochs)
    @test length(prepared) == 2
    @test length(prepared[1]) == 2
    @test length(prepared[2]) == 2
    @test prepared[1][1].condition == 1
    @test prepared[1][2].condition == 2
    @test prepared[2][1].condition == 1
    @test prepared[2][2].condition == 2

    # 2. Channel and interval selections
    ch1 = EegFun.channel_labels(p1_c1)[1]
    prep_sub = EegFun.prepare_decoding(all_epochs; channel_selection = EegFun.channels(ch1), interval_selection = EegFun.times(0.0, 0.2))
    @test length(EegFun.channel_labels(prep_sub[1][1])) == 1
    @test EegFun.channel_labels(prep_sub[1][1])[1] == ch1

    # 3. TimeFreqEpochData preparation
    function make_mock_tf_epoch(sub::Int, cond::Int)
        return EegFun.create_test_tf_epoch_data(participant = sub, file = "tf_cond_$(cond)_sub_$(sub).jld2", condition = cond)
    end

    tf_p1_c1 = make_mock_tf_epoch(1, 1)
    tf_p2_c1 = make_mock_tf_epoch(2, 1)
    tf_p1_c2 = make_mock_tf_epoch(1, 2)
    tf_p2_c2 = make_mock_tf_epoch(2, 2)

    tf_prepared = EegFun.prepare_decoding([tf_p1_c1, tf_p2_c1, tf_p1_c2, tf_p2_c2])
    @test length(tf_prepared) == 2
    @test length(tf_prepared[1]) == 2
    @test tf_prepared[1][1] isa EegFun.TimeFreqEpochData

    # 4. Error cases
    # Empty epochs
    @test_throws EegFun.EegFunError EegFun.prepare_decoding(EegFun.EpochData[])

    # Only 1 condition
    @test_throws EegFun.EegFunError EegFun.prepare_decoding([p1_c1, p2_c1])

    # Participant mismatch between conditions
    p3_c2 = make_mock_epoch(3, 2)
    @test_throws EegFun.EegFunError EegFun.prepare_decoding([p1_c1, p2_c1, p1_c2, p3_c2])

    # Inconsistent structure between conditions (different sample rate)
    p2_c2_diff_fs = make_mock_epoch(2, 2; fs = 512)
    @test_throws EegFun.EegFunError EegFun.prepare_decoding([p1_c1, p2_c1, p1_c2, p2_c2_diff_fs])

    # Empty channel selection
    @test_throws EegFun.EegFunError EegFun.prepare_decoding(all_epochs; channel_selection = EegFun.channels(Symbol[]))

    # 5. Non-1-based / named condition_selection matches .condition / .condition_name rather than positional index
    p1_c3 = make_mock_epoch(1, 3)
    p2_c3 = make_mock_epoch(2, 3)
    p1_c2.condition_name = "faces"
    p2_c2.condition_name = "faces"
    p1_c3.condition_name = "cars"
    p2_c3.condition_name = "cars"
    prep_23 = EegFun.prepare_decoding([p1_c2, p2_c2, p1_c3, p2_c3]; condition_selection = EegFun.conditions([2, 3]))
    @test length(prep_23) == 2 && [ep.condition for ep in prep_23[1]] == [2, 3]
    prep_named = EegFun.prepare_decoding([p1_c2, p2_c2, p1_c3, p2_c3]; condition_selection = EegFun.conditions(["cars", "faces"]))
    @test length(prep_named) == 2 && [ep.condition_name for ep in prep_named[1]] == ["cars", "faces"]

    # 6. In-memory EpochData sharing placeholder .file ("test_data") do not collapse N participants to 1
    inmem_no_col = vcat(
        [EegFun.create_test_epoch_data(n = 50, fs = 100, n_epochs = 4, condition = 1, seed = p) for p = 1:4],
        [EegFun.create_test_epoch_data(n = 50, fs = 100, n_epochs = 4, condition = 2, seed = 10 + p) for p = 1:4],
    )
    prep_pos = EegFun.prepare_decoding(inmem_no_col)
    @test length(prep_pos) == 4
    @test prep_pos[3][1] === inmem_no_col[3] || prep_pos[3][1].data[1].Ch1 == inmem_no_col[3].data[1].Ch1

    # In-memory EpochData with :participant metadata column align shuffled conditions and reject mismatches
    function make_inmem_with_col(sub::Int, cond::Int)
        ep = EegFun.create_test_epoch_data(n = 50, fs = 100, n_epochs = 3, condition = cond, seed = 100 * sub + cond)
        for df in ep.data
            insertcols!(df, 1, :participant => fill(sub, size(df, 1)))
        end
        return ep
    end
    col_c1 = [make_inmem_with_col(1, 1), make_inmem_with_col(2, 1), make_inmem_with_col(3, 1)]
    col_c2_shuf = [make_inmem_with_col(3, 2), make_inmem_with_col(1, 2), make_inmem_with_col(2, 2)]
    prep_col = EegFun.prepare_decoding(vcat(col_c1, col_c2_shuf))
    @test length(prep_col) == 3
    @test [ep[1].data[1].participant[1] for ep in prep_col] == [1, 2, 3]
    @test [ep[2].data[1].participant[1] for ep in prep_col] == [1, 2, 3]
    @test prep_col[1][2].data[1].Ch1 == col_c2_shuf[2].data[1].Ch1
    @test_throws EegFun.EegFunError EegFun.prepare_decoding(
        vcat(col_c1, [make_inmem_with_col(1, 2), make_inmem_with_col(2, 2), make_inmem_with_col(99, 2)]),
    )

    # Duplicate participant IDs within a condition throw instead of silently dropping data
    p1_c1_dup = make_mock_epoch(1, 1)
    p1_c2_dup = make_mock_epoch(1, 2)
    @test_throws EegFun.EegFunError EegFun.prepare_decoding([p1_c1, p1_c1_dup, p2_c1, p1_c2, p1_c2_dup, p2_c2])
end
