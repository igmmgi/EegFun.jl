using Test
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
end
