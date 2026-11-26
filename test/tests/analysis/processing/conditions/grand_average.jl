# Tests for multi-participant grand averaging of ERP and time-frequency data.

using Test
using JLD2
using DataFrames
using EegFun


@testset "Batch Grand Average" begin
    mktempdir() do test_dir
        for p = 1:4
            jldsave(
                joinpath(test_dir, "$(p)_erps_unrejected.jld2");
                data = [EegFun.create_test_erp_data(participant = p, condition = c) for c = 1:3],
            )
        end

        # Basic run, structure, logging & auto-generated directory
        out_dir = joinpath(test_dir, "grand_average")
        EegFun.grand_average("erps_unrejected"; input_dir = test_dir, output_dir = out_dir)
        gas = load(joinpath(out_dir, "grand_average_erps_unrejected.jld2"), "data")
        @test length(gas) == 3
        for ga in gas
            @test ga isa EegFun.ErpData && ga.condition in 1:3
            @test ga.condition_name == "grand_avg_condition_$(ga.condition)"
            @test ga.sample_rate == 1000.0 && nrow(ga.data) == 2501
            @test "time" in names(ga.data) && "condition" ∉ names(ga.data)
            @test all(ch -> ch in names(ga.data), ("Ch1", "Ch2", "Ch3"))
        end
        log_txt = read(joinpath(out_dir, "grand_average.log"), String)
        @test occursin("Batch grand averaging started", log_txt) && occursin("Grand averaging complete", log_txt)

        EegFun.grand_average("erps_unrejected"; input_dir = test_dir)
        @test isdir(joinpath(test_dir, "grand_average_erps_unrejected"))

        # Participant & condition selection + analytical mean verification
        out_p = joinpath(test_dir, "ga_p123")
        EegFun.grand_average(
            "erps_unrejected";
            input_dir = test_dir,
            participant_selection = EegFun.participants([1, 2, 3]),
            output_dir = out_p,
        )
        @test all(ga.n_epochs == 30 for ga in load(joinpath(out_p, "grand_average_erps_unrejected.jld2"), "data"))

        out_c = joinpath(test_dir, "ga_c12")
        EegFun.grand_average("erps_unrejected"; input_dir = test_dir, condition_selection = EegFun.conditions([1, 2]), output_dir = out_c)
        @test sort([ga.condition for ga in load(joinpath(out_c, "grand_average_erps_unrejected.jld2"), "data")]) == [1, 2]

        out_v = joinpath(test_dir, "ga_verify")
        EegFun.grand_average(
            "erps_unrejected";
            input_dir = test_dir,
            participant_selection = EegFun.participants([1, 2]),
            condition_selection = EegFun.conditions([1]),
            output_dir = out_v,
        )
        e1 = load(joinpath(test_dir, "1_erps_unrejected.jld2"), "data")[1]
        e2 = load(joinpath(test_dir, "2_erps_unrejected.jld2"), "data")[1]
        ga_v = only(load(joinpath(out_v, "grand_average_erps_unrejected.jld2"), "data"))
        @test ga_v.n_epochs == e1.n_epochs + e2.n_epochs
        for ch in (:Ch1, :Ch2, :Ch3)
            @test ga_v.data[!, ch] ≈ (e1.data[!, ch] .+ e2.data[!, ch]) ./ 2 atol = 1e-10
        end

        # Error & edge cases (missing dir, nomatch, single participant, unequal condition counts)
        @test_throws Exception EegFun.grand_average("erps_unrejected"; input_dir = joinpath(test_dir, "nonexistent"))
        @test isnothing(EegFun.grand_average("nonexistent_pattern"; input_dir = test_dir))
        @test isnothing(EegFun.grand_average("erps_unrejected"; input_dir = test_dir, participant_selection = EegFun.participants([1])))
        @test isnothing(EegFun.grand_average("erps_unrejected"; input_dir = test_dir, participant_selection = EegFun.participants(999)))

        diff_dir = mkpath(joinpath(test_dir, "diff_counts"))
        for p = 1:3
            conds = p <= 2 ? (1, 2) : (1,)
            jldsave(
                joinpath(diff_dir, "$(p)_erps_unrejected.jld2");
                data = [EegFun.create_test_erp_data(participant = p, condition = c) for c in conds],
            )
        end
        out_dc = joinpath(diff_dir, "out")
        EegFun.grand_average("erps_unrejected"; input_dir = diff_dir, output_dir = out_dc)
        sorted_dc = sort(load(joinpath(out_dc, "grand_average_erps_unrejected.jld2"), "data"), by = x -> x.condition)
        @test [g.n_epochs for g in sorted_dc] == [30, 20]
    end
end

@testset "TF Grand Average" begin

    @testset "Basic in-memory TF grand average" begin
        # Create TF data for multiple "participants" with same condition
        tfs = [
            EegFun.create_test_tf_data(participant = 1, condition = 1),
            EegFun.create_test_tf_data(participant = 2, condition = 1),
            EegFun.create_test_tf_data(participant = 3, condition = 1),
            EegFun.create_test_tf_data(participant = 1, condition = 2),
            EegFun.create_test_tf_data(participant = 2, condition = 2),
            EegFun.create_test_tf_data(participant = 3, condition = 2),
        ]

        grand_avgs = EegFun.grand_average(tfs)

        @test length(grand_avgs) == 2
        @test all(ga -> ga isa EegFun.TimeFreqData, grand_avgs)

        # Verify condition numbering
        conds = sort([ga.condition for ga in grand_avgs])
        @test conds == [1, 2]
    end

    @testset "Structure preservation" begin
        tfs = [EegFun.create_test_tf_data(participant = 1, condition = 1), EegFun.create_test_tf_data(participant = 2, condition = 1)]

        grand_avgs = EegFun.grand_average(tfs)

        @test length(grand_avgs) == 1
        ga = grand_avgs[1]

        # Verify structure matches input
        @test ga.sample_rate == tfs[1].sample_rate
        @test ga.method == tfs[1].method
        @test nrow(ga.data_power) == nrow(tfs[1].data_power)
        @test nrow(ga.data_phase) == nrow(tfs[1].data_phase)

        # Same channels
        @test EegFun.channel_labels(ga) == EegFun.channel_labels(tfs[1])

        # Same freq/time structure
        @test sort(unique(ga.data_power.freq)) == sort(unique(tfs[1].data_power.freq))
        @test sort(unique(ga.data_power.time)) == sort(unique(tfs[1].data_power.time))
    end

    @testset "Power data is averaged" begin
        # Create TF data with known power offsets
        tf1 = EegFun.create_test_tf_data(participant = 1, condition = 1, power_offset = 4.0, noise = 0.0)
        tf2 = EegFun.create_test_tf_data(participant = 2, condition = 1, power_offset = 8.0, noise = 0.0)

        grand_avgs = EegFun.grand_average([tf1, tf2])

        @test length(grand_avgs) == 1
        ga = grand_avgs[1]

        # Grand average power should be mean of 4.0 and 8.0 = 6.0
        for ch in EegFun.channel_labels(ga)
            @test all(abs.(ga.data_power[!, ch] .- 6.0) .< 1e-10)
        end
    end

    @testset "Condition selection" begin
        tfs = [
            EegFun.create_test_tf_data(participant = 1, condition = 1),
            EegFun.create_test_tf_data(participant = 2, condition = 1),
            EegFun.create_test_tf_data(participant = 1, condition = 2),
            EegFun.create_test_tf_data(participant = 2, condition = 2),
            EegFun.create_test_tf_data(participant = 1, condition = 3),
            EegFun.create_test_tf_data(participant = 2, condition = 3),
        ]

        grand_avgs = EegFun.grand_average(tfs, condition_selection = EegFun.conditions([1, 3]))

        @test length(grand_avgs) == 2
        conds = sort([ga.condition for ga in grand_avgs])
        @test conds == [1, 3]
    end

    @testset "Single participant skipped" begin
        tfs = [
            EegFun.create_test_tf_data(participant = 1, condition = 1),
            # Only 1 participant for condition 1 → should skip
        ]

        grand_avgs = EegFun.grand_average(tfs)
        @test isempty(grand_avgs)
    end

    @testset "Empty input throws" begin
        @test_throws Exception EegFun.grand_average(EegFun.TimeFreqData[])
        @test_throws Exception EegFun.grand_average(EegFun.ErpData[])
    end

    @testset "In-memory ERP grand average" begin
        erps = [
            EegFun.create_test_erp_data(participant = 1, condition = 1),
            EegFun.create_test_erp_data(participant = 2, condition = 1),
            EegFun.create_test_erp_data(participant = 1, condition = 2),
            EegFun.create_test_erp_data(participant = 2, condition = 2),
        ]
        grand_avgs = EegFun.grand_average(erps)
        @test length(grand_avgs) == 2
        @test all(ga -> ga isa EegFun.ErpData, grand_avgs)
        @test sort([ga.condition for ga in grand_avgs]) == [1, 2]
    end

    @testset "Batch TF grand average" begin
        tf_test_dir = mktempdir()
        for p = 1:3
            tfs = [EegFun.create_test_tf_data(participant = p, condition = 1), EegFun.create_test_tf_data(participant = p, condition = 2)]
            jldsave(joinpath(tf_test_dir, "$(p)_tf_data.jld2"); data = tfs)
        end

        tf_out = joinpath(tf_test_dir, "tf_grand_avg")
        EegFun.grand_average("tf_data", input_dir = tf_test_dir, output_dir = tf_out)

        out_file = joinpath(tf_out, "grand_average_tf_data.jld2")
        @test isfile(out_file)
        loaded = load(out_file, "data")
        @test length(loaded) == 2
        @test all(x -> x isa EegFun.TimeFreqData, loaded)
        @test sort([x.condition for x in loaded]) == [1, 2]

        rm(tf_test_dir, recursive = true)
    end

    @testset "Batch ERP grand_average matches .condition field on non-contiguous conditions" begin
        nc_dir = mktempdir()
        nc_out = joinpath(nc_dir, "out")
        for p = 1:2
            e2 = EegFun.create_test_erp_data(participant = p, condition = 2)
            e5 = EegFun.create_test_erp_data(participant = p, condition = 5)
            e2.data.Ch1 .= 10.0 * p
            e5.data.Ch1 .= 100.0 * p
            jldsave(joinpath(nc_dir, "$(p)_erps_nc.jld2"); data = [e2, e5])
        end
        EegFun.grand_average("erps_nc"; input_dir = nc_dir, output_dir = nc_out)
        gas = load(joinpath(nc_out, "grand_average_erps_nc.jld2"), "data")
        @test sort([g.condition for g in gas]) == [2, 5]
        ga2 = only(filter(g -> g.condition == 2, gas))
        ga5 = only(filter(g -> g.condition == 5, gas))
        @test ga2.data.Ch1[1] ≈ 15.0
        @test ga5.data.Ch1[1] ≈ 150.0
        rm(nc_dir, recursive = true, force = true)
    end
end
