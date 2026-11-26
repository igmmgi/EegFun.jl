# Tests for merging and combining experimental conditions.

using Test
using DataFrames
using JLD2
using Statistics
using EegFun

@testset "Batch Combine Conditions" begin
    isdefined(@__MODULE__, :run_batch_jld2_contract_tests!) ||
        Base.include(@__MODULE__, joinpath(pkgdir(EegFun), "test", "tests", "test_helpers.jl"))

    run_batch_jld2_contract_tests!(
        (pat; kwargs...) -> EegFun.condition_combine(pat, [[1, 2], [3, 4]]; kwargs...);
        make_data = _ -> EegFun.create_test_epoch_data_vector(conditions = 1:4),
        pattern = "epochs_unrejected",
        n_participants = 3,
        supports_condition_selection = false,
        log_file = "condition_combine.log",
        log_snippets = ["condition_combine", "epochs_unrejected", "Found"],
        expected_auto_dir = "combined_epochs_unrejected_1-2_3-4",
    )

    mktempdir() do test_dir
        orig = EegFun.create_test_epoch_data_vector(conditions = 1:4)
        jldsave(joinpath(test_dir, "1_epochs_unrejected.jld2"); data = orig)

        # Grouping shapes: standard, single, overlapping, duplicate removal, all-in-one, repeated groups
        for (tag, groups, expected_counts) in (
            ("basic", [[1, 2], [3, 4]], [20, 20]),
            ("single", [[1], [2], [3], [4]], [10, 10, 10, 10]),
            ("overlap", [[1, 2], [2, 3], [3, 4]], [20, 20, 20]),
            ("dup", [[1, 1, 2], [3, 3, 4]], [20, 20]),
            ("dup1", [[1, 1, 1], [2, 3]], [10, 20]),
            ("repeat", [[1], [1], [2]], [10, 10, 10]),
            ("all", [[1, 2, 3, 4]], [40]),
            ("many", [[1], [2], [3], [4], [1], [2], [3], [4]], fill(10, 8)),
        )
            odir = joinpath(test_dir, "grp_$tag")
            res = EegFun.condition_combine("epochs_unrejected", groups; input_dir = test_dir, output_dir = odir)
            @test res.success == 1 && res.errors == 0
            comb = load(joinpath(odir, "1_epochs_unrejected.jld2"), "data")
            @test [length(c.data) for c in comb] == expected_counts
            @test [c.condition for c in comb] == collect(1:length(expected_counts))
            if tag == "basic"
                @test comb[1].layout.data == orig[1].layout.data && comb[1].sample_rate == orig[1].sample_rate
                @test hasproperty(comb[1].data[1], :time) && hasproperty(comb[1].data[1], :epoch)
            end
        end

        # Validation & out-of-range / corrupt / empty-data failure handling
        @test_throws Exception EegFun.condition_combine("epochs_unrejected", []; input_dir = test_dir)
        @test_throws Exception EegFun.condition_combine("epochs_unrejected", [[1, 2], []]; input_dir = test_dir)
        @test_throws MethodError EegFun.condition_combine("epochs_unrejected", [["1", "2"], [3, 4]]; input_dir = test_dir)

        for (tag, bad_groups) in (("oor", [[1, 5]]), ("neg", [[1, -1], [2, 3]]), ("zero", [[0, 1], [2, 3]]), ("large", [[1000, 1001]]))
            res_b =
                EegFun.condition_combine("epochs_unrejected", bad_groups; input_dir = test_dir, output_dir = joinpath(test_dir, "bad_$tag"))
            @test res_b.success == 0 && res_b.errors == 1
        end

        corrupt_dir = mkpath(joinpath(test_dir, "corrupt"))
        write(joinpath(corrupt_dir, "1_epochs_corrupt.jld2"), "not a valid JLD2 file")
        res_c =
            EegFun.condition_combine("epochs_corrupt", [[1, 2]]; input_dir = corrupt_dir, output_dir = joinpath(test_dir, "out_corrupt"))
        @test res_c.success == 0 && res_c.errors == 1

        empty_ep_dir = mkpath(joinpath(test_dir, "empty_ep"))
        jldsave(joinpath(empty_ep_dir, "1_epochs_empty.jld2"); data = EegFun.EpochData[])
        res_e = EegFun.condition_combine("epochs_empty", [[1, 2]]; input_dir = empty_ep_dir, output_dir = joinpath(test_dir, "out_empty"))
        @test res_e.success == 0 && res_e.errors == 1
    end
end

@testset "In-Memory Combine Conditions" begin

    @testset "Basic combining" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:4)

        combined = EegFun.condition_combine(epochs, [[1, 2], [3, 4]])

        @test length(combined) == 2
        @test combined[1].condition == 1
        @test combined[2].condition == 2

        # Each group combines 2 conditions of 10 epochs each
        @test length(combined[1].data) == 20
        @test length(combined[2].data) == 20
    end

    @testset "Single condition per group" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:3)

        combined = EegFun.condition_combine(epochs, [[1], [2], [3]])

        @test length(combined) == 3
        for i = 1:3
            @test length(combined[i].data) == 10  # unchanged
            @test combined[i].condition == i
        end
    end

    @testset "All conditions in one group" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:4)

        combined = EegFun.condition_combine(epochs, [[1, 2, 3, 4]])

        @test length(combined) == 1
        @test length(combined[1].data) == 40
    end

    @testset "Overlapping groups" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:3)

        combined = EegFun.condition_combine(epochs, [[1, 2], [2, 3]])

        @test length(combined) == 2
        @test length(combined[1].data) == 20
        @test length(combined[2].data) == 20
    end

    @testset "Condition name combining" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:4)

        combined = EegFun.condition_combine(epochs, [[1, 2], [3, 4]])

        # Combined names should join originals with _
        @test contains(combined[1].condition_name, "_")
        @test contains(combined[2].condition_name, "_")
    end

    @testset "Metadata preservation" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2)
        original_sr = epochs[1].sample_rate
        original_layout = epochs[1].layout

        combined = EegFun.condition_combine(epochs, [[1, 2]])

        @test combined[1].sample_rate == original_sr
        @test combined[1].layout.data == original_layout.data
        @test hasproperty(combined[1].data[1], :time)
    end

    @testset "Error: condition out of range" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2)

        @test_throws Exception EegFun.condition_combine(epochs, [[1, 5]])
    end

    @testset "Error: negative condition" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2)

        @test_throws Exception EegFun.condition_combine(epochs, [[1, -1]])
    end

    @testset "Error: zero condition" begin
        epochs = EegFun.create_test_epoch_data_vector(conditions = 1:2)

        @test_throws Exception EegFun.condition_combine(epochs, [[0, 1]])
    end

    @testset "Condition field matching and independent DataFrame copies" begin
        # Non-contiguous .condition numbers (e.g., [2, 5] at vector positions 1 and 2)
        ep_c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 2, n = 20)
        ep_c5 = EegFun.create_test_epoch_data(condition = 5, n_epochs = 3, n = 20)
        orig_val = ep_c2.data[1].Ch1[1]

        combined = EegFun.condition_combine([ep_c2, ep_c5], [[2, 5]])
        @test length(combined) == 1
        @test length(combined[1].data) == 5
        @test [df.epoch[1] for df in combined[1].data] == [1, 2, 3, 4, 5]

        # Mutating combined DataFrame does not modify ep_c2 DataFrame
        combined[1].data[1].Ch1[1] = orig_val + 999.0
        @test ep_c2.data[1].Ch1[1] == orig_val

        # Inconsistent structure throws EegFunError
        ep_mismatch = EegFun.create_test_epoch_data(condition = 6, n_epochs = 2, n = 40)
        @test_throws EegFun.EegFunError EegFun.condition_combine([ep_c2, ep_mismatch], [[2, 6]])

        # Multi-participant in-memory Vector{EpochData} combines per participant
        mp1_c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 2, n = 20)
        mp1_c1.file = "1_epochs.jld2"
        mp1_c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 3, n = 20)
        mp1_c2.file = "1_epochs.jld2"
        mp2_c1 = EegFun.create_test_epoch_data(condition = 1, n_epochs = 4, n = 20)
        mp2_c1.file = "2_epochs.jld2"
        mp2_c2 = EegFun.create_test_epoch_data(condition = 2, n_epochs = 2, n = 20)
        mp2_c2.file = "2_epochs.jld2"
        mp_comb = EegFun.condition_combine([mp1_c1, mp1_c2, mp2_c1, mp2_c2], [[1, 2]])
        @test length(mp_comb) == 2
        @test [EegFun.n_epochs(d) for d in mp_comb] == [5, 6]

        # Batch condition_combine works on filenames without literal "epochs" (e.g. "trials")
        mktempdir() do tmp_trials
            jldsave(joinpath(tmp_trials, "1_trials_clean.jld2"); data = [ep_c2, ep_c5])
            res_tr = EegFun.condition_combine("trials_clean", [[2, 5]]; input_dir = tmp_trials)
            @test res_tr.success == 1
            comb_tr = EegFun.read_data(joinpath(tmp_trials, "combined_trials_clean_2-5", "1_trials_clean.jld2"))
            @test length(comb_tr) == 1 && EegFun.n_epochs(comb_tr[1]) == 5
        end
    end
end
