# Shared batch JLD2 contract helpers and reusable test utilities for the EegFun test suite.

using Test
using JLD2
using DataFrames
using EegFun

function run_batch_jld2_contract_tests!(
    run_fn::Function;
    make_data::Function,
    pattern::String,
    out_filename::Function = p -> "$(p)_$(pattern).jld2",
    n_participants::Int = 2,
    supports_condition_selection::Bool = true,
    log_file::Union{String,Nothing} = nothing,
    log_snippets::Vector{String} = String[],
    expected_auto_dir::Union{String,Nothing} = nothing,
    supports_partial_failure::Bool = true,
)
    mktempdir() do test_dir
        for p = 1:n_participants
            jldsave(joinpath(test_dir, "$(p)_$(pattern).jld2"); data = make_data(p))
        end

        @testset "Basic run, logging & return structure" begin
            out_dir = joinpath(test_dir, "out_basic")
            res = run_fn(pattern; input_dir = test_dir, output_dir = out_dir)
            @test !isnothing(res)
            @test hasfield(typeof(res), :success) && hasfield(typeof(res), :errors)
            @test res.success isa Integer && res.errors isa Integer
            @test res.success == n_participants && res.errors == 0
            @test isdir(out_dir)
            for p = 1:n_participants
                @test isfile(joinpath(out_dir, out_filename(p)))
            end
            if !isnothing(log_file)
                lf = joinpath(out_dir, log_file)
                @test isfile(lf)
                log_txt = read(lf, String)
                for snippet in log_snippets
                    @test occursin(snippet, log_txt)
                end
            end
        end

        @testset "Participant selection" begin
            out_p1 = joinpath(test_dir, "out_p1")
            res1 = run_fn(pattern; input_dir = test_dir, output_dir = out_p1, participant_selection = EegFun.participants(1))
            @test res1.success == 1 && res1.errors == 0
            @test isfile(joinpath(out_p1, out_filename(1)))
            @test !isfile(joinpath(out_p1, out_filename(2)))

            out_pm = joinpath(test_dir, "out_pm")
            sel_ids = [1, n_participants]
            resm = run_fn(pattern; input_dir = test_dir, output_dir = out_pm, participant_selection = EegFun.participants(sel_ids))
            @test resm.success == length(unique(sel_ids)) && resm.errors == 0
            @test isfile(joinpath(out_pm, out_filename(1)))
            @test isfile(joinpath(out_pm, out_filename(n_participants)))
            if n_participants >= 3
                @test !isfile(joinpath(out_pm, out_filename(2)))
            end
        end

        if supports_condition_selection
            @testset "Condition & combined selection" begin
                out_c1 = joinpath(test_dir, "out_c1")
                res_c1 = run_fn(pattern; input_dir = test_dir, output_dir = out_c1, condition_selection = EegFun.conditions(1))
                @test res_c1.success == n_participants
                loaded_c1 = load(joinpath(out_c1, out_filename(1)), "data")
                @test length(loaded_c1) == 1 && loaded_c1[1].condition == 1

                out_c12 = joinpath(test_dir, "out_c12")
                res_c12 = run_fn(pattern; input_dir = test_dir, output_dir = out_c12, condition_selection = EegFun.conditions([1, 2]))
                @test res_c12.success == n_participants
                @test length(load(joinpath(out_c12, out_filename(1)), "data")) == 2

                out_comb = joinpath(test_dir, "out_comb")
                res_comb = run_fn(
                    pattern;
                    input_dir = test_dir,
                    output_dir = out_comb,
                    participant_selection = EegFun.participants(1),
                    condition_selection = EegFun.conditions(1),
                )
                @test res_comb.success == 1
                @test length(load(joinpath(out_comb, out_filename(1)), "data")) == 1
                @test !isfile(joinpath(out_comb, out_filename(2)))

                out_oor = joinpath(test_dir, "out_oor")
                res_oor = run_fn(pattern; input_dir = test_dir, output_dir = out_oor, condition_selection = EegFun.conditions(99))
                @test res_oor.success == n_participants && res_oor.errors == 0
            end
        end

        @testset "Missing dir, empty match, existing dir & overwrite" begin
            @test_throws Exception run_fn(pattern; input_dir = joinpath(test_dir, "nonexistent_dir"))
            @test isnothing(run_fn("nonexistent_$(pattern)_xyz"; input_dir = test_dir, output_dir = joinpath(test_dir, "nomatch")))

            empty_dir = mkpath(joinpath(test_dir, "empty_dir"))
            @test isnothing(run_fn(pattern; input_dir = empty_dir))

            exist_dir = mkpath(joinpath(test_dir, "existing_out"))
            touch(joinpath(exist_dir, "dummy.txt"))
            write(joinpath(exist_dir, out_filename(1)), "stale_placeholder")
            res_ow = run_fn(pattern; input_dir = test_dir, output_dir = exist_dir)
            @test res_ow.success == n_participants
            @test isfile(joinpath(exist_dir, "dummy.txt"))
            @test !isnothing(load(joinpath(exist_dir, out_filename(1)), "data"))
        end

        if supports_partial_failure
            @testset "Partial failures" begin
                pdir = mkpath(joinpath(test_dir, "partial_in"))
                jldsave(joinpath(pdir, "1_$(pattern).jld2"); data = make_data(1))
                jldsave(joinpath(pdir, "2_$(pattern).jld2"); data = "invalid_data")
                pout = joinpath(test_dir, "partial_out")
                res_pf = run_fn(pattern; input_dir = pdir, output_dir = pout)
                @test res_pf.success == 1 && res_pf.errors == 1
                @test isfile(joinpath(pout, out_filename(1)))
                @test !isfile(joinpath(pout, out_filename(2)))
            end
        end

        if !isnothing(expected_auto_dir)
            run_fn(pattern; input_dir = test_dir)
            @test isdir(joinpath(test_dir, expected_auto_dir))
        end
    end
    return nothing
end
