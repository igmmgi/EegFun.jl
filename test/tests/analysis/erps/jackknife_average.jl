# Tests for leave-one-out jackknife subaveraging and Smulders pseudovalues.

using Test
using JLD2
using DataFrames
using Statistics
using EegFun


@testset "Jackknife Average" begin
    @testset "In-memory jackknife_average function" begin
        @testset "Basic jackknife averaging" begin
            # Create test ERP data for 4 participants
            erps = [EegFun.create_test_erp_data(participant = i, condition = 1) for i = 1:4]

            # Create jackknife averages
            jackknife_results = EegFun.jackknife_average(erps)

            # Should have 4 jackknife averages (one per participant)
            @test length(jackknife_results) == 4

            # Each jackknife should be an ErpData object
            for jk in jackknife_results
                @test jk isa EegFun.ErpData
                @test jk.sample_rate == 1000.0
                @test nrow(jk.data) == 2501
            end
        end

        @testset "Jackknife calculation verification" begin
            # Create simple test data with known values
            erps = [EegFun.create_test_erp_data(participant = i, condition = 1, fs = 50) for i = 1:3]

            jackknife_results = EegFun.jackknife_average(erps)

            # Jackknife 1 should be average of participants 2 and 3
            # Jackknife 2 should be average of participants 1 and 3
            # Jackknife 3 should be average of participants 1 and 2

            # Verify first jackknife (excluding participant 1)
            jk1 = jackknife_results[1]
            expected_ch = (erps[2].data.Ch1 .+ erps[3].data.Ch1) ./ 2
            @test all(abs.(jk1.data.Ch1 .- expected_ch) .< 1e-10)

            # Verify second jackknife (excluding participant 2)
            jk2 = jackknife_results[2]
            expected_ch = (erps[1].data.Ch1 .+ erps[3].data.Ch1) ./ 2
            @test all(abs.(jk2.data.Ch1 .- expected_ch) .< 1e-10)

            # Verify third jackknife (excluding participant 3)
            jk3 = jackknife_results[3]
            expected_ch = (erps[1].data.Ch1 .+ erps[2].data.Ch1) ./ 2
            @test all(abs.(jk3.data.Ch1 .- expected_ch) .< 1e-10)
        end

        @testset "Jackknife with different channels" begin
            # Test that all channels are processed correctly
            erps = [EegFun.create_test_erp_data(participant = i, condition = 1, fs = 50, n_channels = 3) for i = 1:3]

            jackknife_results = EegFun.jackknife_average(erps)

            # Verify all channels are present in jackknife results
            for jk in jackknife_results
                @test hasproperty(jk.data, :Ch1)
                @test hasproperty(jk.data, :Ch2)
                @test hasproperty(jk.data, :Ch3)
            end
        end

        @testset "Error handling: insufficient participants" begin
            # Only 1 participant - should throw error
            erps = [EegFun.create_test_erp_data(participant = 1, condition = 1)]
            @test_throws Exception EegFun.jackknife_average(erps)
        end

        @testset "Error handling: mismatched sample rates" begin
            # Create ERPs with different sample rates
            erp1 = EegFun.create_test_erp_data(participant = 1, condition = 1)
            erp2 = EegFun.create_test_erp_data(participant = 2, condition = 1)
            erp2 = EegFun.ErpData(
                erp2.file,
                erp2.condition,
                erp2.condition_name,
                erp2.data,
                erp2.layout,
                500,
                erp2.analysis_info,
                erp2.n_epochs,
            )  # Different sample rate

            @test_throws Exception EegFun.jackknife_average([erp1, erp2])
        end

        @testset "Error handling: mismatched time points" begin
            # Create ERPs with different numbers of time points
            erp1 = EegFun.create_test_erp_data(participant = 1, condition = 1, fs = 100)
            erp2 = EegFun.create_test_erp_data(participant = 2, condition = 1, fs = 50)

            @test_throws Exception EegFun.jackknife_average([erp1, erp2])
        end

        @testset "Metadata preservation" begin
            erps = [EegFun.create_test_erp_data(participant = i, condition = 1) for i = 1:3]

            jackknife_results = EegFun.jackknife_average(erps)

            # Check that metadata is preserved
            for (i, jk) in enumerate(jackknife_results)
                @test hasproperty(jk.data, :time)
                # condition and condition_name are now in struct, not DataFrame
                @test hasproperty(jk, :condition)
                @test hasproperty(jk, :condition_name)
                @test occursin("jackknife", jk.condition_name)
            end
        end

        @testset "n_epochs calculation" begin
            # Create ERPs with specific n_epochs
            erps = [EegFun.create_test_erp_data(participant = i, condition = 1) for i = 1:4]
            for (i, erp) in enumerate(erps)
                erps[i] = EegFun.ErpData(
                    erp.file,
                    erp.condition,
                    erp.condition_name,
                    erp.data,
                    erp.layout,
                    erp.sample_rate,
                    erp.analysis_info,
                    i * 10,
                )
            end

            jackknife_results = EegFun.jackknife_average(erps)

            # Jackknife 1 excludes participant 1 (10 epochs), includes 2 (20) + 3 (30) + 4 (40) = 90
            @test jackknife_results[1].n_epochs == 90

            # Jackknife 2 excludes participant 2 (20 epochs), includes 1 (10) + 3 (30) + 4 (40) = 80
            @test jackknife_results[2].n_epochs == 80

            # Jackknife 3 excludes participant 3 (30 epochs), includes 1 (10) + 2 (20) + 4 (40) = 70
            @test jackknife_results[3].n_epochs == 70

            # Jackknife 4 excludes participant 4 (40 epochs), includes 1 (10) + 2 (20) + 3 (30) = 60
            @test jackknife_results[4].n_epochs == 60
        end
    end

    @testset "Batch jackknife_average function" begin
        mktempdir() do test_dir
            for p = 1:4
                jldsave(joinpath(test_dir, "$(p)_single_lrp.jld2"); data = EegFun.create_test_erp_data(participant = p, condition = 1))
                jldsave(
                    joinpath(test_dir, "$(p)_multi_lrp.jld2");
                    data = [EegFun.create_test_erp_data(participant = p, condition = c) for c in (1, 2)],
                )
            end

            # Single-condition & multi-condition batch jackknife, logging & auto-dir
            out_single = joinpath(test_dir, "jk_single")
            EegFun.jackknife_average("single_lrp"; input_dir = test_dir, output_dir = out_single)
            for p = 1:4
                @test isfile(joinpath(out_single, "$(p)_single_lrp.jld2"))
            end
            jk1 = load(joinpath(out_single, "1_single_lrp.jld2"), "data")
            lrp2 = load(joinpath(test_dir, "2_single_lrp.jld2"), "data")
            lrp3 = load(joinpath(test_dir, "3_single_lrp.jld2"), "data")
            lrp4 = load(joinpath(test_dir, "4_single_lrp.jld2"), "data")
            @test jk1 isa EegFun.ErpData && nrow(jk1.data) == 2501 && jk1.sample_rate == 1000.0
            @test jk1.data.Ch1 ≈ (lrp2.data.Ch1 .+ lrp3.data.Ch1 .+ lrp4.data.Ch1) ./ 3 atol = 1e-10
            log_txt = read(joinpath(out_single, "jackknife.log"), String)
            @test occursin("Batch jackknife averaging started", log_txt) && occursin("Jackknife averaging complete", log_txt)

            out_multi = joinpath(test_dir, "jk_multi")
            EegFun.jackknife_average("multi_lrp"; input_dir = test_dir, output_dir = out_multi)
            jkm1 = load(joinpath(out_multi, "1_multi_lrp.jld2"), "data")
            @test jkm1 isa Vector{EegFun.ErpData} && [d.condition for d in jkm1] == [1, 2]

            # Participant & condition filtering, auto-dir & error cases
            out_pf = joinpath(test_dir, "jk_pf")
            EegFun.jackknife_average(
                "multi_lrp";
                input_dir = test_dir,
                participant_selection = EegFun.participants([1, 2, 3]),
                output_dir = out_pf,
            )
            @test isfile(joinpath(out_pf, "1_multi_lrp.jld2")) && !isfile(joinpath(out_pf, "4_multi_lrp.jld2"))

            out_cf = joinpath(test_dir, "jk_cf")
            EegFun.jackknife_average("multi_lrp"; input_dir = test_dir, condition_selection = EegFun.conditions([1]), output_dir = out_cf)
            jk_c1 = load(joinpath(out_cf, "1_multi_lrp.jld2"), "data")
            @test jk_c1 isa EegFun.ErpData && jk_c1.condition == 1

            EegFun.jackknife_average("multi_lrp"; input_dir = test_dir)
            @test isdir(joinpath(test_dir, "jackknife_multi_lrp"))
            @test isnothing(EegFun.jackknife_average("nonexistent"; input_dir = test_dir))
            @test isnothing(EegFun.jackknife_average("multi_lrp"; input_dir = test_dir, participant_selection = EegFun.participants([1])))
        end

        @testset "Excluded file metadata, missing-condition participant alignment, and pseudovalues" begin
            e1 = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 2)
            e2 = EegFun.create_test_erp_data(participant = 2, condition = 1, n_channels = 2)
            e3 = EegFun.create_test_erp_data(participant = 3, condition = 1, n_channels = 2)
            e1.file, e2.file, e3.file = "sub01_custom", "sub02_custom", "sub03_custom"
            e1.data.Ch1 .= 10.0
            e2.data.Ch1 .= 20.0
            e3.data.Ch1 .= 60.0

            jk_sub = EegFun.jackknife_average([e1, e2, e3]; pseudovalues = false)
            @test [d.file for d in jk_sub] == ["sub01_custom", "sub02_custom", "sub03_custom"]
            @test jk_sub[1].data.Ch1[1] ≈ 40.0

            jk_pv = EegFun.jackknife_average([e1, e2, e3]; pseudovalues = true)
            @test [d.file for d in jk_pv] == ["sub01_custom", "sub02_custom", "sub03_custom"]
            @test [d.data.Ch1[1] for d in jk_pv] ≈ [10.0, 20.0, 60.0] atol = 1e-10
            @test mean([d.data.Ch1[1] for d in jk_pv]) ≈ 30.0 atol = 1e-10

            mktempdir() do align_dir
                out_align = joinpath(align_dir, "out")
                out_align_pv = joinpath(align_dir, "out_pv")
                p1_c1 = EegFun.create_test_erp_data(participant = 1, condition = 1, n_channels = 2)
                p1_c2 = EegFun.create_test_erp_data(participant = 1, condition = 2, n_channels = 2)
                p2_c2 = EegFun.create_test_erp_data(participant = 2, condition = 2, n_channels = 2)
                p3_c1 = EegFun.create_test_erp_data(participant = 3, condition = 1, n_channels = 2)
                p3_c2 = EegFun.create_test_erp_data(participant = 3, condition = 2, n_channels = 2)
                p1_c1.data.Ch1 .= 100.0
                p3_c1.data.Ch1 .= 200.0

                jldsave(joinpath(align_dir, "1_align.jld2"); data = [p1_c1, p1_c2])
                jldsave(joinpath(align_dir, "2_align.jld2"); data = [p2_c2])
                jldsave(joinpath(align_dir, "3_align.jld2"); data = [p3_c1, p3_c2])

                EegFun.jackknife_average("align"; input_dir = align_dir, output_dir = out_align)
                EegFun.jackknife_average("align"; input_dir = align_dir, output_dir = out_align_pv, pseudovalues = true)
                @test load(joinpath(out_align_pv, "1_align.jld2"), "data")[1].data.Ch1[1] ≈ 100.0 atol = 1e-10
                out1 = load(joinpath(out_align, "1_align.jld2"), "data")
                out2 = load(joinpath(out_align, "2_align.jld2"), "data")
                out3 = load(joinpath(out_align, "3_align.jld2"), "data")

                @test out2 isa EegFun.ErpData && out2.condition == 2
                @test out1 isa Vector{EegFun.ErpData} && length(out1) == 2
                @test out3 isa Vector{EegFun.ErpData} && length(out3) == 2
                @test only(filter(x->x.condition==1, out3)).data.Ch1[1] ≈ 100.0 atol = 1e-10
            end
        end
    end

    @testset "In-memory multi-condition grouping and condition_selection" begin
        multi_erps = EegFun.ErpData[]
        for (p, v1, v2) in ((1, 10.0, 1000.0), (2, 20.0, 2000.0), (3, 30.0, 3000.0))
            for (c, cname, v) in ((1, "cond_1", v1), (2, "cond_2", v2))
                e = EegFun.create_test_erp_data(participant = p, condition = c, n_channels = 2)
                e.condition_name = cname
                e.data.Ch1 .= v
                push!(multi_erps, e)
            end
        end

        jk_all = EegFun.jackknife_average(multi_erps)
        @test length(jk_all) == 6
        @test filter(x->x.condition==1, jk_all)[1].data.Ch1[1] ≈ 25.0 atol = 1e-10
        @test filter(x->x.condition==2, jk_all)[1].data.Ch1[1] ≈ 2500.0 atol = 1e-10

        jk_sel2 = EegFun.jackknife_average(multi_erps; condition_selection = EegFun.conditions(2))
        @test length(jk_sel2) == 3 && all(x -> x.condition == 2, jk_sel2) && isapprox(jk_sel2[1].data.Ch1[1], 2500.0; atol = 1e-10)

        jk_named = EegFun.jackknife_average(multi_erps; condition_selection = EegFun.conditions("cond_2"))
        @test length(jk_named) == 3 && isapprox(jk_named[2].data.Ch1[1], 2000.0; atol = 1e-10)
    end
end
