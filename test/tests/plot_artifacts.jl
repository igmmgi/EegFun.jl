using Test
using DataFrames
using Makie
using EegFun

@testset "plot_artifacts" begin
    dat, bad_idx = EegFun.create_test_epoch_data_with_artifacts()

    # Needs to be marked first
    EegFun.is_extreme_value!(dat, 100)
    artifacts = EegFun.detect_bad_epochs_automatic(dat)

    @testset "plot_artifact_rejection" begin
        dat_rejected = EegFun.reject_epochs(dat, artifacts)
        res = EegFun.plot_artifact_rejection(dat, dat_rejected, artifacts; colormap = :plasma, display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
        @test length(res.axes) == 2
        @test any(b -> b isa Label && occursin("Epoch 1 / $(EegFun.n_epochs(dat))", b.text[]), res.fig.content)
        @test any(p -> p isa Makie.VLines || p isa Makie.HLines, res.axes[1].scene.plots)
    end

    @testset "plot_artifact_detection" begin
        res = EegFun.plot_artifact_detection(dat, artifacts; colormap = :viridis, display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
        @test length(res.axes) == 1
        @test any(b -> b isa Label && occursin("Epoch 1 / $(EegFun.n_epochs(dat))", b.text[]), res.fig.content)
        @test any(p -> p isa Makie.VLines || p isa Makie.HLines, res.axes[1].scene.plots)
    end

    @testset "plot_artifact_repair" begin
        dat_repair, bad_idx2 = EegFun.create_test_epoch_data_with_artifacts()
        EegFun.is_extreme_value!(dat_repair, 100)
        artifacts_repair = EegFun.detect_bad_epochs_automatic(dat_repair)
        dat_repaired = EegFun.repair_artifacts(dat_repair, artifacts_repair)

        res = EegFun.plot_artifact_repair(dat_repair, dat_repaired, artifacts_repair; display_plot = false)
        @test res.fig isa Figure
        @test res.axes isa Vector{Axis}
        @test length(res.axes) == 2
        @test any(b -> b isa Label && occursin("Epoch 1 / $(EegFun.n_epochs(dat_repair))", b.text[]), res.fig.content)
        @test any(p -> p isa Makie.VLines || p isa Makie.HLines, res.axes[1].scene.plots)
    end

    @testset "plot_joint_probability" begin
        jp_df = DataFrame(channel = [:Ch1, :Ch2, :Ch3], jp = [0.2, 2.4, -0.5])
        res_jp = EegFun.plot_joint_probability(jp_df; plot_title = nothing, display_plot = false)
        @test res_jp.fig isa Figure
        @test res_jp.axes isa Vector{Axis}
        @test length(res_jp.axes) == 1
        @test res_jp.axes[1].title[] == ""
    end
end
