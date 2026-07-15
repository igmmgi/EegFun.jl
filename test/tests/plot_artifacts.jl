using Test
using Makie
using EegFun

@testset "plot_artifacts" begin
    dat, bad_idx = EegFun.create_test_epoch_data_with_artifacts()

    # Needs to be marked first
    EegFun.is_extreme_value!(dat, 100)
    artifacts = EegFun.detect_bad_epochs_automatic(dat)

    @testset "plot_artifact_rejection" begin
        dat_rejected = EegFun.reject_epochs(dat, artifacts)
        fig = EegFun.plot_artifact_rejection(dat, dat_rejected, artifacts; display_plot = false)
        @test fig isa Figure
    end

    @testset "plot_artifact_detection" begin
        fig = EegFun.plot_artifact_detection(dat, artifacts; display_plot = false)
        @test fig isa Figure
    end

    @testset "plot_artifact_repair" begin
        dat_repair, bad_idx2 = EegFun.create_test_epoch_data_with_artifacts()
        EegFun.is_extreme_value!(dat_repair, 100)
        artifacts_repair = EegFun.detect_bad_epochs_automatic(dat_repair)
        dat_repaired = EegFun.repair_artifacts(dat_repair, artifacts_repair)

        # We need continuous data for some repair plots, so we'll test that next
        fig = EegFun.plot_artifact_repair(dat_repair, dat_repaired, artifacts_repair; display_plot = false)
        @test fig isa Figure
    end
end
