# Main test runner for EegFun.jl mirroring the `src/` and `ext/` module hierarchy.

using Test
using CairoMakie
using DimensionalData
using MixedModels
using StatsModels
using EegFun

println("Running EegFun.jl Test Suite")
println("="^40)

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    @testset "Method ambiguity" begin
        @test isempty(Test.detect_ambiguities(EegFun; recursive = true))
    end
    Aqua.test_all(
        EegFun;
        ambiguities = false,
        deps_compat = (check_extras = false, ignore = [:Dates, :LinearAlgebra, :Logging, :Printf, :Random, :SparseArrays, :TOML, :Test]),
    )
end

@testset "EegFun" begin

    include("tests/code_organization.jl")

    @testset "Types" begin
        include("tests/types/types.jl")
    end

    @testset "Utilities" begin
        include("tests/utils/batch_utils.jl")
        include("tests/utils/data.jl")
        include("tests/utils/files.jl")
        include("tests/utils/gpu_kernels_cpu.jl")
        include("tests/utils/logging.jl")
        include("tests/utils/misc.jl")
        include("tests/utils/print.jl")
        include("tests/utils/subset.jl")
        include("tests/utils/tables.jl")
        include("tests/utils/utils.jl")
    end

    @testset "Layouts" begin
        include("tests/layouts/layout.jl")
    end

    @testset "I/O" begin
        include("tests/io/read_data.jl")
        include("tests/io/eeglab_import.jl")
        include("tests/io/bids_import.jl")
        include("tests/io/bids_export.jl")
        include("tests/io/xdf_import.jl")
        include("tests/io/fieldtrip_import.jl")
    end

    @testset "Config & Pipelines" begin
        include("tests/config/config.jl")
        include("tests/pipelines/pipeline.jl")
        include("tests/pipelines/report.jl")
    end

    @testset "Analysis — Processing" begin
        include("tests/analysis/processing/artifact_detection/artifact_detection.jl")
        include("tests/analysis/processing/artifact_detection/autoreject.jl")
        include("tests/analysis/processing/signal/baseline.jl")
        include("tests/analysis/processing/signal/filter.jl")
        include("tests/analysis/processing/signal/cleanline.jl")
        include("tests/analysis/processing/signal/mirror.jl")
        include("tests/analysis/processing/signal/resample.jl")
        include("tests/analysis/processing/channels/channel_average.jl")
        include("tests/analysis/processing/channels/channel_delete.jl")
        include("tests/analysis/processing/channels/channel_difference.jl")
        include("tests/analysis/processing/channels/channel_metrics.jl")
        include("tests/analysis/processing/channels/channel_repair.jl")
        include("tests/analysis/processing/channels/channel_summary.jl")
        include("tests/analysis/processing/channels/rereference.jl")
        include("tests/analysis/processing/channels/csd.jl")
        include("tests/analysis/processing/epochs/triggers.jl")
        include("tests/analysis/processing/epochs/epochs.jl")
        include("tests/analysis/processing/epochs/apply.jl")
        include("tests/analysis/processing/conditions/condition_combine.jl")
        include("tests/analysis/processing/conditions/condition_difference.jl")
        include("tests/analysis/processing/conditions/condition_average.jl")
        include("tests/analysis/processing/conditions/grand_average.jl")
        include("tests/analysis/processing/ica/ica.jl")
    end

    @testset "Analysis — ERPs & RIDE" begin
        include("tests/analysis/erps/erp_measurements.jl")
        include("tests/analysis/erps/gfp.jl")
        include("tests/analysis/erps/jackknife_average.jl")
        include("tests/analysis/erps/lrp.jl")
        include("tests/analysis/erps/realign.jl")
        include("tests/analysis/ride/ride.jl")
    end

    @testset "Analysis — Statistics" begin
        include("tests/analysis/statistics/statistics_erp.jl")
        include("tests/analysis/statistics/statistics_tf.jl")
    end

    @testset "Analysis — Time-Frequency" begin
        include("tests/analysis/time_frequency/freq_spectrum.jl")
        include("tests/analysis/time_frequency/tf_utils.jl")
        include("tests/analysis/time_frequency/time_frequency.jl")
    end

    @testset "Analysis — Decoding & RSA" begin
        include("tests/analysis/decoding/decoding_preparation.jl")
        include("tests/analysis/decoding/decoding_erp.jl")
        include("tests/analysis/decoding/decoding_tf.jl")
        include("tests/analysis/decoding/decoding_statistics.jl")
        include("tests/analysis/rsa/rsa.jl")
    end

    @testset "Plotting" begin
        include("tests/plots/utils/layout_system.jl")
        include("tests/plots/utils/makie_recipes.jl")
        include("tests/plots/utils/plot_misc.jl")
        include("tests/plots/processing/plot_artifacts.jl")
        include("tests/plots/processing/plot_channel_summary.jl")
        include("tests/plots/processing/plot_correlation_heatmap.jl")
        include("tests/plots/processing/plot_databrowser.jl")
        include("tests/plots/processing/plot_filter.jl")
        include("tests/plots/processing/plot_ica.jl")
        include("tests/plots/processing/plot_layout.jl")
        include("tests/plots/processing/plot_triggers.jl")
        include("tests/plots/erps/plot_epochs.jl")
        include("tests/plots/erps/plot_epochs_rejection.jl")
        include("tests/plots/erps/plot_erp.jl")
        include("tests/plots/erps/plot_erp_image.jl")
        include("tests/plots/erps/plot_erp_stats.jl")
        include("tests/plots/erps/plot_gfp.jl")
        include("tests/plots/erps/plot_topography.jl")
        include("tests/plots/erps/plot_topography_3d.jl")
        include("tests/plots/erps/plot_topography_stats.jl")
        include("tests/plots/time_frequency/plot_power_spectrum.jl")
        include("tests/plots/time_frequency/plots_tf.jl")
        include("tests/plots/decoding/plot_decoding.jl")
        include("tests/plots/rsa/plot_rsa.jl")
    end

    @testset "Extensions" begin
        include("tests/ext/dimensional_data_ext.jl")
        include("tests/ext/mme_validation.jl")
    end

    @testset "Integration" begin
        include("tests/integration_and_plotting.jl")
    end

end

println("\nAll tests completed!")
