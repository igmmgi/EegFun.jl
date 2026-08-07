using Test
using EegFun
using GLMakie
println("Running EegFun.jl Test Suite")
println("="^40)

@testset "Code quality (Aqua.jl)" begin
    using Aqua
    Aqua.test_all(
        EegFun;
        ambiguities = false,
        persistent_tasks = false,
        deps_compat = (check_extras = false, ignore = [:Dates, :LinearAlgebra, :Logging, :Printf, :Random, :SparseArrays, :TOML, :Test]),
    )
end

@testset "EegFun" begin

    include("tests/apply.jl")
    include("tests/artifact_detection.jl")
    include("tests/baseline.jl")
    include("tests/batch_utils.jl")
    include("tests/channel_average.jl")
    include("tests/channel_delete.jl")
    include("tests/channel_difference.jl")
    include("tests/channel_metrics.jl")
    include("tests/channel_repair.jl")
    include("tests/channel_summary.jl")
    include("tests/condition_combine.jl")
    include("tests/condition_average.jl")
    include("tests/condition_difference.jl")
    include("tests/config.jl")
    include("tests/pipeline.jl")
    include("tests/data.jl")
    include("tests/eeglab_import.jl")
    include("tests/bids_import.jl")
    include("tests/xdf_import.jl")
    include("tests/fieldtrip_import.jl")
    include("tests/read_data.jl")
    include("tests/epochs.jl")
    include("tests/erp_measurements.jl")
    include("tests/files.jl")
    include("tests/filter.jl")
    include("tests/gfp.jl")
    include("tests/grand_average.jl")
    include("tests/ica.jl")
    include("tests/jackknife_average.jl")
    include("tests/layout.jl")
    include("tests/layout_system.jl")
    include("tests/logging.jl")
    include("tests/lrp.jl")
    include("tests/mirror.jl")
    include("tests/misc.jl")
    include("tests/plot_channel_summary.jl")
    include("tests/plot_layout.jl")
    include("tests/plot_topography_3d.jl")
    include("tests/plot_misc.jl")
    include("tests/plot_triggers.jl")
    include("tests/plot_erp.jl")
    include("tests/plot_epochs.jl")
    include("tests/plot_gfp.jl")
    include("tests/plot_erp_image.jl")
    include("tests/plot_erp_stats.jl")
    include("tests/plot_artifacts.jl")
    include("tests/plot_correlation_heatmap.jl")
    include("tests/plot_filter.jl")
    include("tests/plot_ica.jl")
    include("tests/plot_decoding.jl")
    include("tests/plot_rsa.jl")
    include("tests/plot_power_spectrum.jl")
    include("tests/plot_databrowser.jl")
    include("tests/print.jl")
    include("tests/realign.jl")
    include("tests/rereference.jl")
    include("tests/resample.jl")
    include("tests/statistics_erp.jl")
    include("tests/statistics_tf.jl")
    include("tests/plots_tf.jl")
    include("tests/decoding_tf.jl")
    include("tests/decoding_erp.jl")
    include("tests/decoding_statistics.jl")
    include("tests/rsa.jl")
    include("tests/subset.jl")
    include("tests/time_frequency.jl")
    include("tests/tf_utils.jl")
    include("tests/triggers.jl")
    include("tests/types.jl")
    include("tests/utils.jl")
    include("tests/bids_export.jl")

end

println("\nAll tests completed!")
