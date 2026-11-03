using Test
using Makie
using EegFun

@testset "plot_ica" begin
    # Needs a real setup with ICA components
    dat = EegFun.create_test_continuous_data(n = 1000, n_channels = 4)
    dat.layout = EegFun.create_test_layout(n_channels = 4)

    # Run mock ICA (since real ICA takes a while, we test what we can)
    ica_res = EegFun.run_ica(dat)

    @testset "plot_ica_component_activation" begin
        result = EegFun.plot_ica_component_activation(dat, ica_res; display_plot = false)
        @test result.fig isa Figure
    end

    @testset "plot_ica_component_spectrum" begin
        result = EegFun.plot_ica_component_spectrum(dat, ica_res; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "plot_topography ICA" begin
        result = EegFun.plot_topography(ica_res; display_plot = false)
        @test result.fig isa Figure
        @test result.axes isa Vector{Axis}
    end

    @testset "_get_artifact_categories and plot_artifact_components" begin
        artifacts = EegFun.ArtifactComponents(Dict(:vEOG => [1], :hEOG => [2]), [3], Int[], [4])
        cats = EegFun._get_artifact_categories(artifacts)
        @test cats isa Vector{Tuple{Symbol,String,Vector{Int}}}
        @test cats[1] == (:vEOG, "vEOG", [1])
        @test cats[2] == (:hEOG, "hEOG", [2])
        @test cats[3] == (:ecg, "Ecg", [3])
        @test cats[5] == (:channel_noise, "Channel Noise", [4])

        res_art = EegFun.plot_artifact_components(ica_res, artifacts; display_plot = false)
        @test res_art.fig isa Figure
        @test res_art.axes isa Vector{Axis}
    end

    @testset "plot_ecg_component_features does not mutate caller DataFrame" begin
        using DataFrames
        metrics_df = DataFrame(
            Component = [1, 2, 3, 4],
            mean_ibi_s = [0.8, 0.9, NaN, 1.0],
            std_ibi_s = [0.04, 0.11, NaN, 0.25],
            peak_ratio = [0.95, 0.80, 0.10, 0.40],
        )
        cols_before = copy(names(metrics_df))
        res_ecg = EegFun.plot_ecg_component_features(
            [1],
            metrics_df;
            zscore_threshold_color = :purple,
            zscore_threshold_style = :dot,
            display_plot = false,
        )
        @test names(metrics_df) == cols_before
        @test !hasproperty(metrics_df, :heart_rate_bpm)
        @test res_ecg.fig isa Figure
        @test res_ecg.axes isa Vector{Axis}
    end

    @testset "ica.mean centering in plot_ica_component_activation" begin
        dat_shifted = copy(dat)
        ica_shifted = deepcopy(ica_res)
        ch1 = ica_shifted.layout.data.label[1]
        dat_shifted.data[!, ch1] .+= 500.0
        ica_shifted.mean[1] += 500.0

        res_orig = EegFun.plot_ica_component_activation(dat, ica_res; component_selection = EegFun.components(1), display_plot = false)
        res_shift =
            EegFun.plot_ica_component_activation(dat_shifted, ica_shifted; component_selection = EegFun.components(1), display_plot = false)
        line_orig = first(filter(p -> p isa Makie.Lines, res_orig.axes[1].scene.plots))
        line_shift = first(filter(p -> p isa Makie.Lines, res_shift.axes[1].scene.plots))
        y_orig = [pt[2] for pt in line_orig[1][]]
        y_shift = [pt[2] for pt in line_shift[1][]]
        @test y_orig ≈ y_shift
    end
end
