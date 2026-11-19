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

    @testset "empty component_selection, >10 subset pagination, and ECG threshold line preservation" begin
        @test_throws EegFun.EegFunError EegFun.plot_ica_component_activation(
            dat,
            ica_res;
            component_selection = EegFun.components(Int[]),
            display_plot = false,
        )

        dat15 = EegFun.create_test_continuous_data(n = 300, n_channels = 15, seed = 9)
        ica15 = EegFun.run_ica(dat15; params = EegFun.IcaPrms(max_iter = 5))
        res_sub12 = EegFun.plot_ica_component_activation(dat15, ica15; component_selection = EegFun.components(1:12), display_plot = false)
        @test length(res_sub12.axes) == 10

        metrics_df = DataFrame(Component = [1, 2], mean_ibi_s = [0.8, 0.9], std_ibi_s = [0.03, 0.25], peak_ratio = [0.95, 0.40])
        res_ecg = EegFun.plot_ecg_component_features([1], metrics_df; min_peak_ratio = 0.72, max_ibi_std_s = 0.18, display_plot = false)
        lbls = filter(b -> b isa Label, res_ecg.fig.content)
        @test any(l -> occursin("0.72", l.text[]) && occursin("0.18", l.text[]), lbls)
    end

    @testset "use_global_scale shares contour levels across ICA topographies" begin
        ica_scaled = deepcopy(ica_res)
        ica_scaled.mixing[:, 1] .*= 100.0
        ica_scaled.mixing[:, 2] .*= 0.1

        res_ind = EegFun.plot_topography(
            ica_scaled;
            component_selection = EegFun.components([1, 2]),
            use_global_scale = false,
            gridscale = 15,
            display_plot = false,
        )
        res_glob = EegFun.plot_topography(
            ica_scaled;
            component_selection = EegFun.components([1, 2]),
            use_global_scale = true,
            gridscale = 15,
            display_plot = false,
        )
        cf_ind1 = only(filter(p -> p isa Makie.Contourf, res_ind.axes[1].scene.plots))
        cf_ind2 = only(filter(p -> p isa Makie.Contourf, res_ind.axes[2].scene.plots))
        cf_glob1 = only(filter(p -> p isa Makie.Contourf, res_glob.axes[1].scene.plots))
        cf_glob2 = only(filter(p -> p isa Makie.Contourf, res_glob.axes[2].scene.plots))
        @test collect(cf_ind1.levels[]) != collect(cf_ind2.levels[])
        @test collect(cf_glob1.levels[]) == collect(cf_glob2.levels[])

        artifacts_gs = EegFun.ArtifactComponents(Dict(:vEOG => [1], :hEOG => [2]), Int[], Int[], Int[])
        art_glob = EegFun.plot_artifact_components(ica_scaled, artifacts_gs; use_global_scale = true, gridscale = 15, display_plot = false)
        acf1 = only(filter(p -> p isa Makie.Contourf, art_glob.axes[1].scene.plots))
        acf2 = only(filter(p -> p isa Makie.Contourf, art_glob.axes[2].scene.plots))
        @test collect(acf1.levels[]) == collect(acf2.levels[])

        act_glob = EegFun.plot_ica_component_activation(
            dat,
            ica_scaled;
            component_selection = EegFun.components([1, 2]),
            use_global_scale = true,
            gridscale = 15,
            display_plot = false,
        )
        topo_axes = filter(a -> a isa Axis && a ∉ act_glob.axes, act_glob.fig.content)
        vcf1 = only(filter(p -> p isa Makie.Contourf, topo_axes[1].scene.plots))
        vcf2 = only(filter(p -> p isa Makie.Contourf, topo_axes[2].scene.plots))
        @test collect(vcf1.levels[]) == collect(vcf2.levels[])
    end
end
