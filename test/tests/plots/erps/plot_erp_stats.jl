# Tests for ERP condition-comparison and significant-interval statistical plotting.

using Test
using CairoMakie
using EegFun

@testset "plot_erp_stats" begin
    dat1 = EegFun.create_test_erp_data(participant = 1, n_channels = 3)
    dat2 = EegFun.create_test_erp_data(participant = 2, n_channels = 3)
    dat3 = EegFun.create_test_erp_data(participant = 1, condition = 2, n_channels = 3)
    dat4 = EegFun.create_test_erp_data(participant = 2, condition = 2, n_channels = 3)

    all_dat = [dat1, dat2, dat3, dat4]

    stat_data = EegFun.prepare_stats(all_dat; design = :paired)
    res = EegFun.analytic_test(stat_data)

    @testset "basic functionality" begin
        result = EegFun.plot_erp_stats(res; display_plot = false)
        @test result.fig isa Figure
        @test length(result.axes) > 0
    end

    @testset "channel selection" begin
        result = EegFun.plot_erp_stats(res; channel_selection = EegFun.channels([1]), display_plot = false)
        @test result.fig isa Figure

        # Symbol channel selection and String channel_plot_order
        res_sym = EegFun.plot_erp_stats(res; channel_selection = :Ch1, channel_plot_order = ["Ch1"], display_plot = false)
        @test res_sym.fig isa Figure

        # String channel selection and String channel_plot_order
        res_str = EegFun.plot_erp_stats(res; channel_selection = ["Ch1", "Ch2"], channel_plot_order = ["Ch2", "Ch1"], display_plot = false)
        @test res_str.fig isa Figure
    end

    @testset "PermutationResult critical_t shapes and significance_position" begin
        for thresh_method in (:parametric, :nonparametric_common, :nonparametric_individual)
            perm_res = EegFun.permutation_test(stat_data; n_permutations = 10, threshold_method = thresh_method)
            out_perm =
                EegFun.plot_erp_stats(perm_res; plot_critical_t = true, significance_position = :top, legend = false, display_plot = false)
            @test out_perm.fig isa Figure
            @test out_perm.axes isa Vector{Axis}
            @test length(out_perm.axes) > 0
        end
    end

    @testset "plot_erp_measurements and plot_erp_measurement_gui" begin
        # Multi-participant vector returns one figure NamedTuple per participant
        res_meas_multi = EegFun.plot_erp_measurements(all_dat, "mean_amplitude"; analysis_interval = (0.0, 0.2), display_plot = false)
        @test res_meas_multi isa Vector{NamedTuple}
        @test length(res_meas_multi) == 2
        @test res_meas_multi[1].fig isa Figure
        @test res_meas_multi[1].axes isa Vector{Axis}

        # Single-participant condition overlay aligns overlays with channel_plot_order
        p1_c1 = copy(dat1)
        p1_c1.data[!, :Ch1] .= 10.0
        p1_c1.data[!, :Ch2] .= -20.0
        res_ord = EegFun.plot_erp_measurements(
            [p1_c1],
            "mean_amplitude";
            layout = :grid,
            channel_selection = EegFun.channels([:Ch1, :Ch2]),
            channel_plot_order = [:Ch2, :Ch1],
            analysis_interval = (0.0, 0.2),
            display_plot = false,
        )
        @test res_ord.axes[1].title[] == "Ch2"
        @test res_ord.axes[2].title[] == "Ch1"
        txt_ax1 = only(filter(p -> p isa Makie.Text, res_ord.axes[1].scene.plots))
        txt_ax2 = only(filter(p -> p isa Makie.Text, res_ord.axes[2].scene.plots))
        @test occursin("-20.00", string(txt_ax1.text[]))
        @test occursin("10.00", string(txt_ax2.text[]))

        gui_res = EegFun.plot_erp_measurement_gui(dat1; analysis_type = "max_peak_latency", display_plot = false)
        @test gui_res.fig isa Figure
        @test gui_res.axes isa Vector{Axis}
        @test haskey(gui_res, :screen)
    end
end
